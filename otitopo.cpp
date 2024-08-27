#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include "pzlog.h"
#include "pzgmesh.h"
#include "TPZGenGrid2D.h"
#include "TPZVTKGeoMesh.h"
#include "pzcmesh.h"
#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <TPZNullMaterial.h>
#include <pzbuildmultiphysicsmesh.h>
#include <pzskylstrmatrix.h>
#include <pzskylstrmatrix.h>
#include <pzstepsolver.h>
#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <TPZSimpleTimer.h>
#include "pzvisualmatrix.h"
#include "TPZSYSMPMatrix.h"
#include "TPZVTKGenerator.h"
//#include <Elasticity/TPZElasticity3D.h>
#include <Elasticity/TPZElasticity2DOtiTopo.h>
#include "TPZAnalyticSolution.h"
#include "TPZGeoMeshTools.h"
#include <TPZGmshReader.h>
#include "tpzchangeel.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "pzcompelwithmem.h"
#include "TPZCompElH1.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "TPZElementMatrixT.h"
#include "filterstruct.h"
#include "TPZRefPatternTools.h"

using namespace std;

enum EMatid {ENone,EDomain,EDispX,EDispY,EDispXY,EPtDispX,EPtDispY,EPtDispXY,EForceX,EForceY,EForceXY,EPtForceX,EPtForceY,EPtForceXY};
const int global_nthread = 8;

TPZGeoMesh* CreateGMesh(int ndivx, int ndivy);
TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);
void CreateBCs(TPZGeoMesh* gmesh);
TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord, TElasticity2DAnalytic *elas);
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh, const std::string filename);

void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc);

const bool LoadMemoryIntoElementSolution(TPZLinearAnalysis& an, TPZCompMesh *&cmesh, bool isUpdate, TPZVec<FilterStruct> &filterVec, TPZVec<std::set<int>> &neighVec);

void GetSolVec(TPZInterpolationSpace* intel, TPZFMatrix<STATE>& u);

REAL calcVol(TPZCompMesh *cmesh);
void CreateFilterVec(TPZGeoMesh* gmesh, TPZVec<FilterStruct> &filterVec, REAL rmin);
const bool RefineElements(const REAL rhovartol, TPZGeoMesh* gmesh, TPZVec<std::set<int>> &neighVec, TPZVec<bool>& isCompVec, TPZVec<REAL> &rhovecnew);
void CheckRefinedNeighbors(TPZGeoMesh* gmesh, TPZVec<REAL> &rhovecnew);
void UpdateSonsMemory(TPZCompMesh* cmesh, TPZVec<int64_t>& subindexes, REAL rhofather);
void InitializeElemSolOfRefElements(TPZCompMesh* cmesh);
void CreateNeighVec(TPZVec<std::set<int>>& neighVec, TPZGeoMesh* gmesh);

REAL gVolInit = 0.;

int main() {
    
    std::cout << "--------- Starting simulation ---------" << std::endl;
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    
    const int pord = 1;
    const int niterations = 140;
    
    int ndivx = 25, ndivy = 50;
    
    const REAL filterRadius = 1.5;

    bool isReadFromGmsh = true;
    TPZGeoMesh* gmesh = nullptr;
    if (isReadFromGmsh) {
//        gmesh = ReadMeshFromGmsh("beam_coarse.msh");
        // gmesh = ReadMeshFromGmsh("beam.msh");
        gmesh = ReadMeshFromGmsh("beam80_40.msh");        
    }
    else{
        gmesh = CreateGMesh(ndivx,ndivy);
    }

    int firstindex = 0;
    for(int i = 0 ; i < gmesh->NElements() ; i++){
        TPZGeoEl* gel = gmesh->ElementVec()[i];
        if(!gel) continue;
        if (gel->Dimension() == gmesh->Dimension()){
            firstindex = gel->Index();
            break;
        }
    }

    // Create a geoel for the first node in the mesh
    TPZGeoElBC gelbc(gmesh->Element(firstindex),3,50);

    std::set<int> matidref;
    matidref.insert(50);
    // gRefDBase.InitializeRefPatterns(gmesh->Dimension());
    // TPZRefPatternTools::RefineDirectional(gmesh, matidref);

    
    TPZVec<TPZGeoEl*> sons;
    // gmesh->Element(firstindex)->Divide(sons);    

    // FilterStruct filter1(firstindex);
    // filter1.ComputeNeighIndexHf(gmesh, 0.6);
    // filter1.Print(std::cout); 

    std::ofstream out("gmesh.vtk");     
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

    TElasticity2DAnalytic *elas = new TElasticity2DAnalytic;
    elas->gE = 1.;//206.8150271873455;
    elas->gPoisson = 0.3;
    elas->fProblemType = TElasticity2DAnalytic::EStretchx;
    TPZCompMesh* cmeshH1 = CreateH1CMesh(gmesh,pord,elas);

    if (0)
    {
        // get second connect of first son in sons
        TPZGeoEl* gel = sons[0];
        TPZCompEl* cel = gel->Reference();
        TPZInterpolationSpace* intel = dynamic_cast<TPZInterpolationSpace*>(cel);
        if (!intel) DebugStop();
        TPZConnect &c = intel->Connect(2);
        int64_t cind = intel->ConnectIndex(2);
        TPZFMatrix<STATE> depmat(1,1,0.25);        
        int indextosave = -1;
        for (int is = 0; is < sons.size(); is++) {
            TPZGeoEl* gelnext = sons[is];
            TPZCompEl* celnext = gelnext->Reference();
            TPZCompElWithMem <TPZCompElH1<pzshape::TPZShapeQuad> > *celmem = dynamic_cast<TPZCompElWithMem <TPZCompElH1<pzshape::TPZShapeQuad> >*>(cel);
            TPZInterpolationSpace* intelnext = dynamic_cast<TPZInterpolationSpace*>(celnext);
            TPZMatWithMem<TPZOtiTopoDensity>* mat = dynamic_cast<TPZMatWithMem<TPZOtiTopoDensity>*>(celmem->Material());
            TPZVec<int64_t> indices;
            celmem->GetMemoryIndices(indices);
            if(is == 0){
                indextosave = indices[0];
            }
            else{
                indices[0] = indextosave;
                celmem->SetMemoryIndices(indices);
            }
            // TPZOtiTopoDensity &densstruct = mat->MemItem(indices[0]);
            // const STATE dens = densstruct.fDen;
            if (!intelnext) DebugStop();
            TPZConnect &c2 = intelnext->Connect(is);
            int64_t c2ind = intelnext->ConnectIndex(is);
            c.AddDependency(cind, c2ind, depmat, 0, 0, 1, 1);
        }
        cmeshH1->CleanUpUnconnectedNodes();
        c.Print(*cmeshH1);
    }

    // Compute FilterStruct for all the elements in the mesh
    TPZVec<FilterStruct> filterVec;
    CreateFilterVec(gmesh, filterVec, filterRadius);      

    // Create vector of sets with size equal to the number of geometric elements
    TPZVec<std::set<int>> neighVec;
    CreateNeighVec(neighVec, gmesh);

    TPZLinearAnalysis an(cmeshH1);
    LoadMemoryIntoElementSolution(an,cmeshH1,false,filterVec,neighVec);        
    gVolInit = calcVol(cmeshH1);
    const string filenameinit = "initpostproc";
    PrintResults(an,cmeshH1,filenameinit);
        
    for (int i = 0; i < niterations; i++) {
        std::cout << "***************** Optimization Step " << i << " *****************" << std::endl;
        SolveProblemDirect(an,cmeshH1);
        const bool isNewMesh = LoadMemoryIntoElementSolution(an,cmeshH1,true,filterVec,neighVec);
        cout << "Number of elements in cmesh after refine and outside: " << cmeshH1->NElements() << endl;
        if(isNewMesh){
            CreateFilterVec(gmesh, filterVec, filterRadius);
            CreateNeighVec(neighVec, gmesh);
        }
        // std::cout << "--------- PostProcess ---------" << std::endl;
        // PrintResults(an,cmeshH1);
    }
    
    delete cmeshH1;
    delete gmesh;
        
    std::cout << "--------- Simulation finished ---------" << std::endl;
}

TPZGeoMesh* CreateGMesh(int ndivx, int ndivy) {
    TPZGeoMesh* gmesh = new TPZGeoMesh;
    
    MMeshType meshType = MMeshType::EQuadrilateral;
    int dim = 2;
    TPZManVector<REAL,3> minX = {0,0,0};
    TPZManVector<REAL,3> maxX = {100,200,0};
    int nMats = 2*dim+1;
    
    constexpr bool createBoundEls{true};
    TPZVec<int> matIds(nMats,ENone);
    matIds[0] = EDomain;
    
    TPZManVector<int,2> ndivvec = {ndivx,ndivy};
    gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,matIds, ndivvec, meshType,createBoundEls);
    
    TPZManVector<REAL,2> xfixed1 = {0.,0.,0.}, xfixed2 = {0.,200.,0.}, xforce = {100.,100.,0.};
    SetPointBC(gmesh, xfixed1, EPtDispXY);
    SetPointBC(gmesh, xfixed2, EPtDispXY);
    SetPointBC(gmesh, xforce, EPtForceY);

    return gmesh;
}


TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord, TElasticity2DAnalytic *elas) {

    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    const int dim = gmesh->Dimension();
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pord);
    cmesh->SetAllCreateFunctionsContinuousWithMem();
    
    const STATE E = elas->gE, nu = elas->gPoisson;
    TPZManVector<STATE> force = {0,0,0};
    TPZElasticity2DOtiTopo *mat = new TPZElasticity2DOtiTopo(EDomain, E, nu, 0., 0., true);
    mat->SetExactSol(elas->ExactSolution(), 2);
//    mat->SetForcingFunction(elas->ForceFunc(), 4);
    cmesh->InsertMaterialObject(mat);
    
    TPZFMatrix<STATE> val1(2,2,0.);
    TPZManVector<STATE> val2(2,0.);
    
    const int diri = 0, neu = 1, mixed = 2, normaltrac = 4;
//    auto* BCCond0 = mat->CreateBC(mat, EBC, diri, val1, val2);
    
//    BCCond0->SetForcingFunctionBC(elas->ExactSolution(), 4);
//    cmesh->InsertMaterialObject(BCCond0);
    
    // Pointer bcs
    val1(1,1) = mat->BigNumber();
    auto* BCCondFixed1 = mat->CreateBC(mat, EPtDispXY, mixed, val1, val2);
    cmesh->InsertMaterialObject(BCCondFixed1);
    
    val1.Zero();
    val1(0,0) = mat->BigNumber();
    auto* BCCondSym = mat->CreateBC(mat, EDispX, mixed, val1, val2);
    cmesh->InsertMaterialObject(BCCondSym);
    
    val1.Zero();
    val2[1] = -1.0;
    auto* BCCondPoint = mat->CreateBC(mat, EPtForceY, neu, val1, val2);
    cmesh->InsertMaterialObject(BCCondPoint);
    
    cmesh->AutoBuild();
    //este método é chamado para construir automaticamente a malha computacional com base nas configurações e objetos de material e condição de contorno definidos anteriormente.
    
    return cmesh;
    //a função retorna o ponteiro para a malha computacional cmesh
}

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
//Declara a função SolveProblemDirect do tipo void (logo não retornará valor).
{

//    TPZSkylineStructMatrix<STATE> matskl(cmesh);
    TPZSSpStructMatrix<STATE> matskl(cmesh);
    matskl.SetNumThreads(global_nthread);
    an.SetStructuralMatrix(matskl);
    
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);//ELU //ECholesky // ELDLt
    an.SetSolver(step);
    
    //assembles the system
    std::cout << "--------- Assemble ---------" << std::endl;
    TPZSimpleTimer time_ass;
    an.Assemble();
    std::cout << "Total time = " << time_ass.ReturnTimeDouble()/1000. << " s" << std::endl;

    std::cout << "--------- Solve ---------" << std::endl;
    TPZSimpleTimer time_sol;
    an.Solve();
    std::cout << "Total time = " << time_sol.ReturnTimeDouble()/1000. << " s" << std::endl;
    
    return;
}

void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh, const std::string filename)
//Define a função do tipo void chamada PrintResults, que recebe como parâmetros TPZLinearAnalysis &an e  TPZCompMesh *cmesh
{
 
    std::cout << "--------- Post Process ---------" << std::endl;
    //printa para o usuário "--------- Post Process ---------" indicando o início da fase de pós-processamento. 
    TPZSimpleTimer postProc("Post processing time");
    //declara uma variável chamada postProc, do tipo TPZSimpleTimer, chamando um construtor com uma string como argumento, igual a "Post processing time".
    //inicializa um temporizador chamado postProc que será usado para medir o tempo gasto no pós-processamento.
    const std::string plotfile = "postprocess";
    //define o nome base do arquivo de saída para o pós-processamento. O nome base é "postprocess".
    constexpr int vtkRes{0};
    //define a variável do tipo inteiro denominada vtkRes, do tipo constexpr, que significa que é uma expressão constante, ou seja,  vtkRes é um valor constante e não pode ser alterado. Ainda, {0} indica o valor associado a essa constante, e portanto não será alterado, com valor determinado na hora de compilação.
    //define a resolução para o formato de arquivo VTK. Neste caso, a resolução é definida como 0, o que geralmente significa que a resolução será automática.
    TPZVec<std::string> fields = {
        "Displacement",
//        "Stress",
//        "Strain",
        "topodensity"
    };
//    TPZVec<std::string> fieldsvec = {};
//    an.DefineGraphMesh(cmesh->Dimension(), fields, fieldsvec, "postproc.vtk");
//    an.PostProcess(0, cmesh->Dimension());
    
    //nesse conjunto de linhas de código, temos que TPZVec é uma estrutura do tipo vetor que contém como argumento uma variável chamda "fields" que é uma lista de strings, que, pelo que se chamam, são relacionadas ao deslocamento, deformação e tensão.
    //cria um vetor de strings chamado fields que contém os nomes dos campos que serão pós-processados. Neste caso, os campos incluem "Displacement" (deslocamento), "Stress" (tensão) e "Strain" (deformação). Esses campos representam propriedades do problema que desejamos visualizar após a simulação.
    static auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
    //essa linha de código declara uma variável chamada vtk do tipo auto, o que significa que o compilador irá deduzir o tipo que ela terá a depender do que ela é igual. No caso, ela é igual a função TPZVTKGenerator, de parâmetros cmesh, fields, plotfile, vtkRes.
    //cria um objeto vtk da classe TPZVTKGenerator, que é usado para gerar arquivos VTK a partir dos dados da malha computacional cmesh. Os argumentos passados para o construtor incluem a malha computacional, os campos a serem pós-processados, o nome base do arquivo de saída (plotfile) e a resolução VTK (vtkRes).
    vtk.SetNThreads(0);
    //define o número de threads a serem usadas durante o pós-processamento. A variável global_nthread provavelmente contém o número desejado de threads.
    vtk.Do();
    //inicia o processo de geração dos arquivos VTK. Esta função gera arquivos de saída contendo informações sobre os campos especificados na malha computacional.
    std::cout << "Total time = " << postProc.ReturnTimeDouble()/1000. << " s" << std::endl;
    //imprime o tempo gasto no pós-processamento, convertido para segundos.
    
    return;
    //a função é concluída e retorna.
}

const bool LoadMemoryIntoElementSolution(TPZLinearAnalysis& an, TPZCompMesh *&cmesh, bool isUpdate, TPZVec<FilterStruct> &filterVec, TPZVec<std::set<int>> &neighVec) {
    
    TPZSimpleTimer timer("LoadMemoryIntoElementSolution");
    REAL volAtStep = calcVol(cmesh);
    const REAL volFrac = volAtStep/gVolInit;
    const STATE mindens = 1.e-3;
    REAL c = 0.;
    REAL energy = 0.;
    TPZVec<REAL> dcvec(cmesh->NElements(),0.), dcvecfilter(cmesh->NElements(),0.), rhovec(cmesh->NElements(),0.), 
        rhovecnew(cmesh->NElements(),0.), elvolvec(cmesh->NElements(),0.);
    TPZVec<bool> isCompVec(cmesh->NElements(),false);
    
    TPZFMatrix<STATE> &elementSol = cmesh->ElementSolution();
    const int64_t nel = cmesh->NElements();
    elementSol.Resize(nel, 1);
    for(int64_t i = 0 ; i < nel ; i++) {
        TPZCompEl* cel = cmesh->Element(i);        
        if(!cel) continue;
        TPZGeoEl* gel = cel->Reference();
        if(gel->HasSubElement()) continue;
        TPZCompElWithMem <TPZCompElH1<pzshape::TPZShapeQuad> > *celmem = dynamic_cast<TPZCompElWithMem <TPZCompElH1<pzshape::TPZShapeQuad> >*>(cel);
        TPZCompElWithMem <TPZCompElH1<pzshape::TPZShapeTriang> > *celmemtri = dynamic_cast<TPZCompElWithMem <TPZCompElH1<pzshape::TPZShapeTriang> >*>(cel);
        if(!celmem && !celmemtri) continue;
        TPZVec<int64_t> indices;
        if(celmem) celmem->GetMemoryIndices(indices);
        else celmemtri->GetMemoryIndices(indices);
            
        for(int j = 1 ; j < indices.size() ; j++) {
            if (indices[0] != indices[j]) {
                DebugStop(); // assuming same memory for the whole element!
            }
        }
        TPZMatWithMem<TPZOtiTopoDensity>* mat = nullptr;
        if(celmem) mat = dynamic_cast<TPZMatWithMem<TPZOtiTopoDensity>*>(celmem->Material());
        else mat = dynamic_cast<TPZMatWithMem<TPZOtiTopoDensity>*>(celmemtri->Material());
        
        TPZElasticity2D* matElas = nullptr;
        if(celmem) matElas = dynamic_cast<TPZElasticity2D*>(celmem->Material());
        else matElas = dynamic_cast<TPZElasticity2D*>(celmemtri->Material());

        if(!mat) DebugStop();
        if(!matElas) DebugStop();
        isCompVec[cel->Index()] = true;
        TPZOtiTopoDensity &densstruct = mat->MemItem(indices[0]);
        const STATE dens = densstruct.fDen;
        const int64_t index = cel->Index();
        if (isUpdate) {
            TPZElementMatrixT<STATE> ek, ef;
            cel->CalcStiff(ek, ef);
            TPZInterpolationSpace* intel = dynamic_cast<TPZInterpolationSpace*>(cel);
            if (!intel) DebugStop();
            if( i == 3242)
                std::cout << "aqui" << endl;            
            // const int64_t elneq = intel->NEquations();
            const int64_t elneq = intel->NShapeF() * matElas->NStateVariables();
            TPZFMatrix<STATE> u(elneq,1,0.), ku(elneq,1,0.);

            GetSolVec(intel,u);
            ek.Matrix().Multiply(u, ku);
            STATE E = 0.;
            for (int idof = 0; idof < elneq; idof++) {
                E += u(idof,0) * ku(idof,0);
            }            
            
            const REAL p = 3.;
//            c += pow(dens, p) * E;
            c += E;
            energy += E;
//            dcvec[cel->Index()] = - p * pow(dens, p-1) * E;
            dcvec[index] = - p * E / dens;
            rhovec[index] = densstruct.fDen;
            elvolvec[index] = gel->Volume();
            
//            elementSol(index,0) = updatedValue;
//            densstruct.fDen = updatedValue;
        }
        else{
            const STATE initdens = 1.;
            elementSol(index,0) = initdens;
            densstruct.fDen = initdens;
        }
    }

    // Apply filter
    const bool isUseFilter = true;
    if(isUseFilter){
        for (int i = 0; i < dcvec.size(); i++) {
            if (!isCompVec[i]) {
                continue;
            }
            TPZCompEl* cel = cmesh->Element(i);
            if(!cel) DebugStop();
            
            if (filterVec[i].findex == -1) DebugStop();        

            if (i == 250){
                std::cout << "aqui" << endl;
            }
            const REAL dcdxe = filterVec[i].ComputeFiltereddcdxe(dcvec, rhovec, cmesh);
            // if(dcdxe > 0)
            //     DebugStop(); // filter cannot change sign of dc
            dcvecfilter[i] = dcdxe;
        }
    }
    else{
        dcvecfilter = dcvec;
    }
    
    
    // Loop for Optimality criteria updated
    REAL l1 = 0., l2 = 100000., move = 0.2;
    const REAL targetVolFrac = 0.3;
    while (l2-l1 > 1.e-4) {
        REAL lmid = 0.5*(l1+l2);
        for (int i = 0; i < rhovec.size(); i++) {
            const REAL dens = rhovec[i];
            const REAL dc = dcvecfilter[i];
            const REAL first = std::min(dens+move, dens*sqrt(-dc/lmid));
            const REAL second = std::min(first,1.);
            const REAL third = std::max(dens-move,second);
            const REAL fourth = std::max(mindens,third);
            if (fourth > 1)
                DebugStop();
            if (fourth < mindens)
                DebugStop();

            rhovecnew[i] = fourth;
        }
        REAL newvol = 0.;
        for (int i = 0; i < rhovecnew.size(); i++) {
            if(isCompVec[i]){
                newvol += rhovecnew[i] * elvolvec[i];
            }
        }
        if(newvol - targetVolFrac*gVolInit > 0){
            l1 = lmid;
        }
        else{
            l2 = lmid;
        }
    }
    
    // Update value of rho inside each element
    for(int64_t i = 0 ; i < nel ; i++) {
        TPZCompEl* cel = cmesh->Element(i);
        if(!cel) continue;
        TPZCompElWithMem <TPZCompElH1<pzshape::TPZShapeQuad> > *celmem = dynamic_cast<TPZCompElWithMem <TPZCompElH1<pzshape::TPZShapeQuad> >*>(cel);
        TPZCompElWithMem <TPZCompElH1<pzshape::TPZShapeTriang> > *celmemtri = dynamic_cast<TPZCompElWithMem <TPZCompElH1<pzshape::TPZShapeTriang> >*>(cel);
        if(!celmem && !celmemtri) continue;
        TPZVec<int64_t> indices;

        if(celmem) celmem->GetMemoryIndices(indices);
        else celmemtri->GetMemoryIndices(indices);

        for(int j = 1 ; j < indices.size() ; j++) {
            if (indices[0] != indices[j]) {
                DebugStop(); // assuming same memory for the whole element!
            }
        }
        TPZMatWithMem<TPZOtiTopoDensity>* mat = nullptr;
        if(celmem) mat = dynamic_cast<TPZMatWithMem<TPZOtiTopoDensity>*>(celmem->Material());
        else mat = dynamic_cast<TPZMatWithMem<TPZOtiTopoDensity>*>(celmemtri->Material());
        if(!mat) DebugStop();
        isCompVec[cel->Index()] = true;
        TPZOtiTopoDensity &densstruct = mat->MemItem(indices[0]);
        const STATE dens = densstruct.fDen;
        const int64_t index = cel->Index();
        if (isUpdate) {
            REAL newdens = rhovecnew[cel->Index()];
            elementSol(index,0) = newdens;
            densstruct.fDen = newdens;
        }
    }

    // Print results to Paraview
    std::cout << "--------- PostProcess ---------" << std::endl;
    const string filename = "postproc";
    PrintResults(an,cmesh,filename);


    // Refine elements that have neighbours with variation of rho bigger than tolerance
    const REAL rhovartol = move/2.;
    TPZGeoMesh* gmesh = cmesh->Reference();
    bool isNewMesh = false;
    // isNewMesh = RefineElements(rhovartol, gmesh, neighVec, isCompVec, rhovecnew);    
    
    // Print results to Paraview
    if(isNewMesh){
        std::cout << "--------- PostProcess ---------" << std::endl;
        const string filename = "postproc_afterrefine";
        PrintResults(an,cmesh,filename);
    }

    std::cout << "c = " << c << std::endl;
    std::cout << "energy = " << energy << std::endl;
    std::cout << "\ntotal time for load element = " << timer.ReturnTimeDouble()/1000. << " seconds" << std::endl;

    return isNewMesh;
}


void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc) {
    // look for an element/corner node whose distance is close to start
    TPZGeoNode *gn1 = gr->FindNode(x);
    int64_t iel;
    int64_t nelem = gr->ElementVec().NElements();
    TPZGeoEl *gel;
    for (iel = 0; iel < nelem; iel++) {
        gel = gr->ElementVec()[iel];
        if (!gel) continue;
        int nc = gel->NCornerNodes();
        int c;
        for (c = 0; c < nc; c++) {
            TPZGeoNode *gn = gel->NodePtr(c);
            if (gn == gn1) {
                break;
            }
        }
        if (c < nc) {
            TPZGeoElBC(gel, c, bc);
            return;
        }
    }
}

void GetSolVec(TPZInterpolationSpace* intel, TPZFMatrix<STATE>& u) {
    const int nstate = intel->Material()->NStateVariables();
    const int ncon = intel->NConnects();
    TPZBlock &block = intel->Mesh()->Block();
    TPZFMatrix<STATE> &MeshSol = intel->Mesh()->Solution();
    const int64_t numbersol = MeshSol.Cols();
    if (numbersol != 1) DebugStop(); // I did not think about this case yet, but it can be done
    
    int64_t iv = 0;
    for(int in=0; in<ncon; in++) {
        TPZConnect *df = &intel->Connect(in);
        const int64_t dfseq = df->SequenceNumber();
        const int dfvar = block.Size(dfseq);
        const int64_t pos = block.Position(dfseq);
        for(int jn=0; jn<dfvar; jn++) {
            u(iv++,0) = MeshSol(pos+jn,0);
        }
    }
}

REAL calcVol(TPZCompMesh *cmesh) {
    REAL vol = 0.;
    TPZGeoMesh* gmesh = cmesh->Reference();
    const int64_t nel = cmesh->NElements();
    for(int64_t i = 0 ; i < nel ; i++) {
        TPZCompEl* cel = cmesh->Element(i);
        if(!cel) continue;
        TPZCompElWithMem <TPZCompElH1<pzshape::TPZShapeQuad> > *celmem = dynamic_cast<TPZCompElWithMem <TPZCompElH1<pzshape::TPZShapeQuad> >*>(cel);
        TPZCompElWithMem <TPZCompElH1<pzshape::TPZShapeTriang> > *celmemtri = dynamic_cast<TPZCompElWithMem <TPZCompElH1<pzshape::TPZShapeTriang> >*>(cel);        
        // if(!celmem && cel->Reference()->Dimension() == gmesh->Dimension()) DebugStop();
        if(!celmem && !celmemtri) continue;
        TPZVec<int64_t> indices;
        if (celmem) celmem->GetMemoryIndices(indices);
        else celmemtri->GetMemoryIndices(indices);
        for(int j = 1 ; j < indices.size() ; j++) {
            if (indices[0] != indices[j]) {
                DebugStop(); // assuming same memory for the whole element!
            }
        }
        TPZMatWithMem<TPZOtiTopoDensity>* mat = nullptr;
        if(celmem) mat = dynamic_cast<TPZMatWithMem<TPZOtiTopoDensity>*>(celmem->Material());
        else mat = dynamic_cast<TPZMatWithMem<TPZOtiTopoDensity>*>(celmemtri->Material());
        if(!mat) DebugStop();
        TPZOtiTopoDensity &densstruct = mat->MemItem(indices[0]);
        const STATE dens = densstruct.fDen;
        const int64_t index = cel->Index();
        REAL elvol = dens * cel->Reference()->Volume();
        vol += elvol;
    }
    return vol;
}

TPZGeoMesh*
ReadMeshFromGmsh(std::string file_name)
{
    //read mesh from gmsh
    TPZGeoMesh *gmesh;
    gmesh = new TPZGeoMesh();
    {
        TPZGmshReader reader;
        TPZManVector<std::map<std::string,int>,4> stringtoint(13);
        stringtoint[2]["dom"] = EDomain;
        
        stringtoint[1]["dispx"] = EDispX;
        stringtoint[1]["dispy"] = EDispY;
        stringtoint[1]["dispxy"] = EDispXY;
        stringtoint[0]["ptdispx"] = EPtDispX;
        stringtoint[0]["ptdispy"] = EPtDispY;
        stringtoint[0]["ptdispxy"] = EPtDispXY;

        stringtoint[1]["forcex"] = EForceX;
        stringtoint[1]["forcey"] = EForceY;
        stringtoint[1]["forcexy"] = EForceXY;
        stringtoint[0]["ptforcex"] = EPtForceX;
        stringtoint[0]["ptforcey"] = EPtForceY;
        stringtoint[0]["ptforcexy"] = EPtForceXY;
        reader.SetDimNamePhysical(stringtoint);
        reader.GeometricGmshMesh(file_name,gmesh);
    }

    return gmesh;
}

void CreateFilterVec(TPZGeoMesh* gmesh, TPZVec<FilterStruct> &filterVec, REAL rmin){    
    TPZCompMesh* cmesh = gmesh->Reference();
    if(!cmesh) DebugStop();
    const int nel = cmesh->NElements();
    filterVec.Resize(nel);    
    for(int i = 0 ; i < nel ; i++){
        TPZCompEl* cel = cmesh->ElementVec()[i];
        if(!cel) continue;
        if(cel->Reference()->HasSubElement()) continue;
        if (cel->Dimension() == gmesh->Dimension()){
            filterVec[i].findex = i;
            filterVec[i].ComputeNeighIndexHf(gmesh, rmin); 
        }
    }
}

const bool RefineElements(const REAL rhovartol, TPZGeoMesh* gmesh, TPZVec<std::set<int>> &neighVec, TPZVec<bool>& isCompVec, TPZVec<REAL> &rhovecnew){
    TPZCompMesh* cmesh = gmesh->Reference();
    // Loop over elements and insert on set in case of variation of rho between neighbors in NeighVec is bigger than tolerance
    std::set<int> elset;
    for (int i = 0; i < neighVec.size(); i++) {
        TPZCompEl* cel = cmesh->Element(i);        
        if(!cel) continue;
        if (!isCompVec[i]) {
            continue;
        }
        REAL thisdens = rhovecnew[i];
        for (auto it = neighVec[i].begin(); it != neighVec[i].end(); it++) {
            const int64_t celneighindex = *it;
            TPZCompEl* celneigh = cmesh->Element(celneighindex);
            if(!celneigh) continue;
            if (!isCompVec[celneighindex]) {
                DebugStop();
            }
            REAL neighdens = rhovecnew[celneighindex];
            if (fabs(thisdens - neighdens) > rhovartol) {
                elset.insert(i);
                break;
            }
        }
    }

    // Uniform refine all elements with index in set elset
    bool isNewMesh = false;
    std::cout << "Number of elements in gmesh: " << gmesh->NElements() << std::endl;    
    cout << "Number of elements in cmesh before refine: " << cmesh->NElements() << endl;
    const bool interpolatesol = false;
    for (auto it = elset.begin(); it != elset.end(); it++) {
        TPZCompEl* cel = cmesh->Element(*it);
        if(!cel) DebugStop();
        TPZGeoEl* gel = cel->Reference();
        if(!gel) DebugStop();
        if (gel->Dimension() != gmesh->Dimension()){
            DebugStop();
        }
        TPZVec<TPZGeoEl*> sons;
        TPZVec<int64_t> subindexes;        
        cel->Divide(cel->Index(),subindexes,interpolatesol);
        cout << "Number of elements in cmesh here: " << cmesh->NElements() << endl;
        // update the memory of the sons
        UpdateSonsMemory(cmesh,subindexes,rhovecnew[*it]);
        // for (int64_t ison = 0; ison < subindexes.size(); ison++) {
        //     TPZCompEl* celson = cmesh->Element(subindexes[ison]);
        //     TPZCompElWithMem<TPZCompElH1<pzshape::TPZShapeQuad>> *celmem = dynamic_cast<TPZCompElWithMem<TPZCompElH1<pzshape::TPZShapeQuad>>*>(celson);
        //     if (!celmem) DebugStop();
        //     TPZMatWithMem<TPZOtiTopoDensity>* mat = dynamic_cast<TPZMatWithMem<TPZOtiTopoDensity>*>(celmem->Material());
        //     if(!mat) DebugStop();
        //     TPZOtiTopoDensity &densstruct = mat->MemItem(0);
        //     densstruct.fDen = rhovecnew[*it];
        // }
        // gel->Divide(sons);
    }
    if(elset.size()) CheckRefinedNeighbors(gmesh, rhovecnew);
    std::cout << "Number of elements in gmesh after refine: " << gmesh->NElements() << std::endl;
    cout << "Number of elements in cmesh after refine: " << cmesh->NElements() << endl;

    // Redo data structure
    if(elset.size()) {
        isNewMesh = true;
        cmesh->InitializeBlock();
        cmesh->ExpandSolution();
        cmesh->CleanUpUnconnectedNodes();
        cmesh->AdjustBoundaryElements();

        // delete cmesh;
        // TElasticity2DAnalytic *elas = new TElasticity2DAnalytic;
        // elas->gE = 1.;//206.8150271873455;
        // elas->gPoisson = 0.3;
        // elas->fProblemType = TElasticity2DAnalytic::EStretchx;
        // cmesh = CreateH1CMesh(gmesh,1,elas);
        // an.SetCompMesh(cmesh,true);
    }
    if(elset.size()) InitializeElemSolOfRefElements(cmesh);

    cout << "Number of elements in cmesh after refine: " << cmesh->NElements() << endl;
    if(elset.size() > 0){
        std::ofstream out("gmesh_ref.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
    return isNewMesh;
}

void CheckRefinedNeighbors(TPZGeoMesh* gmesh, TPZVec<REAL> &rhovecnew){
    const int dim = gmesh->Dimension();
    for(int i = 0 ; i < gmesh->NElements() ; i++){
        TPZGeoEl* gel = gmesh->ElementVec()[i];
        if(!gel) continue;
        if(gel->Dimension() != dim) continue;
        if(gel->HasSubElement()) continue;
        const int nedges = gel->NSides() - gel->NCornerNodes() - 1;
        int nneighref = 0;
        // Check if all neighbors through side of dimension 1 have subelements        
        for (int j = gel->FirstSide(dim-1); j < gel->NSides(); j++) {
            TPZGeoElSide gelside(gel,j);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                if (neighbour.Element()->Dimension() != dim) {                    
                    neighbour = neighbour.Neighbour();
                    continue; // boundary condition
                }
                if (neighbour.Element()->HasSubElement()) nneighref++;
                neighbour = neighbour.Neighbour();
            }
        }
        if(nneighref >= nedges-1){
            // Refine the compel of gel
            TPZCompEl* cel = gel->Reference();
            const int elindex = gel->Index();
            std::cout << "Refining element " << gel->Index() << std::endl;
            if(!cel) DebugStop();
            TPZVec<int64_t> subindexes;
            const bool interpolatesol = false;
            cel->Divide(cel->Index(),subindexes,interpolatesol);
            UpdateSonsMemory(gmesh->Reference(),subindexes,rhovecnew[elindex]);
        }

    }
}

void UpdateSonsMemory(TPZCompMesh* cmesh, TPZVec<int64_t>& subindexes, REAL rhofather) {
    for (int64_t ison = 0; ison < subindexes.size(); ison++) {
        TPZCompEl* celson = cmesh->Element(subindexes[ison]);
        TPZCompElWithMem<TPZCompElH1<pzshape::TPZShapeQuad>>* celmem = dynamic_cast<TPZCompElWithMem<TPZCompElH1<pzshape::TPZShapeQuad>>*>(celson);
        TPZCompElWithMem <TPZCompElH1<pzshape::TPZShapeTriang> > *celmemtri = dynamic_cast<TPZCompElWithMem <TPZCompElH1<pzshape::TPZShapeTriang> >*>(celson);
        if (!celmem && !celmemtri) DebugStop();
        TPZMatWithMem<TPZOtiTopoDensity>* mat = nullptr;
        if(celmem) mat = dynamic_cast<TPZMatWithMem<TPZOtiTopoDensity>*>(celmem->Material());
        else mat = dynamic_cast<TPZMatWithMem<TPZOtiTopoDensity>*>(celmemtri->Material());
        if (!mat) DebugStop();
        TPZOtiTopoDensity& densstruct = mat->MemItem(0);
        densstruct.fDen = rhofather;
    }
}

void InitializeElemSolOfRefElements(TPZCompMesh* cmesh) {    
    TPZFMatrix<STATE> &elementSol = cmesh->ElementSolution();
    elementSol.Resize(cmesh->NElements(), 1);
    for(int64_t i = 0 ; i < cmesh->NElements() ; i++) {
        TPZCompEl* cel = cmesh->Element(i);
        if(!cel) continue;
        TPZCompElWithMem <TPZCompElH1<pzshape::TPZShapeQuad> > *celmem = dynamic_cast<TPZCompElWithMem <TPZCompElH1<pzshape::TPZShapeQuad> >*>(cel);
        TPZCompElWithMem <TPZCompElH1<pzshape::TPZShapeTriang> > *celmemtri = dynamic_cast<TPZCompElWithMem <TPZCompElH1<pzshape::TPZShapeTriang> >*>(cel);
        if(!celmem && !celmemtri) continue;
        TPZVec<int64_t> indices;
        if(celmem) celmem->GetMemoryIndices(indices);
        else celmemtri->GetMemoryIndices(indices);
        for(int j = 1 ; j < indices.size() ; j++) {
            if (indices[0] != indices[j]) {
                DebugStop(); // assuming same memory for the whole element!
            }
        }
        TPZMatWithMem<TPZOtiTopoDensity>* mat = nullptr;
        if(celmem) mat = dynamic_cast<TPZMatWithMem<TPZOtiTopoDensity>*>(celmem->Material());
        else mat = dynamic_cast<TPZMatWithMem<TPZOtiTopoDensity>*>(celmemtri->Material());
        if(!mat) DebugStop();
        TPZOtiTopoDensity &densstruct = mat->MemItem(indices[0]);
        const STATE dens = densstruct.fDen;
        elementSol(cel->Index(),0) = dens;
    }    
}

void CreateNeighVec(TPZVec<std::set<int>>& neighVec, TPZGeoMesh* gmesh) {
    TPZCompMesh* cmesh = gmesh->Reference();
    if(!cmesh) DebugStop(); // should be created if this function is being called
    neighVec.Resize(cmesh->NElements());
    // Fill vector by looping of elements and getting all neighbors of each element
    for(int i = 0 ; i < cmesh->NElements() ; i++){
        TPZCompEl* cel = cmesh->ElementVec()[i];
        if(!cel) continue;        
        TPZGeoEl* gel = cel->Reference();        
        if(!gel) DebugStop();
        if(gel->HasSubElement()) DebugStop();
        if (gel->Dimension() == gmesh->Dimension()){
            for (int j = 0; j < gel->NCornerNodes(); j++) {
                TPZGeoElSide gelside(gel,j);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while (neighbour != gelside) {
                    TPZCompEl* celneigh = neighbour.Element()->Reference();
                    if (neighbour.Element()->Dimension() == gmesh->Dimension() &&  neighbour.Element()->Index() != i && celneigh) {                             
                        neighVec[i].insert(celneigh->Index());
                    }
                    neighbour = neighbour.Neighbour();
                }
            }
        }
    }
}