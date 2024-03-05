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
//#include <Elasticity/TPZElasticity2DOtiTopo.h>
#include <DarcyFlow/TPZDarcyFlow.h>
#include "TPZAnalyticSolution.h"
#include "TPZGeoMeshTools.h"
#include <TPZGmshReader.h>
#include "tpzchangeel.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "pzcompelwithmem.h"
#include "TPZCompElH1.h"
#include "pzshapequad.h"
#include "TPZElementMatrixT.h"

enum EMatid {ENone,EDomain,EBottomBC,ERightBC,ELeftBC,ETopBC,EHoleBC};
const int global_nthread = 8;

TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);
void CreateBCs(TPZGeoMesh* gmesh);
TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord);
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh);

int main() {
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    
    std::cout << "--------- Starting simulation ---------" << std::endl;
    const int pord = 3;
//    TPZGeoMesh* gmesh = ReadMeshFromGmsh("holeplate_coarse.msh");
    TPZGeoMesh* gmesh = ReadMeshFromGmsh("holeplate_fine.msh");
    std::ofstream out("gmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    
    TPZCompMesh* cmeshH1 = CreateH1CMesh(gmesh,pord);

    TPZLinearAnalysis an(cmeshH1);
    SolveProblemDirect(an,cmeshH1);
    PrintResults(an,cmeshH1);
    
    delete cmeshH1;
    delete gmesh;
        
    std::cout << "--------- Simulation finished ---------" << std::endl;
}

TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord) {
    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    const int dim = gmesh->Dimension();
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pord);
    cmesh->SetAllCreateFunctionsContinuous();
    
    
    TPZDarcyFlow *darcy = new TPZDarcyFlow(EDomain,dim);
    darcy->SetConstantPermeability(1.);
    
    auto f = [](const TPZVec<REAL> &x, TPZVec<REAL> &result){
        result[0] = 0.32;
    };
    darcy->SetForcingFunction(f, 4);
    cmesh->InsertMaterialObject(darcy);
    
    TPZFMatrix<STATE> val1(3,3,0.);
    TPZManVector<STATE> val2(3,0.);
    
    const int diri = 0, neu = 1, mixed = 2, normaltrac = 4;
    
    auto* BCHole = darcy->CreateBC(darcy, EHoleBC, diri, val1, val2);
    cmesh->InsertMaterialObject(BCHole);
    
    cmesh->AutoBuild();
    
    return cmesh;
}

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh) {

//    TPZSkylineStructMatrix<STATE> matskl(cmesh);
    TPZSSpStructMatrix<STATE> matskl(cmesh);
    matskl.SetNumThreads(global_nthread);
    an.SetStructuralMatrix(matskl);
    
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);//ELU //ECholesky // ELDLt
    an.SetSolver(step);
    
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

void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh) {
 
    std::cout << "--------- Post Process ---------" << std::endl;
    TPZSimpleTimer postProc("Post processing time");
    const std::string plotfile = "postprocess_darcyzinho";
    constexpr int vtkRes{0};
    TPZVec<std::string> fields = {
        "Pressure",
        "Flux"
    };
    static auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
    vtk.SetNThreads(global_nthread);
    vtk.Do();
    std::cout << "Total time = " << postProc.ReturnTimeDouble()/1000. << " s" << std::endl;
    
    return;
}

TPZGeoMesh* ReadMeshFromGmsh(std::string file_name)
{
    TPZGeoMesh *gmesh;
    gmesh = new TPZGeoMesh();
    {
        TPZGmshReader reader;
        // essa interface permite voce mapear os nomes dos physical groups para
        // o matid que voce mesmo escolher
        TPZManVector<std::map<std::string,int>,4> stringtoint(6);
        stringtoint[2]["Domain"] = EDomain;
        stringtoint[1]["BottomBC"] = EBottomBC;
        stringtoint[1]["RightBC"] = ERightBC;
        stringtoint[1]["TopBC"] = ETopBC;
        stringtoint[1]["LeftBC"] = ELeftBC;
        stringtoint[1]["Hole"] = EHoleBC;
        
        reader.SetDimNamePhysical(stringtoint);
        reader.GeometricGmshMesh(file_name,gmesh,false);
    }
    
    return gmesh;
}
