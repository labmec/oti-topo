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
#include <Elasticity/TPZElasticity3D.h>
#include "TPZAnalyticSolution.h"
#include "TPZGeoMeshTools.h"
#include <TPZGmshReader.h>
#include "tpzchangeel.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"

enum EMatid {ENone,EDomain,EBC};
//Cria a classe EMatid, que é do tipo enum, isso significa que ENone, EDomain e EBC são variáveis que não podem ter seu valor alterado diretamente, preservando-as em seu estado original relacionado a classe enum.
const int global_nthread = 16;
//define uma constante chamada global_nthread com o valor 16. Essa constante provavelmente é usada para controlar o número de threads usadas em operações de paralelismo, como na montagem da matriz estrutural ou no pós-processamento.

TPZGeoMesh* CreateGMesh(int ndiv);
//define a função CreateGMesh, que recebe um valor inteiro ndiv que é retornado como um pointer a TPZGeoMesh.
//função que cria e retorna uma malha geométrica com base no número de divisões ndiv.
TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);
//Analogamente, // Define a função ReadMeshFromGmsh, que recebe uma string relacionada ao nome de um arquivo, que é retornado como um pointer a TPZGeoMesh.
//função que lê uma malha geométrica de um arquivo no formato Gmsh e retorna essa malha.
void CreateBCs(TPZGeoMesh* gmesh);
//Cria a função CreateBCs que não retorna valor, cuja entrada é o pointer para TPZGeoMesh.
//função que cria condições de contorno na malha geométrica gmesh.
TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord, TElasticity3DAnalytic *elas);
//CreateH1CMesh é uma função que cria uma malha computacional e que retorna um pointer a TPZCompMesh, que tem os seguintes parâmetros:
//TPZGeoMesh* gmesh é um pointer para TPZGeoMesh
//const int pord, é uma constante da variável inteira pord e, portanto, não pode ser alterada dentro da função
//TElasticity3DAnalytic *elas é um pointer para TElasticity3DAnalytic
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
//Define a função do tipo void (Não retorna valor) SolveProblemDirect, que contém os seguintes parâmetros:
//TPZLinearAnalysis &an, é uma referência a função TPZLinearAnalysis denominada "an". Desse modo, quando essa referência, denominada "an", é alterada, o objeto original de onde ela foi retirada também será alterado, e não uma cópia desse objeto.
//TPZCompMesh *cmesh é um pointer para um objeto do tipo TPZCompMesh.
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
//É definida a função PrintResults do tipo void, cujos parâmetros são os mesmos já explicados na linha anterior.

//Essas funções serão detalhadas e suas funcionalidades serão expostas mais à frente.

int main() {
    
    std::cout << "--------- Starting simulation ---------" << std::endl;
    //printa o texto "--------- Starting simulation ---------", indicando que a simulação está se iniciando.
    const int pord = 1;
    //é declarada a constante pord, cujo valor inteiro atribuído é 1. essa variável provavelmente representa o grau dos polinômios usados para aproximar a solução dentro dos elementos finitos.
    int ndiv = 2;
    //é declarada a variável ndiv (representa o número de divisões de cada uma das dimensões da malha), cujo valor inteiro atribuído é 2.
    TPZGeoMesh* gmesh = CreateGMesh(ndiv);
    //criação da váriável gmesh (geometric mesh) cujo valor será um ponteiro de um objeto da classe TPZGeoMesh, atribuído através da função CreateGMesh(ndiv), que recebe o valor de ndiv como entrada. 
    //o que ocorre na função está detalhado mais a frente no código.    
    std::ofstream out("gmesh.vtk");
    //a variável out está sendo declarada como pertencente à classe std::ofstream (usada para lidar com arquivos de saída). o nome do arquivo que será aberto para escrita é "gmesh.vtk". 
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    //está sendo chamada a função PrintGMeshVTK da classe TPZVTKGeoMesh, cujos argumentos são "gmesh" e "out".
    //essa função formata os dados da malha geométrica "gmesh" para o formato VTK e os escreve no objeto "out", que representa o arquivo "gmesh.vtk"

    TElasticity3DAnalytic *elas = new TElasticity3DAnalytic;
    //Uma variável chamada elas é declarada como um ponteiro para um objeto da classe TElasticity3DAnalytic. Em seguida, é criada uma instância dessa classe usando o operador new e o endereço do objeto é atribuído à variável elas.
    //TElasticity3DAnalytic é relacionada ao problema de elasticidade.
    elas->fE = 250.;//206.8150271873455;
    //Define o valor do atributo fE, pertencente ao objeto cujo pointer é elas, para 250, provavelmente antes era para 206.8150271873455 e tal valor foi alterado.
    elas->fPoisson = 0.;
    //Define o valor de fPoisson (outro atributo do objeto) para 0. (float)
    elas->fProblemType = TElasticity3DAnalytic::EStretchx;
    //fProblemType recebe TElasticity3DAnalytic::EStretchx, um valor do tipo enum associado ao problema de elasticidade.
    TPZCompMesh* cmeshH1 = CreateH1CMesh(gmesh,pord,elas);
    // Define o pointer chamado cmeshH1 relacionado à classe TPZCompMesh, associando-o a função CreateH1CMesh, cujos parâmetros (já antes declarados) são: gmesh, pord, elas.
    // Está detalhada mais à frente no código. Essa função cria uma malha computacional para resolver o problema de elasticidade tridimensional, define materiais, condições de contorno e constrói a malha com base nas configurações especificadas.

    TPZLinearAnalysis an(cmeshH1);
    //um objeto da classe TPZLinearAnalysis chamado an é criado e associado à malha computacional cmeshH1 que foi previamente criada. A TPZLinearAnalysis é uma classe que faz análises lineares.
    SolveProblemDirect(an,cmeshH1);
    //chama uma função chamada SolveProblemDirect para resolver o problema. A função recebe como argumentos o objeto an que representa a análise linear e a malha computacional cmeshH1 que contém as informações sobre o problema a ser resolvido.
    //detalhadams a frene, basicamene configura a montagem e a resolução do sistema de equações para resolver o problema definido na malha computacional cmesh. Ela usa um solver direto Cholesky para resolver o sistema, paraleliza a montagem da matriz estrutural e mede o tempo gasto em cada etapa do processo.
    std::cout << "--------- PostProcess ---------" << std::endl;
    //printa na tela "--------- PostProcess ---------", indicando que a simulação está em processamento.
    PrintResults(an,cmeshH1);
    //chama a função PrintResults para realizar o pós-processamento dos resultados. Essa função provavelmente gera saídas com os resultados da simulação.
    
    // deleting stuff
    delete cmeshH1;
    //libera a memória associada à malha computacional cmeshH1 criada anteriormente. Isso é importante para evitar vazamento de memória.
    delete gmesh;
    //libera a memória associada à malha geométrica gmesh criada anteriormente.
        
    std::cout << "--------- Simulation finished ---------" << std::endl;
    //esta linha imprime uma mensagem indicando que a simulação foi concluída. "--------- Simulation finished ---------"
}

TPZGeoMesh* CreateGMesh(int ndiv) {
    //está sendo definida a função CreateGMesh que retornará um ponteiro para um objeto da classe TPZGeoMesh, a função recebe o argumento inteiro "ndiv".
    TPZGeoMesh* gmesh = new TPZGeoMesh;
    //gmesh é o endereço de um novo objeto criado pertencente à classe TPZGeoMesh.
    
    MMeshType meshType = MMeshType::EHexahedral;
    //a variável meshType está recebendo o valor EHexahedral, indicando que os elementos da malha serão hexaédricos (cubos - 6 faces).
    int dim = 3;
    //variável "dim" é definida com valor inteiro 3 indicando que a malha será tridimensional.
    TPZManVector<REAL,3> minX = {-1,-1,-1};
    TPZManVector<REAL,3> maxX = {1,1,1};
    //Ambos os vetores "minX" e "maxX" possuem 3 elementos que indicam as coordenadas máximas e mínimas, respectivamente, da malha geométrica em cada uma das 3 dimensões.
    //Nesse caso, a caixa da malha geométrica vai de -1 a 1 nas 3 dimensões.
    //TPZManVector é uma classe utilizada para criar vetores desse tipo.
    int nMats = 2*dim+1;
    // número de materiais talvez.
    
    constexpr bool createBoundEls{true};
    //é atribuído "true" à váriavel booleana "createBoundEls", indicando provavelmente que condições de contorno serão inseridas na malha.
    TPZVec<int> matIds(nMats,EBC);
    //é criado o vetor inteiro matIds com "nMats" elemento(s) cujos valores serão todos "EBC".
    matIds[0] = EDomain;
    //o primeiro elemento (índice 0) do vetor matIds está sendo definido como Edomain.
    
    TPZManVector<int,3> ndivvec = {ndiv,ndiv,ndiv};
    //é criado o vetor inteiro ndivvec com 3 elementos, todos com valor ndiv, indicando que cada uma das 3 dimensões da malha será dividida em ndiv (nesse caso 2) partes, ou seja, 2 elementos.
    gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,matIds, ndivvec, meshType,createBoundEls);
    //é chamada a função CreateGeoMeshOnGrid para criar finalmente a malha geométrica com os parâmetros específicados anteriormente. o resultado é atribuído à variável gmesh.

    return gmesh;
    //a função retorna o ponteiro gmesh que aponta para a malha criada.
}


TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord, TElasticity3DAnalytic *elas) {
    // Cria a função CreateH1CMesh que retornará um pointer associado a classe TPZCompMesh, que aceita os seguintes parâmetros: TPZGeoMesh* gmesh, const int pord, TElasticity3DAnalytic *elas.

    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    //Cria um objeto da classe TPZCompMesh chamado cmesh e o inicializa com a malha geométrica gmesh.
    const int dim = gmesh->Dimension();
    //Determina a dimensão do modelo (dim, uma constante) com base na dimensão da malha geométrica gmesh através do método Dimension() do objeto gmesh.
    cmesh->SetDimModel(dim);
    //configurando a dimensão do modelo na malha computacional cmesh com base no valor da variável dim. o método SetDimModel() foi chamado para definir a dimensão do modelo para a malha computacional.
    cmesh->SetDefaultOrder(pord);
    //chama o método SetDefaultOrder() para configurar a ordem a ser utilizada na aproximação dos elementos finitos da malha.
    cmesh->SetAllCreateFunctionsContinuous();
    //o método SetAllCreateFunctionsContinuous() é chamado para configurar as funções de base a serem utilizadas na resolução do problema.
    
    const STATE E = elas->fE, nu = elas->fPoisson;
    //duas constantes E e nu são definidas. Essas constantes representam o módulo de elasticidade (E) e o coeficiente de Poisson (nu), que são propriedades do material elástico.
    TPZManVector<STATE> force = {0,0,0};
    //um vetor chamado force é criado com três elementos iniciais igual a zero. Esse vetor representa as forças aplicadas no problema, mas no código atual, todas as forças são definidas como zero.
    TPZElasticity3D *mat = new TPZElasticity3D(EDomain, E, nu, force, 0., 0., 0.);
    //criado um objeto mat da classe TPZElasticity3D, que representa o material elástico tridimensional. Ele recebe como argumentos o domínio (EDomain), o módulo de elasticidade E, o coeficiente de Poisson nu, o vetor de forças force, e valores iniciais de outros parâmetros que são todos definidos como zero.
    mat->SetExactSol(elas->ExactSolution(), 2);
    //define a solução exata do problema no material mat. A solução exata é obtida a partir do objeto elas do tipo TElasticity3DAnalytic.
    mat->SetForcingFunction(elas->ForceFunc(), 4);
    //define a função de força no material mat. A função de força é obtida do objeto elas do tipo TElasticity3DAnalytic.
    cmesh->InsertMaterialObject(mat);
    //insere o material mat na malha computacional cmesh.
    
    TPZFMatrix<STATE> val1(3,3,0.);
    TPZManVector<STATE> val2(3,0.);
    //são definidas matrizes e vetores vazios que serão usados posteriormente na criação das condições de contorno.
    
    const int diri = 0, neu = 1, mixed = 2, normaltrac = 4;
    //definidas algumas constantes para identificar os tipos de condições de contorno, como diri, neu, mixed, e normaltrac.
    auto* BCCond0 = mat->CreateBC(mat, EBC, diri, val1, val2);
    //cria uma condição de contorno (BCCond0) usando o método CreateBC do material mat. 
    BCCond0->SetForcingFunctionBC(elas->ExactSolution(), 4);
    //é definida uma função de força de condição de contorno usando elas->ExactSolution().
    cmesh->InsertMaterialObject(BCCond0);
    //insere a condição de contorno BCCond0 na malha computacional cmesh.
    
    cmesh->AutoBuild();
    //este método é chamado para construir automaticamente a malha computacional com base nas configurações e objetos de material e condição de contorno definidos anteriormente.
    
    return cmesh;
    //a função retorna o ponteiro para a malha computacional cmesh
}

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
//Declara a função SolveProblemDirect do tipo void (logo não retornará valor).
{

    TPZSkylineStructMatrix<STATE> matskl(cmesh);
    //é criado um objeto matskl da classe TPZSkylineStructMatrix que representa a matriz estrutural do sistema de equações. O construtor TPZSkylineStructMatrix recebe a malha computacional cmesh como argumento.
    matskl.SetNumThreads(global_nthread);
    //define o número de threads que serão usadas para a montagem da matriz. A variável global_nthread provavelmente contém o número desejado de threads.
    an.SetStructuralMatrix(matskl);
    //associa a matriz estrutural matskl à análise linear an. Isso configura a matriz estrutural que será usada durante a resolução do problema.
    
    TPZStepSolver<STATE> step;
    //cria um objeto step da classe TPZStepSolver que será usado como o solver para resolver o sistema de equações.
    step.SetDirect(ECholesky);//ELU //ECholesky // ELDLt
    //define o solver direto a ser utilizado. Neste caso, o solver direto Cholesky (ECholesky) é escolhido. O Cholesky é um método de fatoração LU que pode ser usado para resolver sistemas de equações lineares.
    an.SetSolver(step);
    //associa o solver step à análise linear an. Isso configura o solver que será usado durante a resolução do problema.
    
    //assembles the system
    std::cout << "--------- Assemble ---------" << std::endl;
    //imprime "--------- Assemble ---------" indicando o início da montagem do sistema de equações.
    TPZSimpleTimer time_ass;
    //inicializa um temporizador para medir o tempo gasto na montagem do sistema.
    an.Assemble();
    //monta o sistema de equações usando a matriz estrutural e outras configurações previamente definidas. Esta é a etapa onde as contribuições de todos os elementos finitos são somadas para formar o sistema global de equações.
    std::cout << "Total time = " << time_ass.ReturnTimeDouble()/1000. << " s" << std::endl;
    //imprime o tempo gasto na montagem do sistema.

    std::cout << "--------- Solve ---------" << std::endl;
    //printa "--------- Solve ---------" indicando o início da resolução do sistema de equações.
    TPZSimpleTimer time_sol;
    //inicializa um temporizador para medir o tempo gasto na resolução do sistema.
    an.Solve();
    //resolve o sistema de equações usando o solver direto configurado anteriormente.
    std::cout << "Total time = " << time_sol.ReturnTimeDouble()/1000. << " s" << std::endl;
    //imprime o tempo gasto na resolução do sistema.
    
    return;
    //a função é concluída e retorna.
}

void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
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
        "Stress",
        "Strain",
    };
    //nesse conjunto de linhas de código, temos que TPZVec é uma estrutura do tipo vetor que contém como argumento uma variável chamda "fields" que é uma lista de strings, que, pelo que se chamam, são relacionadas ao deslocamento, deformação e tensão.
    //cria um vetor de strings chamado fields que contém os nomes dos campos que serão pós-processados. Neste caso, os campos incluem "Displacement" (deslocamento), "Stress" (tensão) e "Strain" (deformação). Esses campos representam propriedades do problema que desejamos visualizar após a simulação.
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
    //essa linha de código declara uma variável chamada vtk do tipo auto, o que significa que o compilador irá deduzir o tipo que ela terá a depender do que ela é igual. No caso, ela é igual a função TPZVTKGenerator, de parâmetros cmesh, fields, plotfile, vtkRes.
    //cria um objeto vtk da classe TPZVTKGenerator, que é usado para gerar arquivos VTK a partir dos dados da malha computacional cmesh. Os argumentos passados para o construtor incluem a malha computacional, os campos a serem pós-processados, o nome base do arquivo de saída (plotfile) e a resolução VTK (vtkRes).
    vtk.SetNThreads(global_nthread);
    //define o número de threads a serem usadas durante o pós-processamento. A variável global_nthread provavelmente contém o número desejado de threads.
    vtk.Do();
    //inicia o processo de geração dos arquivos VTK. Esta função gera arquivos de saída contendo informações sobre os campos especificados na malha computacional.
    std::cout << "Total time = " << postProc.ReturnTimeDouble()/1000. << " s" << std::endl;
    //imprime o tempo gasto no pós-processamento, convertido para segundos.
    
    return;
    //a função é concluída e retorna.
}
