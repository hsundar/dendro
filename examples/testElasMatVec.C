
#include "mpi.h"
#include "petsc.h"
#include "sys.h"
#include "parUtils.h"
#include "octUtils.h"
#include "TreeNode.h"
#include "omg.h"
#include "handleStencils.h"
#include "elasticityJac.h"
#include "omgJac.h"
#include "externVars.h"
#include "dendro.h"

#ifdef MPI_WTIME_IS_GLOBAL
#undef MPI_WTIME_IS_GLOBAL
#endif

#ifndef iC
#define iC(fun) {CHKERRQ(fun);}
#endif

#ifdef PETSC_USE_LOG
int elasticityDiagEvent;
int elasticityMultEvent;
int elasticityFinestDiagEvent;
int elasticityFinestMultEvent;

int vecMassDiagEvent;
int vecMassMultEvent;
int vecMassFinestDiagEvent;
int vecMassFinestMultEvent;
#endif

#ifdef PETSC_USE_LOG
//user-defined variables
int Jac1DiagEvent;
int Jac1MultEvent;
int Jac1FinestDiagEvent;
int Jac1FinestMultEvent;

int Jac2DiagEvent;
int Jac2MultEvent;
int Jac2FinestDiagEvent;
int Jac2FinestMultEvent;

int Jac3DiagEvent;
int Jac3MultEvent;
int Jac3FinestDiagEvent;
int Jac3FinestMultEvent;
#endif

double***** LaplacianType1Stencil; 
double**** LaplacianType2Stencil; 
double***** MassType1Stencil; 
double**** MassType2Stencil; 
double****** ShapeFnStencil;
double**** ShapeFnCoeffs;

double**** GradDivType2Stencil; 

double gaussian(double mean, double std_deviation);

int main(int argc, char ** argv ) {	
  int size, rank;
  double startTime, endTime;
  double localTime, globalTime;
  double gSize[3];
  unsigned int local_num_pts = 5000;
  unsigned int dim = 3;
  unsigned int maxDepth = 30;
  unsigned int maxNumPts = 1;
  bool writePts = false;
  bool incCorner = 1;  
  bool compressLut = false;
  double mgLoadFac = 2.0;
  unsigned int numLoops = 100;
  std::vector<ot::TreeNode> linOct;
  std::vector<ot::TreeNode> balOct;
  std::vector<ot::TreeNode> tmpNodes;
  std::vector<double> pts;
  unsigned int ptsLen;
  DendroIntL localSz, totalSz;

  PetscInitialize(&argc, &argv, "options", NULL);
  ot::RegisterEvents();

#ifdef PETSC_USE_LOG
  PetscLogEventRegister(&Jac1DiagEvent,"ODAmatDiag",PETSC_VIEWER_COOKIE);
  PetscLogEventRegister(&Jac1MultEvent,"ODAmatMult",PETSC_VIEWER_COOKIE);
  PetscLogEventRegister(&Jac1FinestDiagEvent,"ODAmatDiagFinest",PETSC_VIEWER_COOKIE);
  PetscLogEventRegister(&Jac1FinestMultEvent,"ODAmatMultFinest",PETSC_VIEWER_COOKIE);

  PetscLogEventRegister(&Jac2DiagEvent,"OMGmatDiag-2",PETSC_VIEWER_COOKIE);
  PetscLogEventRegister(&Jac2MultEvent,"OMGmatMult-2",PETSC_VIEWER_COOKIE);
  PetscLogEventRegister(&Jac2FinestDiagEvent,"OMGmatDiagFinest-2",PETSC_VIEWER_COOKIE);
  PetscLogEventRegister(&Jac2FinestMultEvent,"OMGmatMultFinest-2",PETSC_VIEWER_COOKIE);

  PetscLogEventRegister(&Jac3DiagEvent,"OMGmatDiag-3",PETSC_VIEWER_COOKIE);
  PetscLogEventRegister(&Jac3MultEvent,"OMGmatMult-3",PETSC_VIEWER_COOKIE);
  PetscLogEventRegister(&Jac3FinestDiagEvent,"OMGmatDiagFinest-3",PETSC_VIEWER_COOKIE);
  PetscLogEventRegister(&Jac3FinestMultEvent,"OMGmatMultFinest-3",PETSC_VIEWER_COOKIE);
  PetscLogEventRegister(&vecMassDiagEvent,"vMassDiag",PETSC_VIEWER_COOKIE);
  PetscLogEventRegister(&vecMassMultEvent,"vMassMult",PETSC_VIEWER_COOKIE);
  PetscLogEventRegister(&vecMassFinestDiagEvent,"vMassDiagFinest",PETSC_VIEWER_COOKIE);
  PetscLogEventRegister(&vecMassFinestMultEvent,"vMassMultFinest",PETSC_VIEWER_COOKIE);

  PetscLogEventRegister(&elasticityDiagEvent,"elDiag",PETSC_VIEWER_COOKIE);
  PetscLogEventRegister(&elasticityMultEvent,"elMult",PETSC_VIEWER_COOKIE);
  PetscLogEventRegister(&elasticityFinestDiagEvent,"elDiagFinest",PETSC_VIEWER_COOKIE);
  PetscLogEventRegister(&elasticityFinestMultEvent,"elMultFinest",PETSC_VIEWER_COOKIE);

  int stages[3];
  PetscLogStageRegister(&stages[0],"P2O.");
  PetscLogStageRegister(&stages[1],"Bal");  
  PetscLogStageRegister(&stages[2],"Solve");  
#endif

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //  Usage: <exe> local_num_pts[5000] dim[3] maxDepth[30] maxNumPtsPerOctant[1] 
  //  incCorner[1] compressLut[0] mgLoadFac[2.0] numLoops[100] 

  if(argc > 1) {
    local_num_pts = atoi(argv[1]);
  }  
  if(argc > 2) {
    dim = atoi(argv[2]);
  }
  if(argc > 3) {
    maxDepth = atoi(argv[3]);
  }
  if(argc > 4) {
    maxNumPts = atoi(argv[4]);
  }
  if(argc > 5) {
    incCorner = (bool)(atoi(argv[5]));
  }  
  if(argc > 6) { 
    compressLut = (bool)(atoi(argv[6]));
  }
  if(argc > 7) { 
    mgLoadFac = atof(argv[7]);
  }
  if(argc > 8) {
    numLoops = atoi(argv[8]);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if(!rank){
    std::cout << "Creating points on the fly... "<<std::endl;
  }

  startTime = MPI_Wtime();

  srand48(0x12345678 + 76543*rank);

  pts.resize(3*local_num_pts);
  for(unsigned int i = 0; i < (3*local_num_pts); i++) {
    pts[i]= gaussian(0.5, 0.16);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if(!rank){
    std::cout << " finished creating points  "<<std::endl; // Point size
  }
  endTime = MPI_Wtime();
  localTime = endTime - startTime;
  if(!rank){
    std::cout <<"point generation time: "<<localTime << std::endl;
  }

  ptsLen = pts.size();

  if(writePts) {
    char ptsFileName[100];
    sprintf(ptsFileName, "tempPts%d_%d.pts", rank, size);
    ot::writePtsToFile(ptsFileName, pts);
  }

  for(unsigned int i = 0; i < ptsLen; i+=3) {
    if( (pts[i] > 0.0) &&
        (pts[i+1] > 0.0)
        && (pts[i+2] > 0.0) &&
        ( ((unsigned int)(pts[i]*((double)(1u << maxDepth)))) < (1u << maxDepth))  &&
        ( ((unsigned int)(pts[i+1]*((double)(1u << maxDepth)))) < (1u << maxDepth))  &&
        ( ((unsigned int)(pts[i+2]*((double)(1u << maxDepth)))) < (1u << maxDepth)) ) {
      tmpNodes.push_back( ot::TreeNode((unsigned int)(pts[i]*(double)(1u << maxDepth)),
            (unsigned int)(pts[i+1]*(double)(1u << maxDepth)),
            (unsigned int)(pts[i+2]*(double)(1u << maxDepth)),
            maxDepth,dim,maxDepth) );
    }
  }
  pts.clear();

  MPI_Barrier(MPI_COMM_WORLD);

  if(!rank){
    std::cout <<"Removing bad points... "<<localTime << std::endl;
  }

  par::removeDuplicates<ot::TreeNode>(tmpNodes, false, MPI_COMM_WORLD);
  linOct = tmpNodes;
  tmpNodes.clear();

  MPI_Barrier(MPI_COMM_WORLD);

  if(!rank){
    std::cout <<"Partitioning Input... "<<localTime << std::endl;
  }

  par::partitionW<ot::TreeNode>(linOct, NULL, MPI_COMM_WORLD);

  // reduce and only print the total ...
  localSz = linOct.size();
  par::Mpi_Reduce<DendroIntL>(&localSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

  if(!rank) {
    std::cout<<"# pts= " << totalSz<<std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  pts.resize(3*(linOct.size()));
  ptsLen = (3*(linOct.size()));
  for(int i = 0; i < linOct.size(); i++) {
    pts[3*i] = (((double)(linOct[i].getX())) + 0.5)/((double)(1u << maxDepth));
    pts[(3*i)+1] = (((double)(linOct[i].getY())) +0.5)/((double)(1u << maxDepth));
    pts[(3*i)+2] = (((double)(linOct[i].getZ())) +0.5)/((double)(1u << maxDepth));
  }//end for i
  linOct.clear();

  gSize[0] = 1.0;
  gSize[1] = 1.0;
  gSize[2] = 1.0;

  MPI_Barrier(MPI_COMM_WORLD);	

#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[0]);
#endif
  startTime = MPI_Wtime();
  ot::points2Octree(pts, gSize, linOct, dim, maxDepth, maxNumPts, MPI_COMM_WORLD);
  endTime = MPI_Wtime();
  localTime = endTime - startTime;
#ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif
  par::Mpi_Reduce<double>(&localTime, &globalTime, 1, MPI_MAX, 0, MPI_COMM_WORLD);
  if(!rank){
    std::cout <<"P2n Time: "<<globalTime << "secs " << std::endl;
  }
  pts.clear();

  localSz = linOct.size();
  par::Mpi_Reduce<DendroIntL>(&localSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  if(rank==0) {
    std::cout<<"linOct.size = " << totalSz<<std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);	

#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[1]);
#endif
  startTime = MPI_Wtime();
  ot::balanceOctree (linOct, balOct, dim, maxDepth, incCorner, MPI_COMM_WORLD, NULL, NULL);
  endTime = MPI_Wtime();
  localTime = endTime - startTime;
#ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif
  par::Mpi_Reduce<double>(&localTime, &globalTime, 1, MPI_MAX, 0, MPI_COMM_WORLD);
  if(!rank){
    std::cout <<"Bal Time: "<<globalTime << "secs " << std::endl;
  }
  linOct.clear();

  localSz = balOct.size();
  par::Mpi_Reduce<DendroIntL>(&localSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  if(rank==0) {
    std::cout<<"balOct.size = " << totalSz<<std::endl;
  }

#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[2]);
#endif
  ot::DAMG_Initialize(MPI_COMM_WORLD);

  ot::DAMG       *damg;    
  int       nlevels = 1; //number of multigrid levels
  unsigned int   dof = 3; // degrees of freedom per node  

  ot::DAMGCreateAndSetDA(PETSC_COMM_WORLD, nlevels, NULL, &damg, 
      balOct, dof, mgLoadFac, compressLut, incCorner);

  balOct.clear();

  MPI_Barrier(MPI_COMM_WORLD);
  ot::PrintDAMG(damg);
  MPI_Barrier(MPI_COMM_WORLD);

  createLmatType2(LaplacianType2Stencil);
  createMmatType2(MassType2Stencil);
  createGDmatType2(GradDivType2Stencil);

  ot::DAMGCreateSuppressedDOFs(damg);

  SetElasticityContexts(damg);

  ot::getPrivateMatricesForKSP_Shell = getPrivateMatricesForKSP_Shell_Elas;

  VecZeroEntries(DAMGGetx(damg));

  DAMGSetKSP(damg, CreateElasticityMat, ComputeElasticityMat, ComputeRHS4);

  for(unsigned int i = 0; i < numLoops; i++) {
    MatMult(DAMGGetJ(damg), DAMGGetx(damg), DAMGGetRHS(damg));
  }

  destroyLmatType2(LaplacianType2Stencil);
  destroyMmatType2(MassType2Stencil);
  destroyGDmatType2(GradDivType2Stencil);

  DestroyElasticityContexts(damg);

  DAMGDestroy(damg);

  ot::DAMG_Finalize();

#ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif

  PetscFinalize();
}//end main

double gaussian(double mean, double std_deviation) {
  static double t1 = 0, t2=0;
  double x1, x2, x3, r;

  using namespace std;

  // reuse previous calculations
  if(t1) {
    const double tmp = t1;
    t1 = 0;
    return mean + std_deviation * tmp;
  }
  if(t2) {
    const double tmp = t2;
    t2 = 0;
    return mean + std_deviation * tmp;
  }

  // pick randomly a point inside the unit disk
  do {
    x1 = 2 * drand48() - 1;
    x2 = 2 * drand48() - 1;
    x3 = 2 * drand48() - 1;
    r = x1 * x1 + x2 * x2 + x3*x3;
  } while(r >= 1);

  // Box-Muller transform
  r = sqrt(-2.0 * log(r) / r);

  // save for next call
  t1 = (r * x2);
  t2 = (r * x3);

  return mean + (std_deviation * r * x1);
}//end gaussian


