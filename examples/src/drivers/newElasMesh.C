
/**
  @file newElasSolver.C
  @author Rahul Sampath, rahul.sampath@gmail.com
  */

#include "mpi.h"
#include "petsc.h"
#include "sys.h"
#include <vector>
#include "TreeNode.h"
#include "parUtils.h"
#include "omg.h"
#include "odaUtils.h"
#include "handleStencils.h"
#include "elasticityJac.h"
#include "omgJac.h"
#include <cstdlib>
#include <cstring>
#include "colors.h"
#include "externVars.h"
#include "dendro.h"

static char help[] = "Solves static Navier-Lame equations\
                      using FEM, MG, DA and Matrix-Free methods.\n\n";

//Don't want time to be synchronized. Need to check load imbalance.
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

double**** GradDivType2Stencil; 

double gaussian(double mean, double std_deviation);

int main(int argc, char ** argv ) {	
  int size, rank;
  bool incCorner = 1;  
  unsigned int local_num_pts = 5000;
  double gSize[3];
  unsigned int ptsLen;
  unsigned int maxNumPts = 1;
  unsigned int dim = 3;
  unsigned int maxDepth = 30;
  bool compressLut = false;
  double localTime, globalTime;
  double startTime, endTime;
  DendroIntL localSz, totalSz;
  std::vector<ot::TreeNode> linOct, balOct;
  std::vector<double> pts;
  double mgLoadFac = 2.0;
  char partitionFileNameBase[100];
  bool dumpPartitions = false;

  PetscInitialize(&argc,&argv,"optionsElasticity",help);
  ot::RegisterEvents();

  ot::DAMG_Initialize(MPI_COMM_WORLD);

#ifdef PETSC_USE_LOG
  PetscLogEventRegister("ODAmatDiag",PETSC_VIEWER_COOKIE, &Jac1DiagEvent);
  PetscLogEventRegister("ODAmatMult",PETSC_VIEWER_COOKIE, &Jac1MultEvent);
  PetscLogEventRegister("ODAmatDiagFinest",PETSC_VIEWER_COOKIE, &Jac1FinestDiagEvent);
  PetscLogEventRegister("ODAmatMultFinest",PETSC_VIEWER_COOKIE, &Jac1FinestMultEvent);

  PetscLogEventRegister("OMGmatDiag-2",PETSC_VIEWER_COOKIE, &Jac2DiagEvent);
  PetscLogEventRegister("OMGmatMult-2",PETSC_VIEWER_COOKIE, &Jac2MultEvent);
  PetscLogEventRegister("OMGmatDiagFinest-2",PETSC_VIEWER_COOKIE, &Jac2FinestDiagEvent);
  PetscLogEventRegister("OMGmatMultFinest-2",PETSC_VIEWER_COOKIE, &Jac2FinestMultEvent);

  PetscLogEventRegister("OMGmatDiag-3",PETSC_VIEWER_COOKIE, &Jac3DiagEvent);
  PetscLogEventRegister("OMGmatMult-3",PETSC_VIEWER_COOKIE, &Jac3MultEvent);
  PetscLogEventRegister("OMGmatDiagFinest-3",PETSC_VIEWER_COOKIE, &Jac3FinestDiagEvent);
  PetscLogEventRegister("OMGmatMultFinest-3",PETSC_VIEWER_COOKIE, &Jac3FinestMultEvent);
  PetscLogEventRegister("vMassDiag",PETSC_VIEWER_COOKIE, &vecMassDiagEvent);
  PetscLogEventRegister("vMassMult",PETSC_VIEWER_COOKIE, &vecMassMultEvent);
  PetscLogEventRegister("vMassDiagFinest",PETSC_VIEWER_COOKIE, &vecMassFinestDiagEvent);
  PetscLogEventRegister("vMassMultFinest",PETSC_VIEWER_COOKIE, &vecMassFinestMultEvent);

  PetscLogEventRegister("elDiag",PETSC_VIEWER_COOKIE, &elasticityDiagEvent);
  PetscLogEventRegister("elMult",PETSC_VIEWER_COOKIE, &elasticityMultEvent);
  PetscLogEventRegister("elDiagFinest",PETSC_VIEWER_COOKIE, &elasticityFinestDiagEvent);
  PetscLogEventRegister("elMultFinest",PETSC_VIEWER_COOKIE, &elasticityFinestMultEvent);

  int stages[3];
  PetscLogStageRegister("P2O.",&stages[0]);
  PetscLogStageRegister("Bal",&stages[1]);  
  PetscLogStageRegister("Solve",&stages[2]);  
#endif

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  //  Usage: <exe> local_num_pts[5000] dim[3] maxDepth[30] maxNumPtsPerOctant[1] 
  //  incCorner[1] compressLut[0] mgLoadFac[2.0] partitionFileNameBase  

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
    dumpPartitions = true;
    strcpy(partitionFileNameBase, argv[8]);
  }

  //Generate Points
  if(!rank){
    std::cout << "Creating points on the fly... "<<std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  startTime = MPI_Wtime();

  srand48(0x12345678 + 76543*rank);

  pts.resize(3*local_num_pts);
  for(int i = 0; i < (3*local_num_pts); i++) {
    pts[i]= gaussian(0.5, 0.16);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  endTime = MPI_Wtime();
  localTime = endTime - startTime;

  if(!rank){
    std::cout << " finished creating points  "<<std::endl; // Point size
    std::cout <<"point generation time: "<<localTime << std::endl;
  }

  ptsLen = pts.size();
  std::vector<ot::TreeNode> tmpNodes;
  for(int i = 0;i < ptsLen; i+=3) {
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

  par::removeDuplicates<ot::TreeNode>(tmpNodes,false,MPI_COMM_WORLD);	
  linOct = tmpNodes;
  tmpNodes.clear();
  par::partitionW<ot::TreeNode>(linOct, NULL,MPI_COMM_WORLD);
  // reduce and only print the total ...
  localSz = linOct.size();
  par::Mpi_Reduce<DendroIntL>(&localSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);	
  if(rank==0) {
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

  //Points2Octree....
  MPI_Barrier(MPI_COMM_WORLD);	
#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[0]);
#endif
  startTime = MPI_Wtime();
  ot::points2Octree(pts, gSize, linOct, dim, maxDepth, maxNumPts, MPI_COMM_WORLD);
  endTime = MPI_Wtime();
#ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif
  localTime = endTime - startTime;
  par::Mpi_Reduce<double>(&localTime, &globalTime, 1, MPI_MAX, 0, MPI_COMM_WORLD);
  if(!rank){
    std::cout <<"P2n Time: "<<globalTime << std::endl;
  }
  // reduce and only print the total ...
  localSz = linOct.size();
  par::Mpi_Reduce<DendroIntL>(&localSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  if(rank==0) {
    std::cout<<"# of Unbalanced Octants: " << totalSz<<std::endl;
  }
  pts.clear();

  //Balancing...
  MPI_Barrier(MPI_COMM_WORLD);	
#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[1]);
#endif
  startTime = MPI_Wtime();
  ot::balanceOctree (linOct, balOct, dim, maxDepth, incCorner, MPI_COMM_WORLD, NULL, NULL);
  endTime = MPI_Wtime();
#ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif
  linOct.clear();
  // compute total inp size and output size
  localSz = balOct.size();
  localTime = endTime - startTime;
  par::Mpi_Reduce<DendroIntL>(&localSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  par::Mpi_Reduce<double>(&localTime, &globalTime, 1, MPI_MAX, 0, MPI_COMM_WORLD);

  if(!rank) {
    std::cout << "# of Balanced Octants: "<< totalSz << std::endl;       
    std::cout << "bal Time: "<<globalTime << std::endl;
  }

  //Solve ...
  MPI_Barrier(MPI_COMM_WORLD);	
#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[2]);
#endif

  ot::DAMG       *damg;    
  int       nlevels = 1; //number of multigrid levels
  PetscInt       numRefinements = 0;
  unsigned int   dof = 3; // degrees of freedom per node  

  iC(PetscOptionsGetInt(0,"-numRefinements",&numRefinements,0));
  for(int i = 0; i < numRefinements; i++) {
    std::vector<ot::TreeNode> tmpOct = balOct;
    balOct.clear();
    ot::refineOctree(tmpOct, balOct); 
  }

  PetscInt nlevelsPetscInt = nlevels;
  PetscOptionsGetInt(0, "-nlevels", &nlevelsPetscInt, 0);
  nlevels = nlevelsPetscInt;

  // Note: The user context for all levels will be set separately later.
  MPI_Barrier(MPI_COMM_WORLD);	

  if(!rank) {
    std::cout<<"nlevels initial: "<<nlevels<<std::endl;
  }

  ot::DAMGCreateAndSetDA(PETSC_COMM_WORLD, nlevels, NULL, &damg, 
      balOct, dof, mgLoadFac, compressLut, incCorner);

  balOct.clear();

  if(!rank) {
    std::cout<<"nlevels final: "<<nlevels<<std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);	

  if(!rank) {
    std::cout << "Created DA for all levels."<< std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  ot::PrintDAMG(damg);
  MPI_Barrier(MPI_COMM_WORLD);

  //Dump Partitions
  if(dumpPartitions) {
    for(int i = 0; i < nlevels; i++) {
      char tmpFileName[200];
      sprintf(tmpFileName,"%s_Lev%d.vtk", partitionFileNameBase, i);
      ot::writePartitionVTK(damg[i]->da, tmpFileName);
      if(damg[i]->da_aux) {
        sprintf(tmpFileName,"%s_Lev%d_Aux.vtk", partitionFileNameBase, i);
        ot::writePartitionVTK(damg[i]->da_aux, tmpFileName);
      }
    }
  }

  iC(DAMGDestroy(damg));

  if (!rank) {
    std::cout << GRN << "Destroyed DAMG" << NRM << std::endl;
  }

  endTime = MPI_Wtime();
#ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif

  if (!rank) {
    std::cout << GRN << "Finalizing PETSC" << NRM << std::endl;
  }
  ot::DAMG_Finalize();
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


