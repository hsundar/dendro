
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

int main(int argc, char ** argv ) {	
  int size, rank;
  bool incCorner = 1;  
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
  char ptsFileName[256];

  PetscInitialize(&argc,&argv,"optionsElasticity",help);
  ot::RegisterEvents();

  ot::DAMG_Initialize(MPI_COMM_WORLD);

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if(argc < 2) {
    std::cout<<"  Usage: <exe> ptsFilePrefix dim[3] maxDepth[30]"<<
      " maxNumPtsPerOctant[1] incCorner[1] compressLut[0]"<<
      " mgLoadFac[2.0] partitionFileNameBase "<<std::endl; 
    ot::DAMG_Finalize();
    PetscFinalize();
  }

  if(argc > 1) {
    sprintf(ptsFileName, "%s%d_%d.pts", argv[1], rank, size);
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


  //Read pts from files
  MPI_Barrier(MPI_COMM_WORLD);	
  if(!rank){
    std::cout << " reading  "<<ptsFileName<<std::endl; // Point size
  }
  ot::readPtsFromFile(ptsFileName, pts);
  if(!rank){
    std::cout << " finished reading  "<<ptsFileName<<std::endl; // Point size
  }
  MPI_Barrier(MPI_COMM_WORLD);	

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

