
/**
  @file elasticitySolver.C
  @author Rahul Sampath, rahul.sampath@gmail.com
  */

#include "mpi.h"
#include "petsc.h"
#include "sys.h"
#include <vector>
#include "TreeNode.h"
#include "parUtils.h"
#include "omg.h"
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
double**** ShapeFnCoeffs;

double**** GradDivType2Stencil; 

int main(int argc, char ** argv ) {	
  int size, rank;
  bool incCorner = 1;  
  unsigned int numPts;
  unsigned int solveU = 0;
  unsigned int writeB = 0;  
  char Kstr[20];
  char pFile[50],bFile[50],uFile[50];
  double gSize[3];
  unsigned int ptsLen;
  unsigned int maxNumPts= 1;
  unsigned int dim=3;
  unsigned int maxDepth=30;
  bool compressLut=true;
  double localTime, totalTime;
  double startTime, endTime;
  DendroIntL locSz, totalSz;
  std::vector<ot::TreeNode> linOct, balOct;
  std::vector<double> pts;
  double mgLoadFac = 2.0;

  PetscInitialize(&argc,&argv,"optionsElasticity",help);
  ot::RegisterEvents();

  ot::DAMG_Initialize(MPI_COMM_WORLD);

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

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(argc < 2) {
    std::cerr << "Usage: " << argv[0] << "inpfile  maxDepth[30] solveU[0]\
      writeB[0] dim[3] maxNumPtsPerOctant[1] incCorner[1] compressLut[1] mgLoadFac[2.0] " << std::endl;
    return -1;
  }
  if(argc > 2) {
    maxDepth = atoi(argv[2]);
  }
  if(argc > 3) {
    solveU = atoi(argv[3]);
  }
  if(argc > 4) {
    writeB = atoi(argv[4]);
  }
  if(argc > 5) {
    dim = atoi(argv[5]);
  }
  if(argc > 6) {
    maxNumPts = atoi(argv[6]);
  }
  if(argc > 7) { incCorner = (bool)(atoi(argv[7]));}  
  if(argc > 8) { compressLut = (bool)(atoi(argv[8]));}
  if(argc > 9) { mgLoadFac = atof(argv[9]); }

  strcpy(bFile,argv[1]);
  ot::int2str(rank,Kstr);
  strcat(bFile,Kstr);
  strcat(bFile,"_\0");
  ot::int2str(size,Kstr);
  strcat(bFile,Kstr);
  strcpy(pFile,bFile);
  strcpy(uFile,bFile);
  strcat(bFile,"_Bal.ot\0");
  strcat(pFile,".pts\0");
  strcat(uFile,".sol\0");

  //Points2Octree....
  MPI_Barrier(MPI_COMM_WORLD);	
  if(!rank){
    std::cout << " reading  "<<pFile<<std::endl; // Point size
  }
  ot::readPtsFromFile(pFile, pts);
  if(!rank){
    std::cout << " finished reading  "<<pFile<<std::endl; // Point size
  }
  MPI_Barrier(MPI_COMM_WORLD);	
  ptsLen = pts.size();
  std::vector<ot::TreeNode> tmpNodes;
  for(int i=0;i<ptsLen;i+=3) {
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
  locSz = linOct.size();
  par::Mpi_Reduce<DendroIntL>(&locSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);	
  if(rank==0) {
    std::cout<<"# pts= " << totalSz<<std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);	

  pts.resize(3*(linOct.size()));
  ptsLen = (3*(linOct.size()));
  for(int i=0;i<linOct.size();i++) {
    pts[3*i] = (((double)(linOct[i].getX())) + 0.5)/((double)(1u << maxDepth));
    pts[(3*i)+1] = (((double)(linOct[i].getY())) +0.5)/((double)(1u << maxDepth));
    pts[(3*i)+2] = (((double)(linOct[i].getZ())) +0.5)/((double)(1u << maxDepth));
  }//end for i
  linOct.clear();
  gSize[0] = 1.;
  gSize[1] = 1.;
  gSize[2] = 1.;

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
  par::Mpi_Reduce<double>(&localTime, &totalTime, 1, MPI_MAX, 0, MPI_COMM_WORLD);
  if(!rank){
    std::cout <<"P2n Time: "<<totalTime << std::endl;
  }
  // reduce and only print the total ...
  locSz = linOct.size();
  par::Mpi_Reduce<DendroIntL>(&locSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
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
  if(writeB) { 
    ot::writeNodesToFile(bFile,balOct);
  }
  // compute total inp size and output size
  locSz = balOct.size();
  localTime = endTime - startTime;
  par::Mpi_Reduce<DendroIntL>(&locSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  par::Mpi_Reduce<double>(&localTime, &totalTime, 1, MPI_MAX, 0, MPI_COMM_WORLD);

  if(!rank) {
    std::cout << "# of Balanced Octants: "<< totalSz << std::endl;       
    std::cout << "bal Time: "<<totalTime << std::endl;
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

  MPI_Barrier(MPI_COMM_WORLD);
  if(!rank) {
    std::cout << "Creating stencils..."<< std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  createLmatType2(LaplacianType2Stencil);
  createMmatType2(MassType2Stencil);
  createGDmatType2(GradDivType2Stencil);

  MPI_Barrier(MPI_COMM_WORLD);
  if(!rank) {
    std::cout << "Created all stencils."<< std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  ot::DAMGCreateSuppressedDOFs(damg);

  SetElasticityContexts(damg);

  MPI_Barrier(MPI_COMM_WORLD);
  if(!rank) {
    std::cout << "Set Elasticity Contexts all levels."<< std::endl;
  }

  //Functions for using KSP_Shell (will be used @ the coarsest grid if not all
//processors are active on the coarsest grid)

  ot::getPrivateMatricesForKSP_Shell = getPrivateMatricesForKSP_Shell_Elas;
  
  //Set function pointers so that PC_BlockDiag could be used.

  ot::getDofAndNodeSizeForPC_BlockDiag = getDofAndNodeSizeForElasticityMat;

  ot::computeInvBlockDiagEntriesForPC_BlockDiag = computeInvBlockDiagEntriesForElasticityMat;

  PetscInt rhsType = 2;
  iC(PetscOptionsGetInt(0,"-rhsType",&rhsType,0));

  if(rhsType == 4) {
    ot::DAMGSetKSP(damg, CreateElasticityMat,
        ComputeElasticityMat, ComputeRHS4);
  }else {
    ot::DAMGSetKSP(damg, CreateElasticityMat,
        ComputeElasticityMat, ComputeElasticityRHS);
  }

  PetscReal norm2;
  PetscReal normInf;

  MPI_Barrier(MPI_COMM_WORLD);
  if(!rank) {
    std::cout << "Solving with 0-intial guess"<< std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  startTime = MPI_Wtime();
  iC(ot::DAMGSolve(damg));
  endTime = MPI_Wtime();
  localTime = endTime - startTime;
  par::Mpi_Reduce<double>(&localTime, &totalTime, 1, MPI_MAX, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

  if (!rank) {
    std::cout << GRN << "Done Solve" << NRM << std::endl;
    std::cout << "Solve Time: "<<totalTime << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  destroyLmatType2(LaplacianType2Stencil);
  destroyMmatType2(MassType2Stencil);
  destroyGDmatType2(GradDivType2Stencil);

  if (!rank) {
    std::cout << GRN << "Destroyed Stencils" << NRM << std::endl;
  }

  DestroyElasticityContexts(damg);

  if (!rank) {
    std::cout << GRN << "Destroyed Elasticity Contexts." << NRM << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  iC(DAMGDestroy(damg));

  if (!rank) {
    std::cout << GRN << "Destroyed DAMG" << NRM << std::endl;
  }

  endTime = MPI_Wtime();
#ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif
  balOct.clear();

  if (!rank) {
    std::cout << GRN << "Finalizing PETSC" << NRM << std::endl;
  }
  ot::DAMG_Finalize();
  PetscFinalize();
}//end function

