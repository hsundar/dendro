
#include "mpi.h"
#include "petsc.h"
#include "sys.h"
#include "TreeNode.h"
#include "parUtils.h"
#include "omg.h"
#include "oda.h"
#include "omgJac.h"
#include "handleStencils.h"
#include <cstdlib>
#include <cstring>
#include "colors.h"
#include "externVars.h"
#include "dendro.h"

static char help[] = "Scalar Problem";

//Don't want time to be synchronized. Need to check load imbalance.
#ifdef MPI_WTIME_IS_GLOBAL
#undef MPI_WTIME_IS_GLOBAL
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

  PetscInitialize(&argc,&argv,"options",help);
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

  int stages[3];
  PetscLogStageRegister(&stages[0],"P2O.");
  PetscLogStageRegister(&stages[1],"Bal");  
  PetscLogStageRegister(&stages[2],"Solve");  
#endif

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(argc < 3) {
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
        ( ((unsigned int)(pts[i]*((double)(1u<<maxDepth)))) < (1u<<maxDepth))  &&
        ( ((unsigned int)(pts[i+1]*((double)(1u<<maxDepth)))) < (1u<<maxDepth))  &&
        ( ((unsigned int)(pts[i+2]*((double)(1u<<maxDepth)))) < (1u<<maxDepth)) ) {
      tmpNodes.push_back( ot::TreeNode((unsigned int)(pts[i]*(double)(1u<<maxDepth)),
            (unsigned int)(pts[i+1]*(double)(1u<<maxDepth)),
            (unsigned int)(pts[i+2]*(double)(1u<<maxDepth)),
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
  DendroIntL totalNumPts;
  par::Mpi_Reduce<DendroIntL>(&locSz, &totalNumPts, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);	
  if(rank==0) {
    std::cout<<"# pts= " << totalNumPts<<std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);	

  PetscTruth setMatPropsUsingPts;
  PetscOptionsHasName(0,"-setMatPropsUsingPts",&setMatPropsUsingPts);

  std::vector<ot::TreeNode> matPropNodes;
  if(setMatPropsUsingPts) {
    matPropNodes = linOct;
  }

  pts.resize(3*(linOct.size()));
  ptsLen = (3*(linOct.size()));

  for(int i=0;i<linOct.size();i++) {
    pts[3*i] = (((double)(linOct[i].getX())) + 0.5)/((double)(1u<<maxDepth));
    pts[(3*i)+1] = (((double)(linOct[i].getY())) +0.5)/((double)(1u<<maxDepth));
    pts[(3*i)+2] = (((double)(linOct[i].getZ())) +0.5)/((double)(1u<<maxDepth));
  }//end for i
  linOct.clear();
  gSize[0] = 1.;
  gSize[1] = 1.;
  gSize[2] = 1.;

  PetscTruth usingRegularOctree;
  PetscInt regLev;
  PetscOptionsHasName(0,"-useRegularOctreeAtLevel",&usingRegularOctree);
  PetscOptionsGetInt(0,"-useRegularOctreeAtLevel",&regLev,0);

  if(usingRegularOctree) {
    if(!rank) {
      std::cout<<"Creating Regular Fine Grid Octree."<<std::endl;
    }
    balOct.clear();
    createRegularOctree(balOct,regLev,3,maxDepth,MPI_COMM_WORLD);
  } else {   
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
  }

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

  ot::DAMG        *damg;    
  int       nlevels = 1; //number of multigrid levels
  PetscInt       numRefinements = 0;
  unsigned int       dof =1;// degrees of freedom per node  

  PetscInt nlevelsPetscInt = nlevels;
  PetscOptionsGetInt(0, "-nlevels", &nlevelsPetscInt, 0);
  nlevels = nlevelsPetscInt;

  PetscOptionsGetInt(0,"-numRefinements",&numRefinements,0);
  for(int i = 0; i < numRefinements; i++) {
    std::vector<ot::TreeNode> tmpOct = balOct;
    balOct.clear();
    ot::refineOctree(tmpOct, balOct); 
  }

  MPI_Barrier(MPI_COMM_WORLD);	

  if(!rank) {
    std::cout<<"nlevels initial: "<<nlevels<<std::endl;
  }

  // Note: The user context for all levels will be set separately later.
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

  unsigned int* allSizes = NULL;

  if(!rank) {
    allSizes = new unsigned int[6*size]; 
  }

  for(int lev = 0; lev < nlevels; lev++) {
    for(int auxCtr = 0; auxCtr < 2; auxCtr++) {
      ot::DA* currDa = NULL;

      if(auxCtr == 0) {
        currDa = damg[lev]->da;
      }else if(damg[lev]->da_aux) {
        currDa = damg[lev]->da_aux;
      } else {
        break;
      }

      unsigned int localSizes[6];

      //Pre-ghost size
      localSizes[0] = currDa->getIdxElementBegin();

      //Own element size
      localSizes[2] = currDa->getElementSize();

      unsigned int localTotalSize = currDa->getLocalBufferSize();

      //Post-ghost size
      localSizes[1] = (localTotalSize - (currDa->getIdxPostGhostBegin()));

      //Independent size
      localSizes[3] = currDa->getIndependentSize();

      //PreGhostElementSize
      localSizes[4] = currDa->getPreGhostElementSize();

      //Own Boundary size
      localSizes[5] = ((currDa->getIdxPostGhostBegin()) - (currDa->getIdxElementEnd()));

      MPI_Gather(localSizes, 6, MPI_UNSIGNED, allSizes, 6,
          MPI_UNSIGNED, 0, MPI_COMM_WORLD);

      MPI_Barrier(MPI_COMM_WORLD);

      /*
      if(!rank) {
        std::cout<<"Printing Info for lev: "<<lev<<" auxCtr: "<<auxCtr<<std::endl; 
        for(int i = 0; i < size; i++) {
          std::cout<<"Proc: "<<i<<" PreGh: "<<allSizes[6*i]
            <<" PostGh: "<<allSizes[(6*i)+1]<<" own: "
            <<allSizes[(6*i)+2]<<" Ind: "<<allSizes[(6*i)+3]
            <<" preGhElems: "<<allSizes[(6*i)+4]<<" LocalBnd: "
            <<allSizes[(6*i)+5]<<std::endl;
          fflush(stdout);
        }
      }
      */
      fflush(stdout);
    }//end da or da_aux
  }//end lev

  if(!rank) {
    delete [] allSizes;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if(setMatPropsUsingPts) {
    PetscTruth setMatPropsAtCoarsest;
    PetscOptionsHasName(0,"-setMatPropsAtCoarsest",&setMatPropsAtCoarsest);
    if(setMatPropsAtCoarsest) {
      //Coarsest is only 1 processor. So to align with coarse blocks, everything
      //must be sent to p0
      if(!rank) {
        par::scatterValues<ot::TreeNode>(matPropNodes, linOct,
            totalNumPts, MPI_COMM_WORLD);
      } else {
        par::scatterValues<ot::TreeNode>(matPropNodes, linOct,
            0, MPI_COMM_WORLD);
      }
      matPropNodes.clear();

      std::vector<double> matPropPts(3*(linOct.size()));
      for(int i=0;i<linOct.size();i++) {
        matPropPts[3*i] =
          (((double)(linOct[i].getX())) + 0.5)/((double)(1u << maxDepth));
        matPropPts[(3*i)+1] =
          (((double)(linOct[i].getY())) +0.5)/((double)(1u << maxDepth));
        matPropPts[(3*i)+2] =
          (((double)(linOct[i].getZ())) +0.5)/((double)(1u << maxDepth));
      }//end for i
      linOct.clear();

      PetscReal lapFac = 0.0;
      PetscTruth optFound;
      PetscOptionsGetReal("lap","-MatPropFac",&lapFac,&optFound);
      std::vector<double> lapJumps((matPropPts.size()/3));
      for(int i=0;i<lapJumps.size();i++) {
        lapJumps[i] = lapFac;
      }
      SetCoarseToFineFromPts(damg, matPropPts, lapJumps);
    } else {
      if(damg[nlevels-1]->da->iAmActive()) {
        int npesActive;
        int rankActive;

        MPI_Comm_size(damg[nlevels-1]->da->getCommActive(), &npesActive);        
        MPI_Comm_rank(damg[nlevels-1]->da->getCommActive(), &rankActive);        

        unsigned int avgSize = (totalNumPts/npesActive);
        unsigned int remSize = (totalNumPts % npesActive);

        if(rankActive < remSize) {
          par::scatterValues<ot::TreeNode>(matPropNodes, linOct,
              (avgSize + 1), MPI_COMM_WORLD);
        } else {
          par::scatterValues<ot::TreeNode>(matPropNodes, linOct,
              avgSize, MPI_COMM_WORLD);
        }

        matPropNodes.clear();
        for(unsigned int i = 0; i < linOct.size(); i++) {
          unsigned int xpt, ypt, zpt;
          linOct[i].getAnchor(xpt,ypt,zpt);
          ot::TreeNode tmpNode(xpt,ypt,zpt,(maxDepth+1),3,(maxDepth+1));
          matPropNodes.push_back(tmpNode); 
        }
        linOct.clear();

        std::vector<ot::TreeNode> minsArray(npesActive);

        std::vector<ot::TreeNode> blocks = damg[nlevels-1]->da->getBlocks();

        assert(!blocks.empty());

        ot::TreeNode firstBlock = blocks[0];

        par::Mpi_Allgather<ot::TreeNode>(&firstBlock, &(*(minsArray.begin())), 1,
            damg[nlevels-1]->da->getCommActive());

        int *sendCnt = new int[npesActive];
        int *recvCnt = new int[npesActive];
        int *sendOffsets = new int[npesActive];
        int *recvOffsets = new int[npesActive];

        //compute the total number of nodes being sent to each proc ...
        for (int i = 0; i < npesActive; i++) {
          sendCnt[i]=0;
          recvCnt[i]=0;
        }

        unsigned int blockCtr = 0;
        for(unsigned int tempCtr = 0; tempCtr < matPropNodes.size(); tempCtr++) {
          while( ((blockCtr + 1) < npesActive) &&
              (matPropNodes[tempCtr] >= minsArray[blockCtr+1]) ) {
            blockCtr++;
          }
          if ( (blockCtr + 1) == npesActive) {
            sendCnt[blockCtr] += (matPropNodes.size() - tempCtr);
            break;
          } else {
            sendCnt[blockCtr]++;
          }
        }

        minsArray.clear();

        // communicate with other procs how many you shall be sending and get how
        // many to recieve from whom.
        par::Mpi_Alltoall<int>(sendCnt, recvCnt, 1, damg[nlevels-1]->da->getCommActive());

        unsigned int totalRecv = 0;
        for (unsigned int i = 0; i < npesActive; i++) {
          totalRecv += recvCnt[i];
        }//end for i

        linOct.resize(totalRecv);

        sendOffsets[0] = 0;
        recvOffsets[0] = 0;

        // compute offsets ...
        for (int i=1; i < npesActive; i++) {
          sendOffsets[i] = sendOffsets[i-1] + sendCnt[i-1];
          recvOffsets[i] = recvOffsets[i-1] + recvCnt[i-1];
        }//end for i

        // perform All2Allv
        par::Mpi_Alltoallv_sparse<ot::TreeNode>(&(*(matPropNodes.begin())), sendCnt, sendOffsets,        
            &(*(linOct.begin())), recvCnt, recvOffsets, damg[nlevels-1]->da->getCommActive());

        matPropNodes.clear();

        // clean up ...
        delete [] sendCnt;
        sendCnt = NULL;

        delete [] recvCnt;
        recvCnt = NULL;

        delete [] sendOffsets;
        sendOffsets = NULL;

        delete [] recvOffsets;
        recvOffsets = NULL;

      } else {
        par::scatterValues<ot::TreeNode>(matPropNodes, linOct,
            0, MPI_COMM_WORLD);
        matPropNodes.clear();
      }

      std::vector<double> matPropPts(3*(linOct.size()));
      for(int i=0;i<linOct.size();i++) {
        matPropPts[3*i] =
          (((double)(linOct[i].getX())) + 0.5)/((double)(1u << maxDepth));
        matPropPts[(3*i)+1] =
          (((double)(linOct[i].getY())) +0.5)/((double)(1u << maxDepth));
        matPropPts[(3*i)+2] =
          (((double)(linOct[i].getZ())) +0.5)/((double)(1u << maxDepth));
      }//end for i
      linOct.clear();

      PetscReal lapFac = 0.0;
      PetscTruth optFound;
      PetscOptionsGetReal("lap","-MatPropFac",&lapFac,&optFound);
      std::vector<double> lapJumps((matPropPts.size()/3));
      for(int i=0;i<lapJumps.size();i++) {
        lapJumps[i] = lapFac;
      }
      SetUserContextsFromPts(damg, matPropPts, lapJumps); 
    }
  } else {
    SetUserContexts(damg);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if(!rank) {
    std::cout << "Set User Contexts all levels."<< std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  PetscInt       jacType = 1;
  PetscOptionsGetInt(0,"-jacType",&jacType,0);

  MPI_Barrier(MPI_COMM_WORLD);
  if(!rank) {
    std::cout << "Creating stencils..."<< std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  createLmatType2(LaplacianType2Stencil);

  createMmatType2(MassType2Stencil);

  if(jacType == 3) {
    createLmatType1(LaplacianType1Stencil);
    createMmatType1(MassType1Stencil);
  }
  createShFnMat(ShapeFnStencil);

  MPI_Barrier(MPI_COMM_WORLD);
  if(!rank) {
    std::cout << "Created all stencils."<< std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  PetscInt rhsType = 2;
  PetscOptionsGetInt(0,"-rhsType",&rhsType,0);

  //Function handles
  PetscErrorCode (*ComputeRHSHandle)(ot::DAMG damg,Vec rhs) = NULL;
  PetscErrorCode (*CreateJacobianHandle)(ot::DAMG damg,Mat *B) = NULL;
  PetscErrorCode (*ComputeJacobianHandle)(ot::DAMG damg,Mat J, Mat B) = NULL;

  if(rhsType == 0) {
    ComputeRHSHandle = ComputeRHS0;
  } else if (rhsType == 1) {
    ComputeRHSHandle = ComputeRHS1;
  } else if (rhsType == 2) {
    ComputeRHSHandle = ComputeRHS2;
  } else if (rhsType == 3) {
    ComputeRHSHandle = ComputeRHS3;
  } else if (rhsType == 4) {
    ComputeRHSHandle = ComputeRHS4;
  } else if (rhsType == 5) {
    ComputeRHSHandle = ComputeRHS5;
  } else if (rhsType == 6) {
    ComputeRHSHandle = ComputeRHS6;
  } else if (rhsType == 7) {
    ComputeRHSHandle = ComputeRHS7;
  } else if (rhsType == 8) {
    ComputeRHSHandle = ComputeRHS8;
  } else {
    assert(false);
  }

  if(jacType == 1) {
    CreateJacobianHandle = CreateJacobian1;
    ComputeJacobianHandle = ComputeJacobian1;
  } else if (jacType == 2) {
    CreateJacobianHandle = CreateJacobian2;
    ComputeJacobianHandle = ComputeJacobian2;
  } else if (jacType == 3) {
    CreateJacobianHandle = CreateJacobian3;
    ComputeJacobianHandle = ComputeJacobian3;
    //Skip the finest and the coarsest levels. For the other levels, J and B
    //must be different
    for(int i = 1; i < (nlevels-1); i++) {
      ot::DAMGCreateJMatrix(damg[i], CreateJacobianHandle);
    }
  } else {
    assert(false);
  }

  //Global Function Handles for using KSP_Shell (will be used @ the coarsest grid if not all
  //processors are active on the coarsest grid)
  if (jacType == 1) {
    ot::getPrivateMatricesForKSP_Shell = getPrivateMatricesForKSP_Shell_Jac1;
  } else if (jacType == 2) {
    ot::getPrivateMatricesForKSP_Shell = getPrivateMatricesForKSP_Shell_Jac2;
  } else if (jacType == 3) {
    ot::getPrivateMatricesForKSP_Shell = getPrivateMatricesForKSP_Shell_Jac3;
  } else {
    assert(false);
  }

  ot::DAMGSetKSP(damg, CreateJacobianHandle, ComputeJacobianHandle, ComputeRHSHandle);

  MPI_Barrier(MPI_COMM_WORLD);
  if(!rank) {
    std::cout << "Called DAMGSetKSP"<< std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  PetscReal norm2;
  PetscReal normInf;

  PetscTruth setRandomGuess = PETSC_FALSE;
  PetscOptionsHasName(0,"-setRandomGuess",&setRandomGuess);

  if(setRandomGuess) { 
    ot::DAMGSetInitialGuess(damg,ot::DAMGInitialGuessCurrent);

    PetscRandom rctx;  
    PetscRandomCreate(PETSC_COMM_WORLD,&rctx);
    PetscRandomSetType(rctx,PETSCRAND48);
    PetscInt randomSeed = 12345;
    PetscOptionsGetInt(0,"-randomSeed",&randomSeed,0);
    if(!rank) {
      std::cout<<"Using Random Seed: "<<randomSeed<<std::endl;
    }
    PetscRandomSetSeed(rctx,randomSeed);
    PetscRandomSeed(rctx);
    PetscRandomSetFromOptions(rctx);

    VecSetRandom((DAMGGetx(damg)),rctx);

    PetscRandomDestroy(rctx);

    VecNorm((DAMGGetx(damg)),NORM_INFINITY,&normInf);
    VecNorm((DAMGGetx(damg)),NORM_2,&norm2);

    MPI_Barrier(MPI_COMM_WORLD);
    if(!rank) {
      std::cout << "Solving with random intial guess"<< std::endl;
      std::cout<<"Initial Norm2: "<<norm2<<" NormInf: "<<normInf<<std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }else {
    MPI_Barrier(MPI_COMM_WORLD);
    if(!rank) {
      std::cout << "Solving with 0-intial guess"<< std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  startTime = MPI_Wtime();
  ot::DAMGSolve(damg);
  endTime = MPI_Wtime();

  localTime = (endTime - startTime);
  par::Mpi_Reduce<double>(&localTime, &totalTime, 1, MPI_MAX, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

  if (!rank) {
    std::cout << GRN << "Done Solve" << NRM << std::endl;
    std::cout << "Solve Time: "<<totalTime << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  destroyLmatType2(LaplacianType2Stencil);
  destroyMmatType2(MassType2Stencil);
  if(jacType == 3) {
    destroyLmatType1(LaplacianType1Stencil);
    destroyMmatType1(MassType1Stencil);
  }
  destroyShFnMat(ShapeFnStencil);

  if (!rank) {
    std::cout << GRN << "Destroyed Stencils" << NRM << std::endl;
  }

  DestroyUserContexts(damg);

  if (!rank) {
    std::cout << GRN << "Destroyed User Contexts." << NRM << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  DAMGDestroy(damg);

  if (!rank) {
    std::cout << GRN << "Destroyed DAMG" << NRM << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

#ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif
  balOct.clear();

  ot::DAMG_Finalize();
  if (!rank) {
    std::cout << GRN << "Finalizing ..." << NRM << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  PetscFinalize();
}//end function

