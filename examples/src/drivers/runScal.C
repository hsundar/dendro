
#include "mpi.h"
#include "petsc.h"
#include "sys.h"
#include "octUtils.h"
#include "TreeNode.h"
#include "parUtils.h"
#include <cstdlib>
#include "externVars.h"
#include "dendro.h"

#ifdef MPI_WTIME_IS_GLOBAL
#undef MPI_WTIME_IS_GLOBAL
#endif

int main(int argc, char ** argv ) {	
  int size, rank;
  double startTime, endTime, minTime;
  unsigned int incInt = 1;
  bool incCorner = 1;  
  unsigned int numPts;
  unsigned int writePOut = 0;
  unsigned int writeBOut = 0;
  unsigned int readPtsFile = 0;
  unsigned int readOctFile =0;
  unsigned int runOpt =2;//2 for both, 1 for pt alone and 0 for bal alone.
  char Kstr[20];
  char inpFileName[50], p2nIn[50], p2nOut[50], n2oOut[50], balOut[50];
  double gSize[3];

  PetscInitialize(&argc,&argv,"options",NULL);
  ot::RegisterEvents();

#ifdef PETSC_USE_LOG
  int stages[4];
  PetscLogStageRegister(&stages[0],"Prepare Input.");
  PetscLogStageRegister(&stages[1],"P2O");
  PetscLogStageRegister(&stages[2],"N2O");
  PetscLogStageRegister(&stages[3],"Bal");
#endif

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(argc < 3) {
    std::cerr << "Usage: " << argv[0] << " inpfile complx incInt[1] runOpt[2]"<<
      "readPtsFile[0] readOctFile[0] maxDepth[30] writePOut[0] writeBOut[0]"<<
      " dim[3] maxNumPtsPerOctant[1] incCorner[0] "  << std::endl;
    return -1;
  }
  unsigned int maxNumPts= 1;
  int complx = atoi(argv[2]);	
  unsigned int dim=3;
  unsigned int maxDepth=30;
  if(argc > 3) {
    incInt = atoi(argv[3]);
  }
  if(argc > 4) {
    runOpt = atoi(argv[4]);
  }
  if(argc > 5) {
    readPtsFile = atoi(argv[5]);
  }
  if(argc > 6) {
    readOctFile = atoi(argv[6]);
  }
  if(argc > 7) {
    maxDepth = atoi(argv[7]);
  }
  if(argc > 8) {
    writePOut = atoi(argv[8]);
  }
  if(argc > 9) {
    writeBOut = atoi(argv[9]);
  }
  if(argc > 10) {
    dim = atoi(argv[10]);
  }
  if(argc > 11) {
    maxNumPts = atoi(argv[11]);
  }
  if(argc > 12) { incCorner = (bool)(atoi(argv[12]));}

  strcpy(inpFileName, argv[1]);
  strcpy(balOut,inpFileName);
  ot::int2str(rank,Kstr);
  strcat(balOut,Kstr);
  strcat(balOut,"_\0");
  ot::int2str(size,Kstr);
  strcat(balOut,Kstr);
  strcpy(n2oOut,balOut);
  strcpy(p2nIn,balOut);
  strcpy(p2nOut,balOut);
  strcat(balOut,"_Bal.ot\0");
  strcat(n2oOut,"_Con.ot\0");
  strcat(p2nIn,".pts\0");
  strcat(p2nOut,"_Out.pts\0");
  strcat(inpFileName,".inp\0");	
  std::vector<ot::TreeNode> nodes;
  std::vector<double> pts;
  unsigned int ptsLen;

  // print out the format ...
  /*
  if (!rank) {
    cout << "============================================================================================" << endl;
    cout << "#Points p2n_time p2n_imbal #cellsConstr #cellsBal bal_time bal_imbal maxLevBef minLevBef avgLevBef maxLevAft minLevAft avgLevAft" << endl;
    cout << "============================================================================================" << endl << endl;
  }
  */
#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[0]);
#endif
  if(runOpt == 1 || runOpt == 2)  {
    if(readPtsFile) {
      if(!rank){
        std::cout << " reading  "<<p2nIn<<std::endl; // Point size
      }
      ot::readPtsFromFile(p2nIn, pts);
      if(!rank){
        std::cout << " finished reading  "<<p2nIn<<std::endl; // Point size
      }
      ptsLen = pts.size();
      std::vector<ot::TreeNode> tmpNodes;
      for(int i = 0; i < ptsLen; i += 3) {
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
      }//end for i
      pts.clear();
    //  std::cout<<"Before Removing Duplicates, the number of pts: "<<tmpNodes.size()<<std::endl; 
      par::removeDuplicates<ot::TreeNode>(tmpNodes,false,MPI_COMM_WORLD);	
     // std::cout<<"After Removing Duplicates, the number of pts: "<<tmpNodes.size()<<std::endl; 
      nodes = tmpNodes;
      tmpNodes.clear();
      par::partitionW<ot::TreeNode>(nodes, NULL,MPI_COMM_WORLD);
      // reduce and only print the total ...
      DendroIntL locSz, totalSz;
      locSz = nodes.size();

      par::Mpi_Reduce<DendroIntL>(&locSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
      
      if(rank==0) {
      std::cout<<"rank= "<<rank<<" size= " << size << " total pts= " << totalSz<<std::endl;
      }

      pts.resize(3*(nodes.size()));
      ptsLen = (3*(nodes.size()));
      for(int i=0;i<nodes.size();i++) {
        pts[3*i] = (((double)(nodes[i].getX())) + 0.5)/((double)(1u<<maxDepth));
        pts[(3*i)+1] = (((double)(nodes[i].getY())) +0.5)/((double)(1u<<maxDepth));
        pts[(3*i)+2] = (((double)(nodes[i].getZ())) +0.5)/((double)(1u<<maxDepth));
      }//end for i
      nodes.clear();
    }//end if readPts
    gSize[0] = 1.0;
    gSize[1] = 1.0;
   gSize[2] = 1.0;
    if(writePOut) { 
      if(!rank) {
        std::cout<<"Writing pts to: "<<p2nOut<<std::endl;
      }
      ot::writePtsToFile(p2nOut,pts);
    }
 #ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif
    MPI_Barrier(MPI_COMM_WORLD);	
#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[1]);
#endif
    startTime = MPI_Wtime();
    ot::points2Octree(pts, gSize, nodes, dim, maxDepth, maxNumPts, MPI_COMM_WORLD);
    endTime = MPI_Wtime();
#ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif
    double locTime, totalTime, minTime;
    locTime = endTime - startTime;
    par::Mpi_Reduce<double>(&locTime, &totalTime, 1, MPI_MAX, 0, MPI_COMM_WORLD);
    par::Mpi_Reduce<double>(&locTime, &minTime, 1, MPI_MIN, 0, MPI_COMM_WORLD);
    if(!rank){
      std::cout <<"P2n Time: "<<totalTime << " " << "secs imbalance: "<< totalTime/minTime << std::endl;
    }
    pts.clear();
  }
  //TestBuildOCt from here on:
  if(runOpt == 0 || runOpt == 2) {
    if(readOctFile) {
      ot::readNodesFromFile(n2oOut,nodes);
    }
    std::vector<ot::TreeNode > linOct;	
#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[2]);
#endif
    DendroIntL inputNodes = nodes.size();
    DendroIntL totInp;
    par::Mpi_Reduce<DendroIntL>(&inputNodes, &totInp, 1, MPI_SUM, 0, MPI_COMM_WORLD);
    if(!rank) {
      std::cout << "Nodes after p2o: "<< totInp  << std::endl;       
    }
    ot::completeOctree(nodes,linOct,dim,maxDepth,false, false, false, MPI_COMM_WORLD);
#ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif

    assert(!linOct.empty());
    if(writeBOut) { 
      if(!rank) {
        std::cout<<"Writing octree after p2o to: "<<n2oOut<<std::endl; 
      }
      ot::writeNodesToFile(n2oOut,linOct);
    }
    inputNodes = linOct.size();
    nodes.clear();
    std::vector<ot::TreeNode > balOct;

    par::Mpi_Reduce<DendroIntL>(&inputNodes, &totInp, 1, MPI_SUM, 0, MPI_COMM_WORLD);
    if(!rank) {
      std::cout << "Nodes before balancing: "<< totInp  << std::endl;       
    }
    MPI_Barrier(MPI_COMM_WORLD);	
#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[3]);
#endif
    startTime = MPI_Wtime();

    ot::balanceOctree (linOct, balOct, dim, maxDepth, incCorner, MPI_COMM_WORLD, NULL, NULL);
    endTime = MPI_Wtime();
#ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif
    assert(!balOct.empty());
    if(writeBOut) { 
      if(!rank) {
        std::cout<<"Writing octree after balancing to: "<<balOut<<std::endl; 
      }
      ot::writeNodesToFile(balOut,balOct);
    }
    // compute total inp size and output size
    DendroIntL locBalNodes = balOct.size();
    DendroIntL totBal;
    double balTime, totTime;
    balTime = endTime - startTime;
    par::Mpi_Reduce<DendroIntL>(&locBalNodes, &totBal, 1, MPI_SUM, 0, MPI_COMM_WORLD);
    par::Mpi_Reduce<double>(&balTime, &totTime, 1, MPI_MAX, 0, MPI_COMM_WORLD);
    par::Mpi_Reduce<double>(&balTime, &minTime, 1, MPI_MIN, 0, MPI_COMM_WORLD);

    if(!rank) {
      std::cout << "Nodes after balancing: "<< totBal << std::endl;       
      std::cout << "bal Time: "<<totTime << " " << "secs imbalance: "<< totTime/minTime << std::endl;
    }
    linOct.clear();
    balOct.clear();
  }

  PetscFinalize();

}

