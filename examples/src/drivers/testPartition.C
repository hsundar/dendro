#include "mpi.h"
#include "petsc.h"
#include "sys.h"

#include "parUtils.h"
#include "TreeNode.h"
#include "colors.h"
#include "oda.h"

#include <cstdlib>

#include "externVars.h"
#include "dendro.h"
#include "rotation.h"
#include "treenode2vtk.h"

//Don't want time to be synchronized. Need to check load imbalance.
#ifdef MPI_WTIME_IS_GLOBAL
#undef MPI_WTIME_IS_GLOBAL
#endif

int main(int argc, char **argv) { 
  int size, rank;
  bool incCorner = 1;
  // char bFile[50];
  char ptsFileName[256];
  bool compressLut = false;

  unsigned int ptsLen;
  unsigned int maxNumPts = 1;
  unsigned int dim = 3;
  unsigned int maxDepth = 4;
  double gSize[3];
  //initializeHilbetTable(2);
  initializeHilbetTable(3);
  
  double localTime, totalTime;
  double startTime, endTime;
  DendroIntL localSz, totalSz;
  std::vector<double> pts;
  std::vector<ot::TreeNode> linOct, balOct;

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " inpfile " << std::endl;
    return -1;
  }

  PetscInitialize(&argc, &argv, "options.hs", NULL);
  ot::RegisterEvents();
  ot::DA_Initialize(MPI_COMM_WORLD);

#ifdef PETSC_USE_LOG
  int stages[3];
  PetscLogStageRegister("P2O.", &stages[0]);
  PetscLogStageRegister("Bal", &stages[1]);
  PetscLogStageRegister("ODACreate", &stages[2]);
#endif

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // if(argc > 2) { compressLut = (bool)(atoi(argv[2]));}
  // sprintf(bFile, "%s%d_%d.ot", argv[1], rank, size);
  //=============================================================
  sprintf(ptsFileName, "%s%d_%d.pts", argv[1], rank, size);

  //Read pts from files
  MPI_Barrier(MPI_COMM_WORLD);
  if (!rank) {
    std::cout << " reading  " << ptsFileName << std::endl; // Point size
  }
  ot::readPtsFromFile(ptsFileName, pts);
  if (!rank) {
    std::cout << " finished reading  " << ptsFileName << std::endl; // Point size
  }
  
  
  
  MPI_Barrier(MPI_COMM_WORLD);

  
  
  ptsLen = pts.size();
    // @milinda debug code
  //std::cout << rank << ": read ptsLen" << ptsLen << " points" << std::endl;
  std::cout << rank << ": read pts.size()" << pts.size() << " points" << std::endl;
  
  std::vector<ot::TreeNode> tmpNodes;
  for (int i = 0; i < ptsLen; i += 3) {
   
    
        
    if ((pts[i] > 0.0) &&
        (pts[i + 1] > 0.0)
        && (pts[i + 2] > 0.0) &&
        (((unsigned int)(pts[i] * ((double)(1u << maxDepth)))) < (1u << maxDepth)) &&
        (((unsigned int)(pts[i + 1] * ((double)(1u << maxDepth)))) < (1u << maxDepth)) &&
        (((unsigned int)(pts[i + 2] * ((double)(1u << maxDepth)))) < (1u << maxDepth))) {
      // @milinda debug code
      //std::cout<<"Rank:"<<rank<<" x,y,z :"<<pts[i]<<","<<pts[i+1]<<","<<pts[i+2]<<std::endl;
      
      tmpNodes.push_back(ot::TreeNode((unsigned int)(pts[i] * (double)(1u << maxDepth)),
                                      (unsigned int)(pts[i + 1] * (double)(1u << maxDepth)),
                                      (unsigned int)(pts[i + 2] * (double)(1u << maxDepth)),
                                      maxDepth, dim, maxDepth));
    }
  }
  pts.clear();
  
//   std::cout << rank << ": read " << tmpNodes.size() << " points" << std::endl;
//   treeNodesTovtk(tmpNodes,rank,"vtkTreeNode");
//   
  par::removeDuplicates<ot::TreeNode>(tmpNodes, false, MPI_COMM_WORLD);
  
  std::cout << rank << "afterRemoveDuplicates: " << tmpNodes.size() << " points" << std::endl;
  linOct = tmpNodes;
  tmpNodes.clear();
  par::partitionW<ot::TreeNode>(linOct, NULL, MPI_COMM_WORLD);
  // reduce and only print the total ...
  localSz = linOct.size();
  par::Mpi_Reduce<DendroIntL>(&localSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    std::cout << "# pts= " << totalSz << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  pts.resize(3 * (linOct.size()));
  ptsLen = (3 * (linOct.size()));
  for (int i = 0; i < linOct.size(); i++) {
    pts[3 * i] = (((double)(linOct[i].getX())) + 0.5) / ((double)(1u << maxDepth));
    pts[(3 * i) + 1] = (((double)(linOct[i].getY())) + 0.5) / ((double)(1u << maxDepth));
    pts[(3 * i) + 2] = (((double)(linOct[i].getZ())) + 0.5) / ((double)(1u << maxDepth));
  } //end for i
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
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  treeNodesTovtk(linOct,rank,"bfBalancing");
  
  localTime = endTime - startTime;
  par::Mpi_Reduce<double>(&localTime, &totalTime, 1, MPI_MAX, 0, MPI_COMM_WORLD);
  if (!rank) {
    std::cout << "P2n Time: " << totalTime << std::endl;
  }
  // reduce and only print the total ...
  localSz = linOct.size();
  par::Mpi_Reduce<DendroIntL>(&localSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    std::cout << "# of Unbalanced Octants: " << totalSz << std::endl;
  }
  pts.clear();

  //Balancing...
  MPI_Barrier(MPI_COMM_WORLD);
#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[1]);
#endif
  startTime = MPI_Wtime();
  ot::balanceOctree(linOct, balOct, dim, maxDepth, incCorner, MPI_COMM_WORLD, NULL, NULL);
  endTime = MPI_Wtime();
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  treeNodesTovtk(balOct,rank,"afBalancing");
  
#ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif
  linOct.clear();
  // compute total inp size and output size
  localSz = balOct.size();
  localTime = endTime - startTime;
  par::Mpi_Reduce<DendroIntL>(&localSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  par::Mpi_Reduce<double>(&localTime, &totalTime, 1, MPI_MAX, 0, MPI_COMM_WORLD);

  if (!rank) {
    std::cout << "# of Balanced Octants: " << totalSz << std::endl;
    std::cout << "bal Time: " << totalTime << std::endl;
  }



  //=============================================================
  //ODA ...
  MPI_Barrier(MPI_COMM_WORLD);
#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[2]);
#endif
  startTime = MPI_Wtime();
  assert(!(balOct.empty()));
  ot::DA da(balOct, MPI_COMM_WORLD, MPI_COMM_WORLD, compressLut);
  endTime = MPI_Wtime();
#ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif
  balOct.clear();
  // compute total inp size and output size
  localSz = da.getNodeSize();
  localTime = endTime - startTime;
  par::Mpi_Reduce<DendroIntL>(&localSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  par::Mpi_Reduce<double>(&localTime, &totalTime, 1, MPI_MAX, 0, MPI_COMM_WORLD);

  if (!rank) {
    std::cout << "Total # Vertices: " << totalSz << std::endl;
    std::cout << "Time to build ODA: " << totalTime << std::endl;
  }

  //! Quality of the partition ...
  DendroIntL maxNodeSize, minNodeSize,
     maxBdyNode, minBdyNode,
     maxIndepSize, minIndepSize,
     maxElementSize, minElementSize;

  localSz = da.getNodeSize();
  par::Mpi_Reduce<DendroIntL>(&localSz, &maxNodeSize, 1, MPI_MAX, 0, MPI_COMM_WORLD);
  par::Mpi_Reduce<DendroIntL>(&localSz, &minNodeSize, 1, MPI_MIN, 0, MPI_COMM_WORLD);

  localSz = da.getBoundaryNodeSize();
  par::Mpi_Reduce<DendroIntL>(&localSz, &maxBdyNode, 1, MPI_MAX, 0, MPI_COMM_WORLD);
  par::Mpi_Reduce<DendroIntL>(&localSz, &minBdyNode, 1, MPI_MIN, 0, MPI_COMM_WORLD);

  localSz = da.getElementSize();
  par::Mpi_Reduce<DendroIntL>(&localSz, &maxElementSize, 1, MPI_MAX, 0, MPI_COMM_WORLD);
  par::Mpi_Reduce<DendroIntL>(&localSz, &minElementSize, 1, MPI_MIN, 0, MPI_COMM_WORLD);

  localSz = da.getIndependentSize();
  par::Mpi_Reduce<DendroIntL>(&localSz, &maxIndepSize, 1, MPI_MAX, 0, MPI_COMM_WORLD);
  par::Mpi_Reduce<DendroIntL>(&localSz, &minIndepSize, 1, MPI_MIN, 0, MPI_COMM_WORLD);

  if (!rank) {
    std::cout << "Nodes          \t(" << minNodeSize << ", " << maxNodeSize << ")" << std::endl;
    std::cout << "Boundary Node  \t(" << minBdyNode << ", " << maxBdyNode << ")" << std::endl;
    std::cout << "Element        \t(" << minElementSize << ", " << maxElementSize << ")" << std::endl;
    std::cout << "Independent    \t(" << minIndepSize << ", " << maxIndepSize << ")" << std::endl;
  }

  //! ========================

  if (!rank) {
    std::cout << GRN << "Finalizing ..." << NRM << std::endl;
  }

  ot::DA_Finalize();
  PetscFinalize();
} //end function

