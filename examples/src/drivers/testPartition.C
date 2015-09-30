#include "mpi.h"
#include "petsc.h"
#include "sys.h"

#include "parUtils.h"
#include "TreeNode.h"
#include "colors.h"
#include "oda.h"

#include <cstdlib>
#include <execinfo.h>
#include <unistd.h>

#include "externVars.h"
#include "dendro.h"
#include "rotation.h"
#include "treenode2vtk.h"

//Don't want time to be synchronized. Need to check load imbalance.
#ifdef MPI_WTIME_IS_GLOBAL
#undef MPI_WTIME_IS_GLOBAL
#endif


void handler(int sig) {
  void *array[10];
  int size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}


int main(int argc, char **argv) {
  signal(SIGSEGV, handler);   // install our handler

  int size, rank;
  bool incCorner = 1;
  // char bFile[50];
  char ptsFileName[256];
  bool compressLut = false;

  unsigned int ptsLen;
  unsigned int maxNumPts = 1;
  unsigned int dim = 3;
  unsigned int maxDepth = 8;
  double gSize[3];
  //initializeHilbetTable(2);
  G_MAX_DEPTH=maxDepth;
  G_dim=dim;
  initializeHilbetTable(dim);

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

  sprintf(ptsFileName, "%s%d_%d.pts", argv[1], rank, size);

  //Read pts from files
  if (!rank) {
    std::cout << RED " Reading  " << argv[1] << NRM << std::endl; // Point size
  }
  ot::readPtsFromFile(ptsFileName, pts);

  if (!rank) {
    std::cout << GRN " Finished reading  " << argv[1] << NRM << std::endl; // Point size
  }

  ptsLen = pts.size();

  std::vector<ot::TreeNode> tmpNodes;
  for (int i = 0; i < ptsLen; i += 3) {
    if ((pts[i] > 0.0) &&
        (pts[i + 1] > 0.0)
        && (pts[i + 2] > 0.0) &&
        (((unsigned int) (pts[i] * ((double) (1u << maxDepth)))) < (1u << maxDepth)) &&
        (((unsigned int) (pts[i + 1] * ((double) (1u << maxDepth)))) < (1u << maxDepth)) &&
        (((unsigned int) (pts[i + 2] * ((double) (1u << maxDepth)))) < (1u << maxDepth))) {
      tmpNodes.push_back(ot::TreeNode((unsigned int) (pts[i] * (double) (1u << maxDepth)),
                                      (unsigned int) (pts[i + 1] * (double) (1u << maxDepth)),
                                      (unsigned int) (pts[i + 2] * (double) (1u << maxDepth)),
                                      maxDepth, dim, maxDepth));
    }
  }
  pts.clear();

  // treeNodesTovtk(tmpNodes, rank, "vtkTreeNode");

  // std::cout << rank << "removeDuplicates" << std::endl;
  par::removeDuplicates<ot::TreeNode>(tmpNodes, false, MPI_COMM_WORLD);

  linOct = tmpNodes;
  tmpNodes.clear();
  // std::cout << rank << "partition" << std::endl;
  par::partitionW<ot::TreeNode>(linOct, NULL, MPI_COMM_WORLD);

  //treeNodesTovtk(linOct,rank,"par_1");

  // reduce and only print the total ...
  localSz = linOct.size();
  par::Mpi_Reduce<DendroIntL>(&localSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);

  pts.resize(3 * (linOct.size()));
  ptsLen = (3 * (linOct.size()));
  for (int i = 0; i < linOct.size(); i++) {
    pts[3 * i] = (((double) (linOct[i].getX())) + 0.5) / ((double) (1u << maxDepth));
    pts[(3 * i) + 1] = (((double) (linOct[i].getY())) + 0.5) / ((double) (1u << maxDepth));
    pts[(3 * i) + 2] = (((double) (linOct[i].getZ())) + 0.5) / ((double) (1u << maxDepth));
  } //end for i
  linOct.clear();
  gSize[0] = 1.0;
  gSize[1] = 1.0;
  gSize[2] = 1.0;

  // ========== Points2Octree ===========
  MPI_Barrier(MPI_COMM_WORLD);
  if (!rank) {
    std::cout << BLU << "===============================================" << NRM << std::endl;
    std::cout << RED " Starting Points to Octree" NRM << std::endl;
    std::cout << BLU << "===============================================" << NRM << std::endl;
  }

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
    std::cout << GRN " P2n Time: " YLW << totalTime << NRM << std::endl;
  }
  // reduce and only print the total ...
  localSz = linOct.size();
  par::Mpi_Reduce<DendroIntL>(&localSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    std::cout << GRN " # of Unbalanced Octants: " YLW << totalSz << NRM << std::endl;
  }
  pts.clear();

  //treeNodesTovtk(linOct, rank, "p2o_output");

  // =========== Balancing ============
  MPI_Barrier(MPI_COMM_WORLD);
  if (!rank) {
    std::cout << BLU << "===============================================" << NRM << std::endl;
    std::cout << RED " Starting 2:1 Balance" NRM << std::endl;
    std::cout << BLU << "===============================================" << NRM << std::endl;
  }
#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[1]);
#endif
  startTime = MPI_Wtime();
  ot::balanceOctree(linOct, balOct, dim, maxDepth, incCorner, MPI_COMM_WORLD, NULL, NULL);
  endTime = MPI_Wtime();
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //treeNodesTovtk(balOct,rank,"afBalancing");
  
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

 // treeNodesTovtk(balOct, rank, "bal_output");


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

