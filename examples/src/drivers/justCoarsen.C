
#include "mpi.h"
#include "petsc.h"
#include "sys.h"
#include "octUtils.h"
#include "TreeNode.h"
#include <cstdlib>
#include <cstring>
#include "externVars.h"

//Don't want time to be synchronized. Need to check load imbalance.
#ifdef MPI_WTIME_IS_GLOBAL
#undef MPI_WTIME_IS_GLOBAL
#endif

int main(int argc, char ** argv ) {	
  int size, rank;
  unsigned int writeCOut = 1;
  char Kstr[20];
  char inpFileName[50],n2oOut[50],linOut[50];

  PetscInitialize(&argc,&argv,"options",NULL);
  ot::RegisterEvents();

#ifdef PETSC_USE_LOG
  int stages[1];
  PetscLogStageRegister(&stages[0],"Coarsen");
#else
  MPI_Init(&argc,&argv);
#endif
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(argc < 2) {
    std::cerr << "Usage: " << argv[0] << " inpfile writeCOut[1]"<< std::endl;
    return -1;
  }
  if(argc > 2) {
    writeCOut = atoi(argv[2]);
  }

  strcpy(inpFileName, argv[1]);
  strcpy(linOut,inpFileName);
  ot::int2str(rank,Kstr);
  strcat(linOut,Kstr);
  strcat(linOut,"_\0");
  ot::int2str(size,Kstr);
  strcat(linOut,Kstr);
  strcpy(n2oOut,linOut);
  strcat(linOut,"_C.ot\0");
  strcat(n2oOut,".ot\0");
  std::vector<ot::TreeNode> nodes;
  std::vector<ot::TreeNode> linOct;

      ot::readNodesFromFile(n2oOut,nodes);
      unsigned int dim = nodes[0].getDim();
      unsigned int maxDepth = nodes[0].getMaxDepth();
    MPI_Barrier(MPI_COMM_WORLD);	
#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[0]);
#endif
  ot::coarsenOctree(nodes, linOct, dim, maxDepth, MPI_COMM_WORLD, false, NULL, NULL);
#ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif
    if(writeCOut) { 
      ot::writeNodesToFile(linOut,linOct);
    }

  PetscFinalize();

}

