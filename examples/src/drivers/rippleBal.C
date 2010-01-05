
#include "mpi.h"
#include "petsc.h"
#include "sys.h"
#include "parUtils.h"
#include "octUtils.h"
#include "TreeNode.h"
#include <cstdlib>
#include <cstring>
#include "externVars.h"
#include "dendro.h"

//Don't want time to be synchronized. Need to check load imbalance.
#ifdef MPI_WTIME_IS_GLOBAL
#undef MPI_WTIME_IS_GLOBAL
#endif

int main(int argc, char ** argv ) {	
  int size, rank;
  bool incCorner = 1; 
  bool rePart = 1;
  bool checkBailOut = 1; 
  unsigned int writeBOut = 0;
  unsigned int rippleType = 2;
  char Kstr[20];
  char inpFileName[50],n2oOut[50],balOut[50];

  PetscInitialize(&argc,&argv,"options",NULL);
  ot::RegisterEvents();

#ifdef PETSC_USE_LOG
  int stages[1];
  PetscLogStageRegister("Bal",&stages[0]);
#endif

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(argc < 2) {
    std::cerr << "Usage: " << argv[0] << " inpfile writeBOut[0]"<<
      " incCorner[1] checkBailOut[1] rePart[1] rippleType[2] "  << std::endl;
    return -1;
  }
  if(argc > 2) {
    writeBOut = atoi(argv[2]);
  }
  if(argc > 3) { incCorner = (bool)(atoi(argv[3]));}
  if(argc > 4) { checkBailOut = (bool)(atoi(argv[4]));}
  if(argc > 5) { rePart = (bool)(atoi(argv[5]));}
  if(argc > 6) { rippleType = atoi(argv[6]);}

  strcpy(inpFileName, argv[1]);
  strcpy(balOut,inpFileName);
  ot::int2str(rank,Kstr);
  strcat(balOut,Kstr);
  strcat(balOut,"_\0");
  ot::int2str(size,Kstr);
  strcat(balOut,Kstr);
  strcpy(n2oOut,balOut);
  strcat(balOut,"_Bal.ot\0");
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
  ot::completeOctree(nodes,linOct,dim,maxDepth,false, false, false, MPI_COMM_WORLD);
  if(rippleType == 1) {
    if(!rank) {
      std::cout<<"Type 1 ripple."<<std::endl;
    }
    ot::parallelRippleType1(linOct,incCorner,checkBailOut,rePart,dim, maxDepth, MPI_COMM_WORLD);
  }else if(rippleType == 2) {
    if(!rank) {
      std::cout<<"Type 2 ripple."<<std::endl;
    }
    ot::parallelRippleType2(linOct,incCorner,checkBailOut,rePart,dim, maxDepth, MPI_COMM_WORLD);
  }else {
    if(!rank) {
      std::cout<<"Type 3 ripple."<<std::endl;
    }
    ot::parallelRippleType3(linOct,incCorner,checkBailOut,rePart,dim, maxDepth, MPI_COMM_WORLD);
  }
#ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif

  DendroIntL locOutSize = linOct.size();
  DendroIntL globOutSize;

  par::Mpi_Reduce<DendroIntL>(&locOutSize, &globOutSize, 1, MPI_SUM, 0, MPI_COMM_WORLD);

  if(!rank) {
	  std::cout<<"Final Balanced Octree Size: "<<globOutSize<<std::endl;
  }


  if(writeBOut) { 
    ot::writeNodesToFile(balOut,linOct);
  }

  PetscFinalize();

}

