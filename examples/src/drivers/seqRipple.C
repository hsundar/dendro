
#include "mpi.h"
#include "petsc.h"
#include "sys.h"
#include "TreeNode.h"
#include <iostream>
#include "externVars.h"

int main(int argc, char ** argv ) {	
  bool incCorner = 1; 
  unsigned int writeBOut = 0;
  char n2oOut[256],balOut[256];

  PetscInitialize(&argc,&argv,"options",NULL);
  ot::RegisterEvents();

  if(argc < 2) {
    std::cerr << "Usage: " << argv[0] << " inpfile writeBOut[0]"<<
      " incCorner[1] "  << std::endl;
    return -1;
  }
  if(argc > 2) {
    writeBOut = atoi(argv[2]);
  }
  if(argc > 3) { incCorner = (bool)(atoi(argv[3]));}

  strcpy(n2oOut, argv[1]);
  strcpy(balOut,n2oOut);
  strcat(balOut,"_Bal.ot\0");
  strcat(n2oOut,".ot\0");
  std::vector<ot::TreeNode> nodes;
  std::vector<ot::TreeNode> linOct;


  ot::readNodesFromFile(n2oOut,nodes);
  unsigned int dim = nodes[0].getDim();
  unsigned int maxDepth = nodes[0].getMaxDepth();

  ot::TreeNode root(dim,maxDepth);

  ot::completeSubtree(root, nodes, linOct, dim, maxDepth, false, false);

  std::cout<<"Initial octree size: "<<linOct.size()<<std::endl;

  ot::ripple(linOct, incCorner);

  std::cout<<"Final octree size: "<<linOct.size()<<std::endl;

  if(writeBOut) { 
    ot::writeNodesToFile(balOut,linOct);
  }

  PetscFinalize();
}

