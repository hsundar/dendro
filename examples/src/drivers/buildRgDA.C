
/**
  @file buildRgDA.C
  @author Rahul Sampath, rahul.sampath@gmail.com
  */

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

int main(int argc, char ** argv ) {	
  int size, rank;
  bool compressLut=false;
  unsigned int regLev = 2;
  std::vector<ot::TreeNode> balOct;

  PetscInitialize(&argc,&argv,"options",NULL);
  ot::RegisterEvents();
  ot::DA_Initialize(MPI_COMM_WORLD);

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(argc < 2) {
    std::cerr << "Usage: " << argv[0] << "regLev compressLut[0] " << std::endl;
    return -1;
  }

  regLev = atoi(argv[1]);
  if(argc > 2) { compressLut = (bool)(atoi(argv[2]));}

  if(!rank) {
    std::cout<<"sizeof(DendroIntL) = "<<sizeof(DendroIntL)<<" sizeof(PetscInt) = "<<sizeof(PetscInt)<<std::endl;
  }

  ot::createRegularOctree(balOct, regLev, 3, 30, MPI_COMM_WORLD);

  //ODA ...
  MPI_Barrier(MPI_COMM_WORLD);	
  assert(!(balOct.empty()));
  ot::DA da(balOct, MPI_COMM_WORLD, MPI_COMM_WORLD, compressLut);
  balOct.clear();

  ot::DA_Finalize();
  PetscFinalize();
}//end function

