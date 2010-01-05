
/**
  @file buildDA.C
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

//Don't want time to be synchronized. Need to check load imbalance.
#ifdef MPI_WTIME_IS_GLOBAL
#undef MPI_WTIME_IS_GLOBAL
#endif

int main(int argc, char ** argv ) {	
  int size, rank;
  char bFile[50];
  bool compressLut=false;
  double localTime, totalTime;
  double startTime, endTime;
  DendroIntL locSz, totalSz;
  std::vector<ot::TreeNode> balOct;

  PetscInitialize(&argc,&argv,"options",NULL);
  ot::RegisterEvents();
  ot::DA_Initialize(MPI_COMM_WORLD);

#ifdef PETSC_USE_LOG
  int stages[1];
  PetscLogStageRegister("ODACreate",&stages[0]);
#endif

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(argc < 2) {
    std::cerr << "Usage: " << argv[0] << "inpfile compressLut[0] " << std::endl;
    return -1;
  }
  if(argc > 2) { compressLut = (bool)(atoi(argv[2]));}

  sprintf(bFile,"%s%d_%d.ot",argv[1],rank,size);	

  if(!rank){
    std::cout << " reading  "<<bFile<<std::endl; // Point size
  }
  ot::readNodesFromFile (bFile,balOct);
  if(!rank){
    std::cout << " finished reading  "<<bFile<<std::endl; // Point size
  }
  // compute total inp size and output size
  locSz = balOct.size();
  par::Mpi_Reduce<DendroIntL>(&locSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  par::Mpi_Reduce<double>(&localTime, &totalTime, 1, MPI_MAX, 0, MPI_COMM_WORLD);

  if(!rank) {
    std::cout << "# of Balanced Octants: "<< totalSz << std::endl;       
  }

  //ODA ...
  MPI_Barrier(MPI_COMM_WORLD);	
#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[0]);
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
  locSz = da.getNodeSize();
  localTime = endTime - startTime;
  par::Mpi_Reduce<DendroIntL>(&locSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  par::Mpi_Reduce<double>(&localTime, &totalTime, 1, MPI_MAX, 0, MPI_COMM_WORLD);

  if(!rank) {
    std::cout << "Total # Vertices: "<< totalSz << std::endl;       
    std::cout << "Time to build ODA: "<<totalTime << std::endl;
  }

  if (!rank) {
    std::cout << GRN << "Finalizing ..." << NRM << std::endl;
  }

  ot::DA_Finalize();
  PetscFinalize();
}//end function

