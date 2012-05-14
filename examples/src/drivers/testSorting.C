#include "mpi.h"
#include "petsc.h"
#include "sys.h"
#include <vector>
#include "TreeNode.h"
#include "parUtils.h"
#include "omg.h"
#include "odaUtils.h"
#include <cstdlib>
#include <cstring>
#include "colors.h"
#include "externVars.h"
#include "dendro.h"

static char help[] = "Evaluates the performance of various parallel sorts\n\n";

//Don't want time to be synchronized. Need to check load imbalance.
#ifdef MPI_WTIME_IS_GLOBAL
#undef MPI_WTIME_IS_GLOBAL
#endif

#ifndef iC
#define iC(fun) {CHKERRQ(fun);}
#endif

double diffclock(clock_t clock1, clock_t clock2) {
  double diffticks = clock2 - clock1;
  double diffms = (diffticks*1000)/CLOCKS_PER_SEC;
  return diffms;
}


int main(int argc, char ** argv ) {	
  int size, rank;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  srand (4*rank);

  std::vector<int> i_in, i_out;
  int n = atoi(argv[1]);
  
  for (int i=0; i<n; i++) {
    i_in.push_back(rand());
  }


  //create a vector to sort ...

  //! sample sort
  MPI_Barrier(MPI_COMM_WORLD);
  clock_t ssort_beg = clock();
  par::sampleSort<int> (i_in, i_out, MPI_COMM_WORLD); 
  MPI_Barrier(MPI_COMM_WORLD);
  clock_t ssort_end = clock();
  
  //! bitonic 
  MPI_Barrier(MPI_COMM_WORLD);
  clock_t bsort_beg = clock();
  par::bitonicSort<int> (i_in, MPI_COMM_WORLD) ;
  MPI_Barrier(MPI_COMM_WORLD);
  clock_t bsort_end = clock();

  if (!rank) {
    std::cout << "Sample sort on " << n << " integers took " < diffclock(ssort_beg, ssort_end) << " ms" << std::endl;
    std::cout << "Bitonic sort on " << n << " integers took " < diffclock(bsort_beg, bsort_end) << " ms" << std::endl;
  }


  MPI_Finalize();
  return 1;
}

