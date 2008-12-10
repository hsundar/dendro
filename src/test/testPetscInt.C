
#include <mpi.h>
#include "petsc.h"
#include <iostream>

int main(int argc, char** argv) {
  PetscInitialize(&argc, &argv, NULL, NULL);
  std::cout << " sizeof(PetscInt) = " << (sizeof(PetscInt)) << std::endl;
  PetscFinalize();
}

