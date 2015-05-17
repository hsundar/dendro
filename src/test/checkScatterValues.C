
#include "mpi.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include "petscvec.h"
#include "dtypes.h"
#include "omg.h"
#include "parUtils.h"
#include "externVars.h"
#include "dendro.h"

int main(int argc, char ** argv ) {	

  int rank, npes;

  PetscInitialize(&argc,&argv,NULL,NULL);

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&npes);

  //PETSc Vec Version...
  {
    //Ensure that total size of the vectors are the same.
    srand( time(NULL) );
    int inSz = (rank+1)*( rand() % atoi(argv[1]) );

    std::cout<<rank<<" inSz: "<<inSz<<std::endl;

    Vec in;
    VecCreate(MPI_COMM_WORLD,&in);
    VecSetSizes(in,inSz,PETSC_DECIDE);
    std::cout<<rank<<" Created Vec (in) "<<std::endl;
    if(npes > 1) {
      VecSetType(in,VECMPI);
      std::cout<<rank<<" Type: MPI "<<std::endl;
    }else {
      VecSetType(in,VECSEQ);
      std::cout<<rank<<" Type: SEQ "<<std::endl;
    }

    int offRank;

    if(!rank) {
      srand( time(NULL) );
      offRank = ( rand() % npes );
    }

    par::Mpi_Bcast<int>(&offRank,1,0,MPI_COMM_WORLD);

    std::cout<<rank<<" offRank: "<<offRank<<std::endl;

    int *inSizes = new int[npes];
    par::Mpi_Allgather<int>(&inSz, inSizes, 1, MPI_COMM_WORLD);		

    int *inDisps = new int[npes];
    inDisps[0] = 0;
    for(int i = 1; i < npes; i++) {
      inDisps[i] = inDisps[i-1] + inSizes[i-1];
    }

    //CYCLIC SHIFT.
    int outSz = inSizes[(rank + offRank)%npes];

    std::cout<<rank<<" outSz: "<<outSz<<std::endl;

    Vec out;
    VecCreate(MPI_COMM_WORLD, &out);
    VecSetSizes(out,outSz,PETSC_DECIDE);
    std::cout<<rank<<" Created Vec (out) "<<std::endl;
    if(npes > 1) {
      VecSetType(out,VECMPI);
      std::cout<<rank<<" Type: MPI "<<std::endl;
    }else {
      VecSetType(out,VECSEQ);
      std::cout<<rank<<" Type: SEQ "<<std::endl;
    }

    int *outSizes = new int[npes];
    par::Mpi_Allgather<int>(&outSz, outSizes, 1, MPI_COMM_WORLD);	

    int *outDisps = new int[npes];
    outDisps[0] = 0;
    for(int i = 1; i < npes; i++) {
      outDisps[i] = outDisps[i-1] + outSizes[i-1];
    }

    DendroIntL totalInSz;
    DendroIntL totalOutSz;

    DendroIntL inSzCopy = inSz;
    DendroIntL outSzCopy = outSz;

    par::Mpi_Reduce<DendroIntL>(&inSzCopy, &totalInSz, 1, MPI_SUM,0, MPI_COMM_WORLD);	
    par::Mpi_Reduce<DendroIntL>(&outSzCopy, &totalOutSz, 1, MPI_SUM,0, MPI_COMM_WORLD);

    if(!rank) {
      assert(totalInSz == totalOutSz);
      std::cout<<"Global sizes of in and out match."<<std::endl;
    }


    PetscScalar * inArr;
    VecZeroEntries(in);
    VecGetArray(in,&inArr);

    srand( time(NULL) );
    for(int i = 0; i < inSz; i++) {
      inArr[i] = ((double)(rank+1))*(((double)rand())/((double)RAND_MAX));
    }

    PetscScalar *fullInArr = NULL;	
    if(!rank) {
      fullInArr = new PetscScalar[totalInSz];		
    }

    MPI_Gatherv(inArr,inSz,par::Mpi_datatype<PetscScalar>::value(),
        fullInArr,inSizes,inDisps,
        par::Mpi_datatype<PetscScalar>::value(),0,MPI_COMM_WORLD);

    VecRestoreArray(in,&inArr);
    delete [] inSizes;	
    delete [] inDisps;	

    int *sendSz = NULL;
    int *sendOff = NULL;
    int *recvSz = NULL;
    int *recvOff = NULL;

    ot::scatterValues(in, out, inSz, outSz, sendSz, sendOff, recvSz, recvOff, MPI_COMM_WORLD);

    PetscInt tmpSzIn, tmpSzOut;
    VecGetLocalSize(in, &tmpSzIn);
    VecGetLocalSize(out, &tmpSzOut);

    assert(tmpSzIn == inSz);
    assert(tmpSzOut == outSz);

    MPI_Barrier(MPI_COMM_WORLD);
    if(!rank) {
      std::cout<<"The local sizes are preserved."<<std::endl;
    }

    VecDestroy(&in);

    PetscScalar* outArr;
    VecGetArray(out,&outArr);

    PetscScalar* fullOutArr = NULL;
    if(!rank) {
      fullOutArr = new PetscScalar[totalOutSz];
    }

    MPI_Gatherv(outArr, outSz, par::Mpi_datatype<PetscScalar>::value(),
        fullOutArr, outSizes, outDisps, par::Mpi_datatype<PetscScalar>::value(),
        0,MPI_COMM_WORLD);

    delete [] outSizes;
    delete [] outDisps;		
    VecRestoreArray(out,&outArr);
    VecDestroy(&out);

    if(!rank) {
      for(DendroIntL i = 0; i < totalInSz; i++) {
        assert(fullInArr[i] == fullOutArr[i]);
      }
      delete [] fullInArr;	
      delete [] fullOutArr;
      std::cout<<"Passed Vec Version."<<std::endl;
    }
  }//end PETSc version

  MPI_Barrier(MPI_COMM_WORLD);

  //STL version...
  {
    //Ensure that total size of the vectors are the same.
    srand( time(NULL) );
    int inSz = (rank+1)*( rand() % atoi(argv[1]) );

    std::cout<<rank<<" inSz: "<<inSz<<std::endl;

    std::vector<int> in(inSz);
    std::cout<<rank<<" Created std::vec (in) "<<std::endl;

    int offRank;

    if(!rank) {
      srand( time(NULL) );
      offRank = ( rand() % npes );
    }

    par::Mpi_Bcast<int>(&offRank,1,0,MPI_COMM_WORLD);

    std::cout<<rank<<" offRank: "<<offRank<<std::endl;

    int *inSizes = new int[npes];
    par::Mpi_Allgather<int>(&inSz, inSizes, 1, MPI_COMM_WORLD);		

    int *inDisps = new int[npes];
    inDisps[0] = 0;
    for(int i = 1; i < npes; i++) {
      inDisps[i] = inDisps[i-1] + inSizes[i-1];
    }

    //CYCLIC SHIFT.
    int outSz = inSizes[(rank + offRank)%npes];

    std::cout<<rank<<" outSz: "<<outSz<<std::endl;

    std::vector<int> out(outSz);
    std::cout<<rank<<" Created std::vec (out) "<<std::endl;

    int *outSizes = new int[npes];
    par::Mpi_Allgather<int>(&outSz, outSizes, 1, MPI_COMM_WORLD);	

    int *outDisps = new int[npes];
    outDisps[0] = 0;
    for(int i = 1; i < npes; i++) {
      outDisps[i] = outDisps[i-1] + outSizes[i-1];
    }

    DendroIntL totalInSz;
    DendroIntL totalOutSz;

    DendroIntL inSzCopy = inSz;
    DendroIntL outSzCopy = outSz;

    par::Mpi_Reduce<DendroIntL>(&inSzCopy, &totalInSz, 1, MPI_SUM,0, MPI_COMM_WORLD);	
    par::Mpi_Reduce<DendroIntL>(&outSzCopy, &totalOutSz, 1, MPI_SUM,0, MPI_COMM_WORLD);

    if(!rank) {
      assert(totalInSz == totalOutSz);
      std::cout<<"Global sizes of in and out match."<<std::endl;
    }

    srand( time(NULL) );
    for(int i = 0; i < inSz; i++) {
      in[i] = rand()*(rank+1);
    }

    int *fullInArr = NULL;	
    if(!rank) {
      fullInArr = new int[totalInSz];		
    }

    MPI_Gatherv(&(*in.begin()),inSz,MPI_INT,fullInArr,inSizes,inDisps,MPI_INT,0,MPI_COMM_WORLD);

    delete [] inSizes;	
    delete [] inDisps;	

    par::scatterValues<int>(in, out, outSz, MPI_COMM_WORLD);

    int tmpSzIn = in.size();
    int tmpSzOut = out.size();

    assert(tmpSzIn == inSz);
    assert(tmpSzOut == outSz);

    MPI_Barrier(MPI_COMM_WORLD);
    if(!rank) {
      std::cout<<"The local sizes are preserved."<<std::endl;
    }

    in.clear();

    int *fullOutArr = NULL;
    if(!rank) {
      fullOutArr = new int[totalOutSz];
    }

    MPI_Gatherv(&(*out.begin()),outSz,MPI_INT,fullOutArr,outSizes,outDisps,MPI_INT,0,MPI_COMM_WORLD);

    delete [] outSizes;
    delete [] outDisps;		

    out.clear();

    if(!rank) {
      for(int i = 0; i < totalInSz; i++) {
        assert(fullInArr[i] == fullOutArr[i]);
      }
      delete [] fullInArr;	
      delete [] fullOutArr;
      std::cout<<"Passed STL Version."<<std::endl;
    }
  }//end STL version

  MPI_Barrier(MPI_COMM_WORLD);

  PetscFinalize();

}//end main

