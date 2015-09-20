
/**
  @file handleType1Stencils.C
  @author Rahul Sampath, rahul.sampath@gmail.com 
  */

#include "mpi.h"
#include <cstdio>
#include <cassert>
#include <iostream>
#include "parUtils.h"
#include "handleStencils.h"

/*Type 1 Matrices: Coarse is the parent of the Fine elements.*/

int createMmatType1(double *****& Mmat) {

#ifdef __USE_MG_INIT_TYPE3__
  createMmatType1_Type3(Mmat);
#else
#ifdef __USE_MG_INIT_TYPE2__
  createMmatType1_Type2(Mmat);
#else
  createMmatType1_Type1(Mmat);
#endif
#endif

  return 1;
}

int createMmatType1_Type3(double *****& Mmat) {
  FILE* infile;
  int rank, res;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  char fname[100];
  sprintf(fname,"MassType1Stencils_%d.inp", rank);
  infile = fopen(fname,"r");
  if(!infile) {
    std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
    assert(false);
  }

  typedef double* doublePtr;
  typedef doublePtr* double2Ptr;
  typedef double2Ptr* double3Ptr;
  typedef double3Ptr* double4Ptr;

  Mmat = new double4Ptr[8];
  for(unsigned int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
    Mmat[cNumCoarse] = new double3Ptr[18];
    for(unsigned int eType = 0; eType < 18; eType++) {
      Mmat[cNumCoarse][eType] = new double2Ptr[8];
      for(unsigned int cNumFine = 0; cNumFine < 8; cNumFine++) {
        Mmat[cNumCoarse][eType][cNumFine] = new doublePtr[8];
        for(unsigned int i = 0; i < 8; i++) {
          Mmat[cNumCoarse][eType][cNumFine][i] = new double[8];
          for(unsigned int j = 0; j < 8; j++) {
            res = fscanf(infile,"%lf",&(Mmat[cNumCoarse][eType][cNumFine][i][j]));
          }
        }
      }
    }
  }

  fclose(infile);

  return 1;
}

int createMmatType1_Type2(double *****& Mmat) {
  FILE* infile;
  MPI_Comm comm = MPI_COMM_WORLD;

  int rank, npes, res;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &npes);

  const int THOUSAND = 1000;
  int numGroups = (npes/THOUSAND);
  if( (numGroups*THOUSAND) < npes) {
    numGroups++;
  }

  MPI_Comm newComm;

  bool* isEmptyList = new bool[npes];
  for(int i = 0; i < numGroups; i++) {
    for(int j = 0; (j < (i*THOUSAND)) && (j < npes); j++) {
      isEmptyList[j] = true;
    }
    for(int j = (i*THOUSAND); (j < ((i+1)*THOUSAND)) && (j < npes); j++) {
      isEmptyList[j] = false;
    }
    for(int j = ((i + 1)*THOUSAND); j < npes; j++) {
      isEmptyList[j] = true;
    }
    MPI_Comm tmpComm;
    par::splitComm2way(isEmptyList, &tmpComm, comm);
    if(!(isEmptyList[rank])) {
      newComm = tmpComm;
    }
  }//end for i
  delete [] isEmptyList;

  if((rank % THOUSAND) == 0) {
    char fname[100];
    sprintf(fname,"MassType1Stencils_%d.inp",(rank/THOUSAND));
    infile = fopen(fname,"r");
    if(!infile) {
      std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
      assert(false);
    }
  }

  typedef double* doublePtr;
  typedef doublePtr* double2Ptr;
  typedef double2Ptr* double3Ptr;
  typedef double3Ptr* double4Ptr;

  Mmat = new double4Ptr[8];
  for(unsigned int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
    Mmat[cNumCoarse] = new double3Ptr[18];
    for(unsigned int eType = 0; eType < 18; eType++) {
      Mmat[cNumCoarse][eType] = new double2Ptr[8];
      for(unsigned int cNumFine = 0; cNumFine < 8; cNumFine++) {
        Mmat[cNumCoarse][eType][cNumFine] = new doublePtr[8];
        for(unsigned int i = 0; i < 8; i++) {
          Mmat[cNumCoarse][eType][cNumFine][i] = new double[8];
          if((rank % THOUSAND) == 0) {
            for(unsigned int j = 0; j < 8; j++) {
              res = fscanf(infile,"%lf",&(Mmat[cNumCoarse][eType][cNumFine][i][j]));
            }
          }
        }
      }
    }
  }

  if((rank % THOUSAND) == 0) {
    fclose(infile);
  }

  double * tmpMat = new double[73728];

  if((rank % THOUSAND) == 0) {
    unsigned int ctr = 0;
    for(unsigned int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int cNumFine = 0; cNumFine < 8; cNumFine++) {
          for(unsigned int i = 0; i < 8; i++) {
            for(unsigned int j = 0; j < 8; j++) {
              tmpMat[ctr] = Mmat[cNumCoarse][eType][cNumFine][i][j];
              ctr++;
            }
          }
        }
      }
    }
  }

  par::Mpi_Bcast<double>(tmpMat,73728, 0, newComm);

  if((rank % THOUSAND) != 0) {
    unsigned int ctr = 0;
    for(unsigned int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int cNumFine = 0; cNumFine < 8; cNumFine++) {
          for(unsigned int i = 0; i < 8; i++) {
            for(unsigned int j = 0; j < 8; j++) {
              Mmat[cNumCoarse][eType][cNumFine][i][j] = tmpMat[ctr];
              ctr++;
            }
          }
        }
      }
    }
  }

  delete [] tmpMat;

  return 1;
}//end of function

int createMmatType1_Type1(double *****& Mmat) {
  FILE* infile;
  int rank, res;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if(!rank) {
    char fname[100];
    sprintf(fname,"MassType1Stencils.inp");
    infile = fopen(fname,"r");
    if(!infile) {
      std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
      assert(false);
    }
  }

  typedef double* doublePtr;
  typedef doublePtr* double2Ptr;
  typedef double2Ptr* double3Ptr;
  typedef double3Ptr* double4Ptr;

  Mmat = new double4Ptr[8];
  for(unsigned int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
    Mmat[cNumCoarse] = new double3Ptr[18];
    for(unsigned int eType = 0; eType < 18; eType++) {
      Mmat[cNumCoarse][eType] = new double2Ptr[8];
      for(unsigned int cNumFine = 0; cNumFine < 8; cNumFine++) {
        Mmat[cNumCoarse][eType][cNumFine] = new doublePtr[8];
        for(unsigned int i = 0; i < 8; i++) {
          Mmat[cNumCoarse][eType][cNumFine][i] = new double[8];
          if(!rank) {
            for(unsigned int j = 0; j < 8; j++) {
              res = fscanf(infile,"%lf",&(Mmat[cNumCoarse][eType][cNumFine][i][j]));
            }
          }
        }
      }
    }
  }

  if(!rank) {
    fclose(infile);
  }

  double * tmpMat = new double[73728];

  if(!rank) {
    unsigned int ctr = 0;
    for(unsigned int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int cNumFine = 0; cNumFine < 8; cNumFine++) {
          for(unsigned int i = 0; i < 8; i++) {
            for(unsigned int j = 0; j < 8; j++) {
              tmpMat[ctr] = Mmat[cNumCoarse][eType][cNumFine][i][j];
              ctr++;
            }
          }
        }
      }
    }
  }

  par::Mpi_Bcast<double>(tmpMat,73728, 0, MPI_COMM_WORLD);

  if(rank) {
    unsigned int ctr = 0;
    for(unsigned int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int cNumFine = 0; cNumFine < 8; cNumFine++) {
          for(unsigned int i = 0; i < 8; i++) {
            for(unsigned int j = 0; j < 8; j++) {
              Mmat[cNumCoarse][eType][cNumFine][i][j] = tmpMat[ctr];
              ctr++;
            }
          }
        }
      }
    }
  }

  delete [] tmpMat;

  return 1;
}//end of function

int createLmatType1(double *****& Lmat) {

#ifdef __USE_MG_INIT_TYPE3__
  createLmatType1_Type3(Lmat);
#else
#ifdef __USE_MG_INIT_TYPE2__
  createLmatType1_Type2(Lmat);
#else
  createLmatType1_Type1(Lmat);
#endif
#endif

  return 1;
}

int createLmatType1_Type3(double *****& Lmat) {
  FILE* infile;
  int rank, res;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  char fname[100];
  sprintf(fname,"LapType1Stencils_%d.inp", rank); 	
  infile = fopen(fname,"r");
  if(!infile) {
    std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
    assert(false);
  }

  typedef double* doublePtr;
  typedef doublePtr* double2Ptr;
  typedef double2Ptr* double3Ptr;
  typedef double3Ptr* double4Ptr;

  Lmat = new double4Ptr[8];
  for(unsigned int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
    Lmat[cNumCoarse] = new double3Ptr[18];
    for(unsigned int eType = 0; eType < 18; eType++) {
      Lmat[cNumCoarse][eType] = new double2Ptr[8];
      for(unsigned int cNumFine = 0; cNumFine < 8; cNumFine++) {
        Lmat[cNumCoarse][eType][cNumFine] = new doublePtr[8];
        for(unsigned int i = 0; i < 8; i++) {
          Lmat[cNumCoarse][eType][cNumFine][i] = new double[8];
          for(unsigned int j = 0; j < 8; j++) {
            res = fscanf(infile,"%lf",&(Lmat[cNumCoarse][eType][cNumFine][i][j]));
          }
        }
      }
    }
  }

  fclose(infile);

  return 1;
}

int createLmatType1_Type2(double *****& Lmat) {
  FILE* infile;
  MPI_Comm comm = MPI_COMM_WORLD;

  int rank, npes, res;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &npes);

  const int THOUSAND = 1000;
  int numGroups = (npes/THOUSAND);
  if( (numGroups*THOUSAND) < npes) {
    numGroups++;
  }

  MPI_Comm newComm;

  bool* isEmptyList = new bool[npes];
  for(int i = 0; i < numGroups; i++) {
    for(int j = 0; (j < (i*THOUSAND)) && (j < npes); j++) {
      isEmptyList[j] = true;
    }
    for(int j = (i*THOUSAND); (j < ((i+1)*THOUSAND)) && (j < npes); j++) {
      isEmptyList[j] = false;
    }
    for(int j = ((i + 1)*THOUSAND); j < npes; j++) {
      isEmptyList[j] = true;
    }
    MPI_Comm tmpComm;
    par::splitComm2way(isEmptyList, &tmpComm, comm);
    if(!(isEmptyList[rank])) {
      newComm = tmpComm;
    }
  }//end for i
  delete [] isEmptyList;

  if((rank % THOUSAND) == 0) {
    char fname[100];
    sprintf(fname,"LapType1Stencils_%d.inp",(rank/THOUSAND)); 	
    infile = fopen(fname,"r");
    if(!infile) {
      std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
      assert(false);
    }
  }

  typedef double* doublePtr;
  typedef doublePtr* double2Ptr;
  typedef double2Ptr* double3Ptr;
  typedef double3Ptr* double4Ptr;

  Lmat = new double4Ptr[8];
  for(unsigned int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
    Lmat[cNumCoarse] = new double3Ptr[18];
    for(unsigned int eType = 0; eType < 18; eType++) {
      Lmat[cNumCoarse][eType] = new double2Ptr[8];
      for(unsigned int cNumFine = 0; cNumFine < 8; cNumFine++) {
        Lmat[cNumCoarse][eType][cNumFine] = new doublePtr[8];
        for(unsigned int i = 0; i < 8; i++) {
          Lmat[cNumCoarse][eType][cNumFine][i] = new double[8];
          if((rank % THOUSAND) == 0) {
            for(unsigned int j = 0; j < 8; j++) {
              res = fscanf(infile,"%lf",&(Lmat[cNumCoarse][eType][cNumFine][i][j]));
            }
          }
        }
      }
    }
  }

  if((rank % THOUSAND) == 0) {
    fclose(infile);
  }

  double * tmpMat = new double[73728];

  if((rank % THOUSAND) == 0) {
    unsigned int ctr = 0;
    for(unsigned int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int cNumFine = 0; cNumFine < 8; cNumFine++) {
          for(unsigned int i = 0; i < 8; i++) {
            for(unsigned int j = 0; j < 8; j++) {
              tmpMat[ctr] = Lmat[cNumCoarse][eType][cNumFine][i][j];
              ctr++;
            }
          }
        }
      }
    }
  }

  par::Mpi_Bcast<double>(tmpMat,73728, 0, newComm);

  if((rank % THOUSAND) != 0) {
    unsigned int ctr = 0;
    for(unsigned int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int cNumFine = 0; cNumFine < 8; cNumFine++) {
          for(unsigned int i = 0; i < 8; i++) {
            for(unsigned int j = 0; j < 8; j++) {
              Lmat[cNumCoarse][eType][cNumFine][i][j] = tmpMat[ctr];
              ctr++;
            }
          }
        }
      }
    }
  }

  delete [] tmpMat;

  return 1;
}//end of function


int createLmatType1_Type1(double *****& Lmat) {
  FILE* infile;
  int rank, res;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if(!rank) {
    char fname[100];
    sprintf(fname,"LapType1Stencils.inp"); 	
    infile = fopen(fname,"r");
    if(!infile) {
      std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
      assert(false);
    }
  }

  typedef double* doublePtr;
  typedef doublePtr* double2Ptr;
  typedef double2Ptr* double3Ptr;
  typedef double3Ptr* double4Ptr;

  Lmat = new double4Ptr[8];
  for(unsigned int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
    Lmat[cNumCoarse] = new double3Ptr[18];
    for(unsigned int eType = 0; eType < 18; eType++) {
      Lmat[cNumCoarse][eType] = new double2Ptr[8];
      for(unsigned int cNumFine = 0; cNumFine < 8; cNumFine++) {
        Lmat[cNumCoarse][eType][cNumFine] = new doublePtr[8];
        for(unsigned int i = 0; i < 8; i++) {
          Lmat[cNumCoarse][eType][cNumFine][i] = new double[8];
          if(!rank) {
            for(unsigned int j = 0; j < 8; j++) {
              res = fscanf(infile,"%lf",&(Lmat[cNumCoarse][eType][cNumFine][i][j]));
            }
          }
        }
      }
    }
  }

  if(!rank) {
    fclose(infile);
  }

  double * tmpMat = new double[73728];

  if(!rank) {
    unsigned int ctr = 0;
    for(unsigned int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int cNumFine = 0; cNumFine < 8; cNumFine++) {
          for(unsigned int i = 0; i < 8; i++) {
            for(unsigned int j = 0; j < 8; j++) {
              tmpMat[ctr] = Lmat[cNumCoarse][eType][cNumFine][i][j];
              ctr++;
            }
          }
        }
      }
    }
  }

  par::Mpi_Bcast<double>(tmpMat,73728, 0, MPI_COMM_WORLD);

  if(rank) {
    unsigned int ctr = 0;
    for(unsigned int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int cNumFine = 0; cNumFine < 8; cNumFine++) {
          for(unsigned int i = 0; i < 8; i++) {
            for(unsigned int j = 0; j < 8; j++) {
              Lmat[cNumCoarse][eType][cNumFine][i][j] = tmpMat[ctr];
              ctr++;
            }
          }
        }
      }
    }
  }

  delete [] tmpMat;

  return 1;
}//end of function

int destroyLmatType1(double *****& Lmat ) {
  for(unsigned int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
    for(unsigned int eType = 0; eType < 18; eType++) {
      for(unsigned int cNumFine = 0; cNumFine < 8; cNumFine++) {
        for(unsigned int i = 0; i < 8; i++) {
          delete [] (Lmat[cNumCoarse][eType][cNumFine][i]);
          Lmat[cNumCoarse][eType][cNumFine][i] = NULL;
        }
        delete [] (Lmat[cNumCoarse][eType][cNumFine]);
        Lmat[cNumCoarse][eType][cNumFine] = NULL;
      }
      delete [] (Lmat[cNumCoarse][eType]);
      Lmat[cNumCoarse][eType] = NULL;
    }
    delete [] (Lmat[cNumCoarse]);
    Lmat[cNumCoarse] = NULL;
  }

  delete [] Lmat;
  Lmat = NULL;

  return 1;
}//end of function

int destroyMmatType1(double *****& Mmat ) {
  for(unsigned int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
    for(unsigned int eType = 0; eType < 18; eType++) {
      for(unsigned int cNumFine = 0; cNumFine < 8; cNumFine++) {
        for(unsigned int i = 0; i < 8; i++) {
          delete [] (Mmat[cNumCoarse][eType][cNumFine][i]);
          Mmat[cNumCoarse][eType][cNumFine][i] = NULL;
        }
        delete [] (Mmat[cNumCoarse][eType][cNumFine]);
        Mmat[cNumCoarse][eType][cNumFine] = NULL;
      }
      delete [] (Mmat[cNumCoarse][eType]);
      Mmat[cNumCoarse][eType] = NULL;
    }
    delete [] (Mmat[cNumCoarse]);
    Mmat[cNumCoarse] = NULL;
  }

  delete [] Mmat;
  Mmat = NULL;

  return 1;
}//end of function



