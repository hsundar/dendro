
#include "mpi.h"
#include <cstdio>
#include <iostream>
#include <cassert>
#include "parUtils.h"
#include "handleStencils.h"

int createShFnMat(double******& shFnMat) {

#ifdef __USE_MG_INIT_TYPE3__
  createShFnMat_Type3(shFnMat);
#else
#ifdef __USE_MG_INIT_TYPE2__
  createShFnMat_Type2(shFnMat);
#else
  createShFnMat_Type1(shFnMat);
#endif
#endif

  return 1;
}

int createShFnMat_Type3(double******& shFnMat) {
  FILE* infile;
  int rank, res;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  char fname[250];
  sprintf(fname,"ShFnStencils_%d.inp", rank);
  infile = fopen(fname,"r");
  if(!infile) {
    std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
    assert(false);
  }

  typedef double* doublePtr;
  typedef doublePtr* double2Ptr;
  typedef double2Ptr* double3Ptr;
  typedef double3Ptr* double4Ptr;
  typedef double4Ptr* double5Ptr;

  shFnMat = new double5Ptr[8];
  for(unsigned int cNum = 0; cNum < 8; cNum++) {
    shFnMat[cNum] = new double4Ptr[18];
    for(unsigned int eType = 0; eType < 18; eType++) {
      shFnMat[cNum][eType] = new double3Ptr[8];
      for(unsigned int j = 0; j < 8; j++) {
        shFnMat[cNum][eType][j] = new double2Ptr[3];
        for(unsigned int m = 0; m < 3; m++) {
          shFnMat[cNum][eType][j][m] = new doublePtr[3];
          for(unsigned int n = 0; n < 3; n++) {
            shFnMat[cNum][eType][j][m][n] = new double[3];
            for(unsigned int p = 0; p < 3; p++) {
              res = fscanf(infile,"%lf",&(shFnMat[cNum][eType][j][m][n][p]));
            }//end p
          }//end n
        }//end m
      }//end j
    }//end etype
  }//end cNum

  fclose(infile);

  return 1;
}

int createShFnMat_Type2(double******& shFnMat) {
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
    char fname[250];
    sprintf(fname,"ShFnStencils_%d.inp", (rank/THOUSAND));
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
  typedef double4Ptr* double5Ptr;

  shFnMat = new double5Ptr[8];
  for(unsigned int cNum = 0; cNum < 8; cNum++) {
    shFnMat[cNum] = new double4Ptr[18];
    for(unsigned int eType = 0; eType < 18; eType++) {
      shFnMat[cNum][eType] = new double3Ptr[8];
      for(unsigned int j = 0; j < 8; j++) {
        shFnMat[cNum][eType][j] = new double2Ptr[3];
        for(unsigned int m = 0; m < 3; m++) {
          shFnMat[cNum][eType][j][m] = new doublePtr[3];
          for(unsigned int n = 0; n < 3; n++) {
            shFnMat[cNum][eType][j][m][n] = new double[3];
            if((rank % THOUSAND) == 0) {
              for(unsigned int p = 0; p < 3; p++) {
                res = fscanf(infile,"%lf",&(shFnMat[cNum][eType][j][m][n][p]));
              }//end p
            }
          }//end n
        }//end m
      }//end j
    }//end etype
  }//end cNum

  if((rank % THOUSAND) == 0) {
    fclose(infile);
  }

  double * tmpMat = new double[31104];

  if((rank % THOUSAND) == 0) {
    unsigned int ctr = 0;
    for(unsigned int cNum = 0; cNum < 8; cNum++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int j = 0; j < 8; j++) {
          for(unsigned int m = 0; m < 3; m++) {
            for(unsigned int n = 0; n < 3; n++) {
              for(unsigned int p = 0; p < 3; p++) {
                tmpMat[ctr] = shFnMat[cNum][eType][j][m][n][p];
                ctr++;
              }
            }
          }
        }
      }
    }
  }

  par::Mpi_Bcast<double>(tmpMat,31104, 0, newComm);

  if((rank % THOUSAND) != 0) {
    unsigned int ctr = 0;
    for(unsigned int cNum = 0; cNum < 8; cNum++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int j = 0; j < 8; j++) {
          for(unsigned int m = 0; m < 3; m++) {
            for(unsigned int n = 0; n < 3; n++) {
              for(unsigned int p = 0; p < 3; p++) {
                shFnMat[cNum][eType][j][m][n][p] = tmpMat[ctr];
                ctr++;
              }
            }
          }
        }
      }
    }
  }

  delete [] tmpMat;

  return 1;
}//end fn.

int createShFnMat_Type1(double******& shFnMat) {
  FILE* infile;
  int rank, res;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(!rank) {
    char fname[250];
    sprintf(fname,"ShFnStencils.inp");
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
  typedef double4Ptr* double5Ptr;

  shFnMat = new double5Ptr[8];
  for(unsigned int cNum = 0; cNum < 8; cNum++) {
    shFnMat[cNum] = new double4Ptr[18];
    for(unsigned int eType = 0; eType < 18; eType++) {
      shFnMat[cNum][eType] = new double3Ptr[8];
      for(unsigned int j = 0; j < 8; j++) {
        shFnMat[cNum][eType][j] = new double2Ptr[3];
        for(unsigned int m = 0; m < 3; m++) {
          shFnMat[cNum][eType][j][m] = new doublePtr[3];
          for(unsigned int n = 0; n < 3; n++) {
            shFnMat[cNum][eType][j][m][n] = new double[3];
            if(!rank) {
              for(unsigned int p = 0; p < 3; p++) {
                res = fscanf(infile,"%lf",&(shFnMat[cNum][eType][j][m][n][p]));
              }//end p
            }
          }//end n
        }//end m
      }//end j
    }//end etype
  }//end cNum

  if(!rank) {
    fclose(infile);
  }

  double * tmpMat = new double[31104];

  if(!rank) {
    unsigned int ctr = 0;
    for(unsigned int cNum = 0; cNum < 8; cNum++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int j = 0; j < 8; j++) {
          for(unsigned int m = 0; m < 3; m++) {
            for(unsigned int n = 0; n < 3; n++) {
              for(unsigned int p = 0; p < 3; p++) {
                tmpMat[ctr] = shFnMat[cNum][eType][j][m][n][p];
                ctr++;
              }
            }
          }
        }
      }
    }
  }

  par::Mpi_Bcast<double>(tmpMat,31104, 0, MPI_COMM_WORLD);

  if(rank) {
    unsigned int ctr = 0;
    for(unsigned int cNum = 0; cNum < 8; cNum++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int j = 0; j < 8; j++) {
          for(unsigned int m = 0; m < 3; m++) {
            for(unsigned int n = 0; n < 3; n++) {
              for(unsigned int p = 0; p < 3; p++) {
                shFnMat[cNum][eType][j][m][n][p] = tmpMat[ctr];
                ctr++;
              }
            }
          }
        }
      }
    }
  }

  delete [] tmpMat;

  return 1;
}//end fn.

int createGDmatType2(double ****& GDmat) {

#ifdef __USE_MG_INIT_TYPE3__
  createGDmatType2_Type3(GDmat);
#else
#ifdef __USE_MG_INIT_TYPE2__
  createGDmatType2_Type2(GDmat);
#else
  createGDmatType2_Type1(GDmat);
#endif
#endif

  return 1;
}

int createGDmatType2_Type3(double ****& GDmat) {
  FILE* infile;
  int rank, res;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  char fname[250];
  sprintf(fname,"GDType2Stencils_%d.inp", rank);
  infile = fopen(fname,"r");
  if(!infile) {
    std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
    assert(false);
  }

  typedef double* doublePtr;
  typedef doublePtr* double2Ptr;
  typedef double2Ptr* double3Ptr;

  GDmat = new double3Ptr[8];
  for(unsigned int cNum = 0; cNum < 8; cNum++) {
    GDmat[cNum] = new double2Ptr[18];
    for(unsigned int eType = 0; eType < 18; eType++) {
      GDmat[cNum][eType] = new doublePtr[24];
      for(unsigned int i = 0; i < 24; i++) {
        GDmat[cNum][eType][i] = new double[24];
        for(unsigned int j = 0; j < 24; j++) {
          res = fscanf(infile,"%lf",&(GDmat[cNum][eType][i][j]));
        }
      }
    }
  }

  fclose(infile);

  return 1;
}

int createGDmatType2_Type2(double ****& GDmat) {
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
    char fname[250];
    sprintf(fname,"GDType2Stencils_%d.inp", (rank/THOUSAND));
    infile = fopen(fname,"r");
    if(!infile) {
      std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
      assert(false);
    }
  }

  typedef double* doublePtr;
  typedef doublePtr* double2Ptr;
  typedef double2Ptr* double3Ptr;

  GDmat = new double3Ptr[8];
  for(unsigned int cNum = 0; cNum < 8; cNum++) {
    GDmat[cNum] = new double2Ptr[18];
    for(unsigned int eType = 0; eType < 18; eType++) {
      GDmat[cNum][eType] = new doublePtr[24];
      for(unsigned int i = 0; i < 24; i++) {
        GDmat[cNum][eType][i] = new double[24];
        if((rank % THOUSAND) == 0) {
          for(unsigned int j = 0; j < 24; j++) {
            res = fscanf(infile,"%lf",&(GDmat[cNum][eType][i][j]));
          }
        }
      }
    }
  }

  if((rank % THOUSAND) == 0) {
    fclose(infile);
  }

  double * tmpMat = new double[82944];

  if((rank % THOUSAND) == 0) {
    unsigned int ctr = 0;
    for(unsigned int cNum = 0; cNum < 8; cNum++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int i = 0; i < 24; i++) {
          for(unsigned int j = 0; j < 24; j++) {
            tmpMat[ctr] = GDmat[cNum][eType][i][j];
            ctr++;
          }
        }
      }
    }
  }

  par::Mpi_Bcast<double>(tmpMat,82944, 0, newComm);

  if((rank % THOUSAND) != 0) {
    unsigned int ctr = 0;
    for(unsigned int cNum = 0; cNum < 8; cNum++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int i = 0; i < 24; i++) {
          for(unsigned int j = 0; j < 24; j++) {
            GDmat[cNum][eType][i][j] = tmpMat[ctr];
            ctr++;
          }
        }
      }
    }
  }

  delete [] tmpMat;

  return 1;
}//end fn.


/*Type 2 Matrices: Coarse and Fine are the same.*/
int createGDmatType2_Type1(double ****& GDmat) {
  FILE* infile;
  int rank, res;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(!rank) {
    char fname[250];
    sprintf(fname,"GDType2Stencils.inp");
    infile = fopen(fname,"r");
    if(!infile) {
      std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
      assert(false);
    }
  }

  typedef double* doublePtr;
  typedef doublePtr* double2Ptr;
  typedef double2Ptr* double3Ptr;

  GDmat = new double3Ptr[8];
  for(unsigned int cNum = 0; cNum < 8; cNum++) {
    GDmat[cNum] = new double2Ptr[18];
    for(unsigned int eType = 0; eType < 18; eType++) {
      GDmat[cNum][eType] = new doublePtr[24];
      for(unsigned int i = 0; i < 24; i++) {
        GDmat[cNum][eType][i] = new double[24];
        if(!rank) {
          for(unsigned int j = 0; j < 24; j++) {
            res = fscanf(infile,"%lf",&(GDmat[cNum][eType][i][j]));
          }
        }
      }
    }
  }

  if(!rank) {
    fclose(infile);
  }

  double * tmpMat = new double[82944];

  if(!rank) {
    unsigned int ctr = 0;
    for(unsigned int cNum = 0; cNum < 8; cNum++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int i = 0; i < 24; i++) {
          for(unsigned int j = 0; j < 24; j++) {
            tmpMat[ctr] = GDmat[cNum][eType][i][j];
            ctr++;
          }
        }
      }
    }
  }

  par::Mpi_Bcast<double>(tmpMat,82944, 0, MPI_COMM_WORLD);

  if(rank) {
    unsigned int ctr = 0;
    for(unsigned int cNum = 0; cNum < 8; cNum++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int i = 0; i < 24; i++) {
          for(unsigned int j = 0; j < 24; j++) {
            GDmat[cNum][eType][i][j] = tmpMat[ctr];
            ctr++;
          }
        }
      }
    }
  }

  delete [] tmpMat;

  return 1;
}//end fn.

int createMmatType2(double ****& Mmat) {

#ifdef __USE_MG_INIT_TYPE3__
  createMmatType2_Type3(Mmat);
#else
#ifdef __USE_MG_INIT_TYPE2__
  createMmatType2_Type2(Mmat);
#else
  createMmatType2_Type1(Mmat);
#endif
#endif

  return 1;
}

int createMmatType2_Type3(double ****& Mmat) {
  FILE* infile;
  int rank, res;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  char fname[250];
  sprintf(fname,"MassType2Stencils_%d.inp", rank);
  infile = fopen(fname,"r");
  if(!infile) {
    std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
    assert(false);
  }

  typedef double* doublePtr;
  typedef doublePtr* double2Ptr;
  typedef double2Ptr* double3Ptr;

  Mmat = new double3Ptr[8];
  for(unsigned int cNum = 0; cNum < 8; cNum++) {
    Mmat[cNum] = new double2Ptr[18];
    for(unsigned int eType = 0; eType < 18; eType++) {
      Mmat[cNum][eType] = new doublePtr[8];
      for(unsigned int i = 0; i < 8; i++) {
        Mmat[cNum][eType][i] = new double[8];
        for(unsigned int j = 0; j < 8; j++) {
          res = fscanf(infile,"%lf",&(Mmat[cNum][eType][i][j]));
        }
      }
    }
  }

  fclose(infile);

  return 1;
}

int createMmatType2_Type2(double ****& Mmat) {
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
    char fname[250];
    sprintf(fname,"MassType2Stencils_%d.inp", (rank/THOUSAND));
    infile = fopen(fname,"r");
    if(!infile) {
      std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
      assert(false);
    }
  }

  typedef double* doublePtr;
  typedef doublePtr* double2Ptr;
  typedef double2Ptr* double3Ptr;

  Mmat = new double3Ptr[8];
  for(unsigned int cNum = 0; cNum < 8; cNum++) {
    Mmat[cNum] = new double2Ptr[18];
    for(unsigned int eType = 0; eType < 18; eType++) {
      Mmat[cNum][eType] = new doublePtr[8];
      for(unsigned int i = 0; i < 8; i++) {
        Mmat[cNum][eType][i] = new double[8];
        if((rank % THOUSAND) == 0) {
          for(unsigned int j = 0; j < 8; j++) {
            res = fscanf(infile,"%lf",&(Mmat[cNum][eType][i][j]));
          }
        }
      }
    }
  }

  if((rank % THOUSAND) == 0) {
    fclose(infile);
  }

  double * tmpMat = new double[9216];

  if((rank % THOUSAND) == 0) {
    unsigned int ctr = 0;
    for(unsigned int cNum = 0; cNum < 8; cNum++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int i = 0; i < 8; i++) {
          for(unsigned int j = 0; j < 8; j++) {
            tmpMat[ctr] = Mmat[cNum][eType][i][j];
            ctr++;
          }
        }
      }
    }
  }

  par::Mpi_Bcast<double>(tmpMat,9216, 0, newComm);

  if((rank % THOUSAND) != 0) {
    unsigned int ctr = 0;
    for(unsigned int cNum = 0; cNum < 8; cNum++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int i = 0; i < 8; i++) {
          for(unsigned int j = 0; j < 8; j++) {
            Mmat[cNum][eType][i][j] = tmpMat[ctr];
            ctr++;
          }
        }
      }
    }
  }

  delete [] tmpMat;

  return 1;
}//end of function


int createMmatType2_Type1(double ****& Mmat) {
  FILE* infile;
  int rank, res;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(!rank) {
    char fname[250];
    sprintf(fname,"MassType2Stencils.inp");
    infile = fopen(fname,"r");
    if(!infile) {
      std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
      assert(false);
    }
  }

  typedef double* doublePtr;
  typedef doublePtr* double2Ptr;
  typedef double2Ptr* double3Ptr;

  Mmat = new double3Ptr[8];
  for(unsigned int cNum = 0; cNum < 8; cNum++) {
    Mmat[cNum] = new double2Ptr[18];
    for(unsigned int eType = 0; eType < 18; eType++) {
      Mmat[cNum][eType] = new doublePtr[8];
      for(unsigned int i = 0; i < 8; i++) {
        Mmat[cNum][eType][i] = new double[8];
        if(!rank) {
          for(unsigned int j = 0; j < 8; j++) {
            res = fscanf(infile,"%lf",&(Mmat[cNum][eType][i][j]));
          }
        }
      }
    }
  }

  if(!rank) {
    fclose(infile);
  }

  double * tmpMat = new double[9216];

  if(!rank) {
    unsigned int ctr = 0;
    for(unsigned int cNum = 0; cNum < 8; cNum++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int i = 0; i < 8; i++) {
          for(unsigned int j = 0; j < 8; j++) {
            tmpMat[ctr] = Mmat[cNum][eType][i][j];
            ctr++;
          }
        }
      }
    }
  }

  par::Mpi_Bcast<double>(tmpMat,9216, 0, MPI_COMM_WORLD);

  if(rank) {
    unsigned int ctr = 0;
    for(unsigned int cNum = 0; cNum < 8; cNum++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int i = 0; i < 8; i++) {
          for(unsigned int j = 0; j < 8; j++) {
            Mmat[cNum][eType][i][j] = tmpMat[ctr];
            ctr++;
          }
        }
      }
    }
  }

  delete [] tmpMat;

  return 1;
}//end of function

int createLmatType2(double ****& Lmat) {

#ifdef __USE_MG_INIT_TYPE3__
  createLmatType2_Type3(Lmat);
#else
#ifdef __USE_MG_INIT_TYPE2__
  createLmatType2_Type2(Lmat);
#else
  createLmatType2_Type1(Lmat);
#endif
#endif

  return 1;
}

int createLmatType2_Type3(double ****& Lmat) {
  FILE* infile;
  int rank, res;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  char fname[250];
  sprintf(fname,"LapType2Stencils_%d.inp", rank);
  infile = fopen(fname,"r");
  if(!infile) {
    std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
    assert(false);
  }

  typedef double* doublePtr;
  typedef doublePtr* double2Ptr;
  typedef double2Ptr* double3Ptr;

  Lmat = new double3Ptr[8];
  for(unsigned int cNum = 0; cNum < 8; cNum++) {
    Lmat[cNum] = new double2Ptr[18];
    for(unsigned int eType = 0; eType < 18; eType++) {
      Lmat[cNum][eType] = new doublePtr[8];
      for(unsigned int i = 0; i < 8; i++) {
        Lmat[cNum][eType][i] = new double[8];
        for(unsigned int j = 0; j < 8; j++) {
          res = fscanf(infile,"%lf",&(Lmat[cNum][eType][i][j]));
        }
      }
    }
  }

  fclose(infile);

  return 1;
}

int createLmatType2_Type2(double ****& Lmat) {
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
    char fname[250];
    sprintf(fname,"LapType2Stencils_%d.inp", (rank/THOUSAND));
    infile = fopen(fname,"r");
    if(!infile) {
      std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
      assert(false);
    }
  }

  typedef double* doublePtr;
  typedef doublePtr* double2Ptr;
  typedef double2Ptr* double3Ptr;

  Lmat = new double3Ptr[8];
  for(unsigned int cNum = 0; cNum < 8; cNum++) {
    Lmat[cNum] = new double2Ptr[18];
    for(unsigned int eType = 0; eType < 18; eType++) {
      Lmat[cNum][eType] = new doublePtr[8];
      for(unsigned int i = 0; i < 8; i++) {
        Lmat[cNum][eType][i] = new double[8];
        if((rank % THOUSAND) == 0) {
          for(unsigned int j = 0; j < 8; j++) {
            res = fscanf(infile,"%lf",&(Lmat[cNum][eType][i][j]));
          }
        }
      }
    }
  }

  if((rank % THOUSAND) == 0) {
    fclose(infile);
  }

  double * tmpMat = new double[9216];

  if((rank % THOUSAND) == 0) {
    unsigned int ctr = 0;
    for(unsigned int cNum = 0; cNum < 8; cNum++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int i = 0; i < 8; i++) {
          for(unsigned int j = 0; j < 8; j++) {
            tmpMat[ctr] = Lmat[cNum][eType][i][j];
            ctr++;
          }
        }
      }
    }
  }

  par::Mpi_Bcast<double>(tmpMat,9216, 0, newComm);

  if((rank % THOUSAND) != 0) {
    unsigned int ctr = 0;
    for(unsigned int cNum = 0; cNum < 8; cNum++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int i = 0; i < 8; i++) {
          for(unsigned int j = 0; j < 8; j++) {
            Lmat[cNum][eType][i][j] = tmpMat[ctr];
            ctr++;
          }
        }
      }
    }
  }

  delete [] tmpMat;

  return 1;
}//end of function

int createLmatType2_Type1(double ****& Lmat) {
  FILE* infile;
  int rank, res;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(!rank) {
    char fname[250];
    sprintf(fname,"LapType2Stencils.inp");
    infile = fopen(fname,"r");
    if(!infile) {
      std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
      assert(false);
    }
  }

  typedef double* doublePtr;
  typedef doublePtr* double2Ptr;
  typedef double2Ptr* double3Ptr;

  Lmat = new double3Ptr[8];
  for(unsigned int cNum = 0; cNum < 8; cNum++) {
    Lmat[cNum] = new double2Ptr[18];
    for(unsigned int eType = 0; eType < 18; eType++) {
      Lmat[cNum][eType] = new doublePtr[8];
      for(unsigned int i = 0; i < 8; i++) {
        Lmat[cNum][eType][i] = new double[8];
        if(!rank) {
          for(unsigned int j = 0; j < 8; j++) {
            res = fscanf(infile,"%lf",&(Lmat[cNum][eType][i][j]));
          }
        }
      }
    }
  }

  if(!rank) {
    fclose(infile);
  }

  double * tmpMat = new double[9216];

  if(!rank) {
    unsigned int ctr = 0;
    for(unsigned int cNum = 0; cNum < 8; cNum++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int i = 0; i < 8; i++) {
          for(unsigned int j = 0; j < 8; j++) {
            tmpMat[ctr] = Lmat[cNum][eType][i][j];
            ctr++;
          }
        }
      }
    }
  }

  par::Mpi_Bcast<double>(tmpMat,9216, 0, MPI_COMM_WORLD);

  if(rank) {
    unsigned int ctr = 0;
    for(unsigned int cNum = 0; cNum < 8; cNum++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int i = 0; i < 8; i++) {
          for(unsigned int j = 0; j < 8; j++) {
            Lmat[cNum][eType][i][j] = tmpMat[ctr];
            ctr++;
          }
        }
      }
    }
  }

  delete [] tmpMat;

  return 1;
}//end of function

int createRHSType2(double ***& Lmat ) {

#ifdef __USE_MG_INIT_TYPE3__
  createRHSType2_Type3(Lmat);
#else
#ifdef __USE_MG_INIT_TYPE2__
  createRHSType2_Type2(Lmat);
#else
  createRHSType2_Type1(Lmat);
#endif
#endif

  return 1;
}

int createRHSType2_Type3 (double ***& Lmat ) {
  FILE* infile;
  FILE* debugfile;

  int rank, res;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  char fname[250];
  sprintf(fname,"RHSType2Stencils_%d.inp", rank);
  infile = fopen(fname,"r");
  if(!infile) {
    std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
    assert(false);
  }

  typedef double* doublePtr;
  typedef doublePtr* double2Ptr;

  Lmat = new double2Ptr[8];
  for(unsigned int cNum = 0; cNum < 8; cNum++) {
    Lmat[cNum] = new doublePtr[18];
    for(unsigned int eType = 0; eType < 18; eType++) {
      Lmat[cNum][eType] = new double[8];
      for(unsigned int i = 0; i < 8; i++) {
        res = fscanf(infile,"%lf",&(Lmat[cNum][eType][i]));
      }
    }
  }

  fclose(infile);

  return 1;
}

int createRHSType2_Type2(double ***& Lmat ) {
  FILE* infile;
  FILE* debugfile;
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
    char fname[250];
    sprintf(fname,"RHSType2Stencils_%d.inp", (rank/THOUSAND));
    infile = fopen(fname,"r");
    if(!infile) {
      std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
      assert(false);
    }
  }

  typedef double* doublePtr;
  typedef doublePtr* double2Ptr;

  Lmat = new double2Ptr[8];
  for(unsigned int cNum = 0; cNum < 8; cNum++) {
    Lmat[cNum] = new doublePtr[18];
    for(unsigned int eType = 0; eType < 18; eType++) {
      Lmat[cNum][eType] = new double[8];
      if((rank % THOUSAND) == 0) {
        for(unsigned int i = 0; i < 8; i++) {
          res = fscanf(infile,"%lf",&(Lmat[cNum][eType][i]));
        }
      }
    }
  }

  if((rank % THOUSAND) == 0) {
    fclose(infile);
  }

  const int data_size = 8*18*8;
  double * tmpMat = new double[data_size];

  if((rank % THOUSAND) == 0) {
    unsigned int ctr = 0;
    for(unsigned int cNum = 0; cNum < 8; cNum++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int i = 0; i < 8; i++) {
          tmpMat[ctr] = Lmat[cNum][eType][i];
          ctr++;
        }
      }
    }
  }

  par::Mpi_Bcast<double>(tmpMat,data_size, 0, newComm);

  if((rank % THOUSAND) != 0) {
    unsigned int ctr = 0;
    for(unsigned int cNum = 0; cNum < 8; cNum++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int i = 0; i < 8; i++) {
          Lmat[cNum][eType][i] = tmpMat[ctr];
          ctr++;
        }
      }
    }
  }

  delete [] tmpMat;

  return 1;
}

int createRHSType2_Type1 (double ***& Lmat ) {
  FILE* infile;
  FILE* debugfile;

  int rank, res;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(!rank) {
    char fname[250];
    sprintf(fname,"RHSType2Stencils.inp");
    infile = fopen(fname,"r");
    if(!infile) {
      std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
      assert(false);
    }
  }

  typedef double* doublePtr;
  typedef doublePtr* double2Ptr;

  Lmat = new double2Ptr[8];
  for(unsigned int cNum = 0; cNum < 8; cNum++) {
    Lmat[cNum] = new doublePtr[18];
    for(unsigned int eType = 0; eType < 18; eType++) {
      Lmat[cNum][eType] = new double[8];
      if(!rank) {
        for(unsigned int i = 0; i < 8; i++) {
          res = fscanf(infile,"%lf",&(Lmat[cNum][eType][i]));
        }
      }
    }
  }

  if(!rank) {
    fclose(infile);
  }

  const int data_size = 8*18*8;
  double * tmpMat = new double[data_size];

  if(!rank) {
    unsigned int ctr = 0;
    for(unsigned int cNum = 0; cNum < 8; cNum++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int i = 0; i < 8; i++) {
          tmpMat[ctr] = Lmat[cNum][eType][i];
          ctr++;
        }
      }
    }
  }

  par::Mpi_Bcast<double>(tmpMat,data_size, 0, MPI_COMM_WORLD);

  if(rank) {
    unsigned int ctr = 0;
    for(unsigned int cNum = 0; cNum < 8; cNum++) {
      for(unsigned int eType = 0; eType < 18; eType++) {
        for(unsigned int i = 0; i < 8; i++) {
          Lmat[cNum][eType][i] = tmpMat[ctr];
          ctr++;
        }
      }
    }
  }

  delete [] tmpMat;

  return 1;
}

int destroyRHSType2 (double ***& Lmat )
{  
  for(unsigned int cNum = 0; cNum < 8; cNum++) {
    for(unsigned int eType = 0; eType < 18; eType++) {
      delete [] (Lmat[cNum][eType]);
      Lmat[cNum][eType] = NULL;
    }

    delete [] (Lmat[cNum]);
    Lmat[cNum] = NULL;
  }

  delete [] Lmat;  
  Lmat = NULL;

  return 1;
}

int destroyLmatType2(double ****& Lmat ) {
  for(unsigned int cNum = 0; cNum < 8; cNum++) {
    for(unsigned int eType = 0; eType < 18; eType++) {
      for(unsigned int i = 0; i < 8; i++) {
        delete [] (Lmat[cNum][eType][i]);
        Lmat[cNum][eType][i] = NULL;
      }
      delete [] (Lmat[cNum][eType]);
      Lmat[cNum][eType] = NULL;
    }
    delete [] (Lmat[cNum]);
    Lmat[cNum] = NULL;
  }

  delete [] Lmat;
  Lmat = NULL;

  return 1;
}//end of function

int destroyMmatType2(double ****& Mmat ) {
  for(unsigned int cNum = 0; cNum < 8; cNum++) {
    for(unsigned int eType = 0; eType < 18; eType++) {
      for(unsigned int i = 0; i < 8; i++) {
        delete [] (Mmat[cNum][eType][i]);
        Mmat[cNum][eType][i] = NULL;
      }
      delete [] (Mmat[cNum][eType]);
      Mmat[cNum][eType] = NULL;
    }
    delete [] (Mmat[cNum]);
    Mmat[cNum] = NULL;
  }

  delete [] Mmat;
  Mmat = NULL;

  return 1;
}//end of function

int destroyGDmatType2(double ****& GDmat) {
  for(unsigned int cNum = 0; cNum < 8; cNum++) {
    for(unsigned int eType = 0; eType < 18; eType++) {
      for(unsigned int i = 0; i < 24; i++) {
        delete [] (GDmat[cNum][eType][i]);
        GDmat[cNum][eType][i] = NULL;
      }
      delete [] (GDmat[cNum][eType]);
      GDmat[cNum][eType] = NULL;
    }
    delete [] (GDmat[cNum]);
    GDmat[cNum] = NULL;
  }

  delete [] GDmat;
  GDmat = NULL;

  return 1;
}//end fn.

int destroyShFnMat(double******& shFnMat) {
  for(unsigned int cNum = 0; cNum < 8; cNum++) {
    for(unsigned int eType = 0; eType < 18; eType++) {
      for(unsigned int j = 0; j < 8; j++) {
        for(unsigned int m = 0; m < 3; m++) {
          for(unsigned int n = 0; n < 3; n++) {
            delete [] (shFnMat[cNum][eType][j][m][n]);
            shFnMat[cNum][eType][j][m][n] = NULL;
          }
          delete [] (shFnMat[cNum][eType][j][m]);
          shFnMat[cNum][eType][j][m] = NULL;
        }
        delete [] (shFnMat[cNum][eType][j]);
        shFnMat[cNum][eType][j] = NULL;
      }
      delete [] (shFnMat[cNum][eType]);
      shFnMat[cNum][eType] = NULL;
    }
    delete [] (shFnMat[cNum]);
    shFnMat[cNum] = NULL;
  }

  delete [] shFnMat;
  shFnMat = NULL;

  return 1;
}//end fn.



