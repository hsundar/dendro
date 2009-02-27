
/**
  @file omg.C
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  */

#include "dtypes.h"
#include "petscpc.h"
#include "petscmg.h"
#include "petscmat.h"
#include "private/pcimpl.h"
#include "omg.h"
#include "oda.h"
#include "odaUtils.h" 
#include "parUtils.h"
#include <iostream>
#include "dendro.h"

#ifndef iC
#define iC(fun) {CHKERRQ(fun);}
#endif

#ifdef __DEBUG__
#ifndef __DEBUG_MG__
#define __DEBUG_MG__
#endif
#endif

namespace ot {

  extern double ***** RmatType1Stencil;
  extern double **** RmatType2Stencil;
  extern unsigned short**** VtxMap1; 
  extern unsigned short***** VtxMap2; 
  extern unsigned short***** VtxMap3; 
  extern unsigned short****** VtxMap4; 

  extern void (*getPrivateMatricesForKSP_Shell)(Mat mat,
      Mat *AmatPrivate, Mat *PmatPrivate, MatStructure* pFlag);

  //Public Functions

  int DAMGCreateSuppressedDOFs(DAMG* damg) {
    int       nlevels = damg[0]->nlevels; //number of multigrid levels
    unsigned int dof = damg[0]->dof;
    for(int i = 0; i < nlevels; i++) {
      unsigned int sz = (dof*(damg[i]->da->getLocalBufferSize()));
      if(sz) {
        damg[i]->suppressedDOF = new unsigned char[sz];
      }

      if(damg[i]->da_aux) {
        unsigned int sz2 = (dof*(damg[i]->da_aux->getLocalBufferSize()));
        if(sz2) {
          damg[i]->suppressedDOFaux = new unsigned char[sz2];
        }
      }
    }
    return 1;
  }

  PetscErrorCode DAMG_Initialize(MPI_Comm comm) {
    PROF_MG_INIT_BEGIN 

      ot::DA_Initialize(comm);

#ifdef __USE_MG_INIT_TYPE3__
    ot::DAMG_InitPrivateType3(comm);
#else
#ifdef __USE_MG_INIT_TYPE2__
    ot::DAMG_InitPrivateType2(comm);
#else
    ot::DAMG_InitPrivateType1(comm);
#endif
#endif

    PROF_MG_INIT_END  
  }

  void DAMG_InitPrivateType3(MPI_Comm comm) {

    int rank;
    MPI_Comm_rank(comm, &rank);

    IreadRmatType1Stencil(RmatType1Stencil, rank);
    IreadRmatType2Stencil(RmatType2Stencil, rank);
    IreadVtxMaps(VtxMap1, VtxMap2, VtxMap3, VtxMap4, rank);

  }

  void DAMG_InitPrivateType2(MPI_Comm comm) {

    int rank, npes;
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

    if( (rank % THOUSAND) == 0) {
      IreadRmatType1Stencil(RmatType1Stencil, (rank/THOUSAND));
      IreadRmatType2Stencil(RmatType2Stencil, (rank/THOUSAND));
      IreadVtxMaps(VtxMap1, VtxMap2, VtxMap3, VtxMap4, (rank/THOUSAND));
    } else {
      //Other processors simply allocate the required amount of memory
      typedef double**** double4Ptr;
      typedef double*** double3Ptr;
      typedef double** double2Ptr;
      typedef double* doublePtr;

      RmatType1Stencil = new double4Ptr[8];
      for(int i = 0; i < 8; i++) {
        RmatType1Stencil[i] = new double3Ptr[8];
        for(int j = 0; j < 8; j++) {
          RmatType1Stencil[i][j] = new double2Ptr[18];
          for(int k = 0; k < 18; k++) {
            RmatType1Stencil[i][j][k] = new doublePtr[8];
            for(int l = 0; l < 8; l++) {
              RmatType1Stencil[i][j][k][l] = new double[8];
            }//end for l
          }//end for k
        }//end for j
      }//end for i

      RmatType2Stencil  = new double3Ptr[8];
      for(int j = 0; j < 8; j++) {
        RmatType2Stencil[j] = new double2Ptr[18];
        for(int k = 0; k < 18; k++) {
          RmatType2Stencil[j][k] = new doublePtr[8];
          for(int l = 0; l < 8; l++) {
            RmatType2Stencil[j][k][l] = new double[8];
          }//end for l
        }//end for k
      }//end for j

      typedef unsigned short***** us5Ptr;
      typedef unsigned short**** us4Ptr;
      typedef unsigned short*** us3Ptr;
      typedef unsigned short** us2Ptr;
      typedef unsigned short* usPtr;

      VtxMap1 = new us3Ptr[8];
      VtxMap2 = new us4Ptr[8];
      VtxMap3 = new us4Ptr[7];
      VtxMap4 = new us5Ptr[7];
      for(int i = 0; i < 8; i++) {
        VtxMap1[i] = new us2Ptr[8];
        VtxMap2[i] = new us3Ptr[8];
        for(int j = 0; j < 8; j++) {
          VtxMap1[i][j] = new usPtr[18];
          VtxMap2[i][j] = new us2Ptr[8];
          for(int k = 0; k < 18; k++) {
            VtxMap1[i][j][k] = new unsigned short[8];
          }
          for(int k = 0; k < 8; k++) {
            VtxMap2[i][j][k] = new usPtr[18];
            for(int l = 0; l < 18; l++) {
              VtxMap2[i][j][k][l] = new unsigned short[8];
            }//end for l
          }//end for k
        }//end for j
      }//end for i

      for(int i = 0; i < 7; i++) {
        VtxMap3[i] = new us3Ptr[2];
        VtxMap4[i] = new us4Ptr[2];
        for(int j = 0; j < 2; j++) {
          VtxMap3[i][j] = new us2Ptr[8];
          VtxMap4[i][j] = new us3Ptr[8]; 
          for(int k = 0; k < 8; k++) {
            VtxMap3[i][j][k] = new usPtr[18];
            VtxMap4[i][j][k] = new us2Ptr[8];
            for(int l = 0; l < 18; l++) {
              VtxMap3[i][j][k][l] =new unsigned short[8];
            }
            for(int l = 0; l < 8; l++) {
              VtxMap4[i][j][k][l] = new usPtr[18];
              for(int m = 0; m < 18; m++) {
                VtxMap4[i][j][k][l][m] = new unsigned short[8];
              }
            }//end for l
          }//end for k
        }//end for j
      }//end for i

    }//end if processor reads

    //Processor 0 in each comm  Broadcasts to other processors in the comm
    //doubles ...
    //RmatType1[8][8][18][8][8]: 73728
    //RmatType2[8][18][8][8]: 9216
    double * tmpRmats = new double [82944];

    if((rank % THOUSAND) == 0) {
      unsigned int ctr = 0;
      for(int i = 0; i < 8; i++) {
        for(int j = 0; j < 8; j++) {
          for(int k = 0; k < 18; k++) {
            for(int l = 0; l < 8; l++) {
              for(int m = 0; m < 8; m++) {
                tmpRmats[ctr] = RmatType1Stencil[i][j][k][l][m];
                ctr++;
              }//end for m
            }//end for l
          }//end for k
        }//end for j
      }//end for i
      for(int j = 0; j < 8; j++) {
        for(int k = 0; k < 18; k++) {
          for(int l = 0; l < 8; l++) {
            for(int m = 0; m < 8; m++) {
              tmpRmats[ctr] = RmatType2Stencil[j][k][l][m];
              ctr++;
            }//end for m
          }//end for l
        }//end for k
      }//end for j
    }

    par::Mpi_Bcast<double>(tmpRmats, 82944, 0, newComm);

    if((rank % THOUSAND) != 0) {
      unsigned int ctr = 0;
      for(int i = 0; i < 8; i++) {
        for(int j = 0; j < 8; j++) {
          for(int k = 0; k < 18; k++) {
            for(int l = 0; l < 8; l++) {
              for(int m = 0; m < 8; m++) {
                RmatType1Stencil[i][j][k][l][m] = tmpRmats[ctr];
                ctr++;
              }//end for m
            }//end for l
          }//end for k
        }//end for j
      }//end for i
      for(int j = 0; j < 8; j++) {
        for(int k = 0; k < 18; k++) {
          for(int l = 0; l < 8; l++) {
            for(int m = 0; m < 8; m++) {
              RmatType2Stencil[j][k][l][m] = tmpRmats[ctr];
              ctr++;
            }//end for m
          }//end for l
        }//end for k
      }//end for j
    }

    delete [] tmpRmats;

    //unsigned shorts...
    //map1[8][8][18][8]: 9216
    //map2[8][8][8][18][8]: 73728
    //map3[7][2][8][18][8]: 16128
    //map4[7][2][8][8][18][8]: 129024
    unsigned short * tmpVtxMaps = new unsigned short[228096];

    if((rank % THOUSAND) == 0) {
      unsigned int ctr = 0;
      for(int i = 0; i < 8; i++) {
        for(int j = 0; j < 8; j++) {
          for(int k = 0; k < 18; k++) {
            for(int l = 0; l < 8; l++) {
              tmpVtxMaps[ctr] = VtxMap1[i][j][k][l]; 
              ctr++;
            }//end for l
          }//end for k
        }//end for j
      }//end for i
      for(int i = 0; i < 8; i++) {
        for(int j = 0; j < 8; j++) {
          for(int k = 0; k < 8; k++) {
            for(int l = 0; l < 18; l++) {
              for(int m = 0; m < 8; m++) {
                tmpVtxMaps[ctr] = VtxMap2[i][j][k][l][m]; 
                ctr++;
              }//end for m
            }//end for l
          }//end for k
        }//end for j
      }//end for i
      for(int i = 0; i < 7; i++) {
        for(int j = 0; j < 2; j++) {
          for(int k = 0; k < 8; k++) {
            for(int l = 0; l < 18; l++) {
              for(int m = 0; m < 8; m++) {
                tmpVtxMaps[ctr] = VtxMap3[i][j][k][l][m]; 
                ctr++;
              }//end for m
            }//end for l
          }//end for k
        }//end for j
      }//end for i
      for(int i = 0; i < 7; i++) {
        for(int j = 0; j < 2; j++) {
          for(int k = 0; k < 8; k++) {
            for(int l = 0; l < 8; l++) {
              for(int m = 0; m < 18; m++) {
                for(int n = 0; n < 8; n++) {
                  tmpVtxMaps[ctr] = VtxMap4[i][j][k][l][m][n]; 
                  ctr++;
                }//end for n
              }//end for m
            }//end for l
          }//end for k
        }//end for j
      }//end for i
    }

    par::Mpi_Bcast<unsigned short>(tmpVtxMaps, 228096, 0, newComm);

    if((rank % THOUSAND) != 0) {
      unsigned int ctr = 0;
      for(int i = 0; i < 8; i++) {
        for(int j = 0; j < 8; j++) {
          for(int k = 0; k < 18; k++) {
            for(int l = 0; l < 8; l++) {
              VtxMap1[i][j][k][l] = tmpVtxMaps[ctr]; 
              ctr++;
            }//end for l
          }//end for k
        }//end for j
      }//end for i
      for(int i = 0; i < 8; i++) {
        for(int j = 0; j < 8; j++) {
          for(int k = 0; k < 8; k++) {
            for(int l = 0; l < 18; l++) {
              for(int m = 0; m < 8; m++) {
                VtxMap2[i][j][k][l][m] = tmpVtxMaps[ctr]; 
                ctr++;
              }//end for m
            }//end for l
          }//end for k
        }//end for j
      }//end for i
      for(int i = 0; i < 7; i++) {
        for(int j = 0; j < 2; j++) {
          for(int k = 0; k < 8; k++) {
            for(int l = 0; l < 18; l++) {
              for(int m = 0; m < 8; m++) {
                VtxMap3[i][j][k][l][m] = tmpVtxMaps[ctr]; 
                ctr++;
              }//end for m
            }//end for l
          }//end for k
        }//end for j
      }//end for i
      for(int i = 0; i < 7; i++) {
        for(int j = 0; j < 2; j++) {
          for(int k = 0; k < 8; k++) {
            for(int l = 0; l < 8; l++) {
              for(int m = 0; m < 18; m++) {
                for(int n = 0; n < 8; n++) {
                  VtxMap4[i][j][k][l][m][n] = tmpVtxMaps[ctr]; 
                  ctr++;
                }//end for n
              }//end for m
            }//end for l
          }//end for k
        }//end for j
      }//end for i
    }

    delete [] tmpVtxMaps;

  }

  void DAMG_InitPrivateType1(MPI_Comm comm) {

    int rank;
    MPI_Comm_rank(comm, &rank);

    if(!rank) {
      //Processor 0 reads the stencils
      readRmatType1Stencil(RmatType1Stencil);
      readRmatType2Stencil(RmatType2Stencil);
      readVtxMaps(VtxMap1, VtxMap2, VtxMap3, VtxMap4);
    } else {
      //Other processors simply allocate the required amount of memory
      typedef double**** double4Ptr;
      typedef double*** double3Ptr;
      typedef double** double2Ptr;
      typedef double* doublePtr;

      RmatType1Stencil = new double4Ptr[8];
      for(int i = 0; i < 8; i++) {
        RmatType1Stencil[i] = new double3Ptr[8];
        for(int j = 0; j < 8; j++) {
          RmatType1Stencil[i][j] = new double2Ptr[18];
          for(int k = 0; k < 18; k++) {
            RmatType1Stencil[i][j][k] = new doublePtr[8];
            for(int l = 0; l < 8; l++) {
              RmatType1Stencil[i][j][k][l] = new double[8];
            }
          }
        }
      }

      RmatType2Stencil  = new double3Ptr[8];
      for(int j = 0; j < 8; j++) {
        RmatType2Stencil[j] = new double2Ptr[18];
        for(int k = 0; k < 18; k++) {
          RmatType2Stencil[j][k] = new doublePtr[8];
          for(int l = 0; l < 8; l++) {
            RmatType2Stencil[j][k][l] = new double[8];
          }
        }
      }

      typedef unsigned short***** us5Ptr;
      typedef unsigned short**** us4Ptr;
      typedef unsigned short*** us3Ptr;
      typedef unsigned short** us2Ptr;
      typedef unsigned short* usPtr;

      VtxMap1 = new us3Ptr[8];
      VtxMap2 = new us4Ptr[8];
      VtxMap3 = new us4Ptr[7];
      VtxMap4 = new us5Ptr[7];
      for(int i = 0; i < 8; i++) {
        VtxMap1[i] = new us2Ptr[8];
        VtxMap2[i] = new us3Ptr[8];
        for(int j = 0; j < 8; j++) {
          VtxMap1[i][j] = new usPtr[18];
          VtxMap2[i][j] = new us2Ptr[8];
          for(int k = 0; k < 18; k++) {
            VtxMap1[i][j][k] = new unsigned short[8];
          }
          for(int k = 0; k < 8; k++) {
            VtxMap2[i][j][k] = new usPtr[18];
            for(int l = 0; l < 18; l++) {
              VtxMap2[i][j][k][l] = new unsigned short[8];
            }
          }
        }
      }//end for i

      for(int i = 0; i < 7; i++) {
        VtxMap3[i] = new us3Ptr[2];
        VtxMap4[i] = new us4Ptr[2];
        for(int j = 0; j < 2; j++) {
          VtxMap3[i][j] = new us2Ptr[8];
          VtxMap4[i][j] = new us3Ptr[8]; 
          for(int k = 0; k < 8; k++) {
            VtxMap3[i][j][k] = new usPtr[18];
            VtxMap4[i][j][k] = new us2Ptr[8];
            for(int l = 0; l < 18; l++) {
              VtxMap3[i][j][k][l] =new unsigned short[8];
            }
            for(int l = 0; l < 8; l++) {
              VtxMap4[i][j][k][l] = new usPtr[18];
              for(int m = 0; m < 18; m++) {
                VtxMap4[i][j][k][l][m] = new unsigned short[8];
              }
            }
          }
        }
      }//end for i

    }//end if p0

    //Processor 0 Broadcasts to other processors
    //doubles ...
    //RmatType1[8][8][18][8][8]: 73728
    //RmatType2[8][18][8][8]: 9216
    double * tmpRmats = new double [82944];

    if(!rank) {
      unsigned int ctr = 0;
      for(int i = 0; i < 8; i++) {
        for(int j = 0; j < 8; j++) {
          for(int k = 0; k < 18; k++) {
            for(int l = 0; l < 8; l++) {
              for(int m = 0; m < 8; m++) {
                tmpRmats[ctr] = RmatType1Stencil[i][j][k][l][m];
                ctr++;
              }//end for m
            }//end for l
          }//end for k
        }//end for j
      }//end for i
      for(int j = 0; j < 8; j++) {
        for(int k = 0; k < 18; k++) {
          for(int l = 0; l < 8; l++) {
            for(int m = 0; m < 8; m++) {
              tmpRmats[ctr] = RmatType2Stencil[j][k][l][m];
              ctr++;
            }//end for m
          }//end for l
        }//end for k
      }//end for j
    }

    par::Mpi_Bcast<double>(tmpRmats, 82944, 0, comm);

    if(rank) {
      unsigned int ctr = 0;
      for(int i = 0; i < 8; i++) {
        for(int j = 0; j < 8; j++) {
          for(int k = 0; k < 18; k++) {
            for(int l = 0; l < 8; l++) {
              for(int m = 0; m < 8; m++) {
                RmatType1Stencil[i][j][k][l][m] = tmpRmats[ctr];
                ctr++;
              }//end for m
            }//end for l
          }//end for k
        }//end for j
      }//end for i
      for(int j = 0; j < 8; j++) {
        for(int k = 0; k < 18; k++) {
          for(int l = 0; l < 8; l++) {
            for(int m = 0; m < 8; m++) {
              RmatType2Stencil[j][k][l][m] = tmpRmats[ctr];
              ctr++;
            }//end for m
          }//end for l
        }//end for k
      }//end for j
    }

    delete [] tmpRmats;

    //unsigned shorts...
    //map1[8][8][18][8]: 9216
    //map2[8][8][8][18][8]: 73728
    //map3[7][2][8][18][8]: 16128
    //map4[7][2][8][8][18][8]: 129024
    unsigned short * tmpVtxMaps = new unsigned short[228096];

    if(!rank) {
      unsigned int ctr = 0;
      for(int i = 0; i < 8; i++) {
        for(int j = 0; j < 8; j++) {
          for(int k = 0; k < 18; k++) {
            for(int l = 0; l < 8; l++) {
              tmpVtxMaps[ctr] = VtxMap1[i][j][k][l]; 
              ctr++;
            }//end for l
          }//end for k
        }//end for j
      }//end for i
      for(int i = 0; i < 8; i++) {
        for(int j = 0; j < 8; j++) {
          for(int k = 0; k < 8; k++) {
            for(int l = 0; l < 18; l++) {
              for(int m = 0; m < 8; m++) {
                tmpVtxMaps[ctr] = VtxMap2[i][j][k][l][m]; 
                ctr++;
              }//end for m
            }//end for l
          }//end for k
        }//end for j
      }//end for i
      for(int i = 0; i < 7; i++) {
        for(int j = 0; j < 2; j++) {
          for(int k = 0; k < 8; k++) {
            for(int l = 0; l < 18; l++) {
              for(int m = 0; m < 8; m++) {
                tmpVtxMaps[ctr] = VtxMap3[i][j][k][l][m]; 
                ctr++;
              }//end for m
            }//end for l
          }//end for k
        }//end for j
      }//end for i
      for(int i = 0; i < 7; i++) {
        for(int j = 0; j < 2; j++) {
          for(int k = 0; k < 8; k++) {
            for(int l = 0; l < 8; l++) {
              for(int m = 0; m < 18; m++) {
                for(int n = 0; n < 8; n++) {
                  tmpVtxMaps[ctr] = VtxMap4[i][j][k][l][m][n]; 
                  ctr++;
                }//end for n
              }//end for m
            }//end for l
          }//end for k
        }//end for j
      }//end for i
    }

    par::Mpi_Bcast<unsigned short>(tmpVtxMaps, 228096, 0, comm);

    if(rank) {
      unsigned int ctr = 0;
      for(int i = 0; i < 8; i++) {
        for(int j = 0; j < 8; j++) {
          for(int k = 0; k < 18; k++) {
            for(int l = 0; l < 8; l++) {
              VtxMap1[i][j][k][l] = tmpVtxMaps[ctr]; 
              ctr++;
            }//end for l
          }//end for k
        }//end for j
      }//end for i
      for(int i = 0; i < 8; i++) {
        for(int j = 0; j < 8; j++) {
          for(int k = 0; k < 8; k++) {
            for(int l = 0; l < 18; l++) {
              for(int m = 0; m < 8; m++) {
                VtxMap2[i][j][k][l][m] = tmpVtxMaps[ctr]; 
                ctr++;
              }//end for m
            }//end for l
          }//end for k
        }//end for j
      }//end for i
      for(int i = 0; i < 7; i++) {
        for(int j = 0; j < 2; j++) {
          for(int k = 0; k < 8; k++) {
            for(int l = 0; l < 18; l++) {
              for(int m = 0; m < 8; m++) {
                VtxMap3[i][j][k][l][m] = tmpVtxMaps[ctr]; 
                ctr++;
              }//end for m
            }//end for l
          }//end for k
        }//end for j
      }//end for i
      for(int i = 0; i < 7; i++) {
        for(int j = 0; j < 2; j++) {
          for(int k = 0; k < 8; k++) {
            for(int l = 0; l < 8; l++) {
              for(int m = 0; m < 18; m++) {
                for(int n = 0; n < 8; n++) {
                  VtxMap4[i][j][k][l][m][n] = tmpVtxMaps[ctr]; 
                  ctr++;
                }//end for n
              }//end for m
            }//end for l
          }//end for k
        }//end for j
      }//end for i
    }

    delete [] tmpVtxMaps;

  }

  PetscErrorCode DAMG_Finalize() {
    PROF_MG_FINAL_BEGIN

      ot::DA_Finalize();

    destroyRmatType1Stencil(RmatType1Stencil);
    destroyRmatType2Stencil(RmatType2Stencil);
    destroyVtxMaps(VtxMap1, VtxMap2, VtxMap3, VtxMap4);

    PROF_MG_FINAL_END  
  }

  PetscErrorCode DAMGDestroy(DAMG* damg)
  {
    PetscErrorCode ierr;
    int       i,nlevels = damg[0]->nlevels;

    PetscFunctionBegin;

    int rank;
    MPI_Comm_rank(damg[0]->comm,&rank);

    if (!damg) SETERRQ(PETSC_ERR_ARG_NULL,"Passing null as DAMG");

    for (i = 1; i < nlevels; i++) {
      if (damg[i]->R) {
        ierr = MatDestroy(damg[i]->R);
        CHKERRQ(ierr);
      }
    }

    for (i = 0; i < nlevels; i++) {
      if (damg[i]->x)       {
        ierr = VecDestroy(damg[i]->x);
        CHKERRQ(ierr);
      }

      if (damg[i]->b)       {
        ierr = VecDestroy(damg[i]->b);
        CHKERRQ(ierr);
      }

      if (damg[i]->r)       {
        ierr = VecDestroy(damg[i]->r);
        CHKERRQ(ierr);
      }

      if(damg[i]->suppressedDOF) {
        delete [] (damg[i]->suppressedDOF);
        damg[i]->suppressedDOF = NULL;
      }

      if(damg[i]->suppressedDOFaux) {
        delete [] (damg[i]->suppressedDOFaux);
        damg[i]->suppressedDOFaux = NULL;
      }

      if (damg[i]->B && (damg[i]->B != damg[i]->J)) {
        ierr = MatDestroy(damg[i]->B);CHKERRQ(ierr);
      }

      if (damg[i]->J)         {
        ierr = MatDestroy(damg[i]->J);CHKERRQ(ierr);
      }

      if (damg[i]->ksp) {
        ierr = KSPDestroy(damg[i]->ksp);CHKERRQ(ierr);
      }

      if (damg[i]->da)      {
        delete damg[i]->da; 
        damg[i]->da = NULL;
      }

      if (damg[i]->da_aux)      {
        delete damg[i]->da_aux; 
        damg[i]->da_aux = NULL;
      }

      delete damg[i];
      damg[i] = NULL;
    }//end for all levels

    delete [] damg;
    damg = NULL;  
    PetscFunctionReturn(0);
  }

  //RS: This is the function used to create the J matrix. 
  //This must be called before calling DAMGSetKSP or DAMGSetSNES if you want to set
  //J and B to be different. If this is not called, by default J and B will be equal.
  PetscErrorCode DAMGCreateJMatrix(DAMG damg, PetscErrorCode (*crjac)(DAMG,Mat*)) {
    PetscErrorCode ierr;
    PetscFunctionBegin;   	
    if(crjac){
      ierr = (*crjac)(damg,&(damg->J)); CHKERRQ(ierr);
    }else{
      SETERRQ(PETSC_ERR_ARG_NULL,"ot::DA can not create a Matrix for you! Pass a function handle.");  
    }
    PetscFunctionReturn(0);
  }

  //crjac is the function to create the B matrix for each level
  //compJac is an optional function to change the entries of the matrix. This
  //could be used if you don't want to re-build the matrix but only want to
  //change its entries. It is not useful for a single stand-alone solve.
  PetscErrorCode DAMGSetKSP(DAMG *damg, PetscErrorCode (*crjac)(DAMG,Mat*), 
      PetscErrorCode (*compJac)(DAMG,Mat,Mat), PetscErrorCode (*rhs)(DAMG,Vec) )	
  {
    PetscErrorCode ierr;
    int       i,nlevels = damg[0]->nlevels;

    //Galerkin coarsening is NOT an option.
    PROF_MG_SET_KSP_BEGIN

      int rank;
    MPI_Comm_rank(damg[0]->comm, &rank);

    if (!damg) SETERRQ(PETSC_ERR_ARG_NULL,"Passing null as DAMG");  
    if (!crjac) SETERRQ(PETSC_ERR_ARG_NULL,"ot::DA can not create a Matrix for you! Pass a function handle.");

    if (!damg[0]->ksp) {
      /* create solvers for each level if they don't already exist*/
      for (i = 0; i < nlevels; i++) {
        if (!damg[i]->B) {		
          ierr = (*crjac)(damg[i],&(damg[i]->B));CHKERRQ(ierr);
        }
        if (!damg[i]->J) {
#ifdef __DEBUG_MG__
          if(!rank) {	
            std::cout<<"J and B are the same for lev: "<<i<<std::endl;
            fflush(stdout);
          }
#endif
          damg[i]->J = damg[i]->B;
        }else {
#ifdef __DEBUG_MG__
          if(!rank) {	
            std::cout<<"J and B are different for lev: "<<i<<std::endl;
            fflush(stdout);
          }
#endif
          assert(damg[i]->J != damg[i]->B);
        }

        ierr = KSPCreate(damg[i]->comm,&(damg[i]->ksp));CHKERRQ(ierr);

        //No prefix for the finest level
        if(i < (nlevels-1)) {
          char prefix[256];
          sprintf(prefix,"damg_levels_%d_",(int)i);
          KSPSetOptionsPrefix(damg[i]->ksp,prefix);
        }

        ierr = DAMGSetUpLevel(damg,damg[i]->ksp,i+1); CHKERRQ(ierr);

        ierr = KSPSetFromOptions(damg[i]->ksp);CHKERRQ(ierr);

        if( (damg[0]->da->getNpesActive()) < (damg[0]->da->getNpesAll()) ) {
          //If all processors are NOT active on the coarsest grid
          //reset the solver on the coarsest grid to be KSP_Shell 
          PC pc;
          const char* clearOptionPrefix;
          char optionName[256];
          if(i == 0) {
            KSPGetOptionsPrefix(damg[0]->ksp, &clearOptionPrefix);

            sprintf(optionName, "-%sksp_type",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            ierr = KSPSetType(damg[0]->ksp, KSPRICHARDSON); CHKERRQ(ierr);

            sprintf(optionName, "-%sksp_knoll",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            ierr = KSPSetInitialGuessKnoll(damg[0]->ksp, PETSC_FALSE);
            CHKERRQ(ierr);

            sprintf(optionName, "-%sksp_richardson_scale",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            ierr = KSPRichardsonSetScale(damg[0]->ksp, 1.0);
            CHKERRQ(ierr);

            sprintf(optionName, "-%sksp_right_pc",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            sprintf(optionName, "-%sksp_symmetric_pc",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            ierr = KSPSetPreconditionerSide(damg[0]->ksp, PC_LEFT);
            CHKERRQ(ierr);

            sprintf(optionName, "-%sksp_norm_type",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            ierr = KSPSetNormType(damg[0]->ksp, KSP_NORM_NO);
            CHKERRQ(ierr);

            sprintf(optionName, "-%sksp_rtol",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            sprintf(optionName, "-%sksp_atol",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            sprintf(optionName, "-%sksp_divtol",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            sprintf(optionName, "-%sksp_max_it",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            ierr = KSPSetTolerances(damg[0]->ksp, PETSC_DEFAULT, PETSC_DEFAULT,
                PETSC_DEFAULT, 1); 
            CHKERRQ(ierr);

            ierr = KSPSetConvergenceTest(damg[0]->ksp, KSPSkipConverged,
                PETSC_NULL);
            CHKERRQ(ierr);

            ierr = KSPSetInitialGuessNonzero(damg[0]->ksp, PETSC_FALSE);
            CHKERRQ(ierr);

            ierr  = KSPGetPC(damg[0]->ksp, &pc); CHKERRQ(ierr);

            PCGetOptionsPrefix(pc, &clearOptionPrefix);
            sprintf(optionName, "-%spc_type",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            ierr  = PCSetType(pc, PCSHELL); CHKERRQ(ierr);
            ierr = PCShellSetName(pc, "PC_KSP_Shell"); CHKERRQ(ierr);

            PC_KSP_Shell* pcShellContext = new PC_KSP_Shell;
            pcShellContext->sol_private = NULL;
            pcShellContext->rhs_private = NULL;
            pcShellContext->ksp_private = NULL;
            pcShellContext->pc = pc;
            pcShellContext->iAmActive = damg[0]->da->iAmActive();
            pcShellContext->commActive = damg[0]->da->getCommActive();
            ierr = PCShellSetContext(pc, pcShellContext); CHKERRQ(ierr);
            ierr = PCShellSetSetUp(pc, PC_KSP_Shell_SetUp); CHKERRQ(ierr);
            ierr = PCShellSetApply(pc, PC_KSP_Shell_Apply); CHKERRQ(ierr);
            ierr = PCShellSetDestroy(pc, PC_KSP_Shell_Destroy); CHKERRQ(ierr);
          } else {
            ierr  = KSPGetPC(damg[i]->ksp, &pc); CHKERRQ(ierr);
          }

          PetscTruth ismg;
          PetscTypeCompare((PetscObject)pc, PCMG, &ismg);

          if(ismg) {
            PetscTruth useRTLMG;
            ierr = PetscOptionsHasName(PETSC_NULL,"-damg_useRTLMG",&useRTLMG);
            CHKERRQ(ierr);
            KSP lksp;
            if(useRTLMG) {
              while(ismg) {
                ierr = PCMGGetSmoother(pc,0,&lksp);CHKERRQ(ierr);
                ierr  = KSPGetPC(lksp, &pc); CHKERRQ(ierr);
                PetscTypeCompare((PetscObject)pc, PCMG, &ismg);
              }
            } else {
              ierr = PCMGGetSmoother(pc,0,&lksp);CHKERRQ(ierr);
            }

            KSPGetOptionsPrefix(lksp, &clearOptionPrefix);

            sprintf(optionName, "-%sksp_type",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            ierr = KSPSetType(lksp, KSPRICHARDSON); CHKERRQ(ierr);              

            sprintf(optionName, "-%sksp_knoll",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            ierr = KSPSetInitialGuessKnoll(lksp, PETSC_FALSE);
            CHKERRQ(ierr);

            sprintf(optionName, "-%sksp_richardson_scale",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            ierr = KSPRichardsonSetScale(lksp, 1.0);
            CHKERRQ(ierr);

            sprintf(optionName, "-%sksp_right_pc",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            sprintf(optionName, "-%sksp_symmetric_pc",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            ierr = KSPSetPreconditionerSide(lksp, PC_LEFT);
            CHKERRQ(ierr);

            sprintf(optionName, "-%sksp_norm_type",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            ierr = KSPSetNormType(lksp, KSP_NORM_NO);
            CHKERRQ(ierr);

            sprintf(optionName, "-%sksp_rtol",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            sprintf(optionName, "-%sksp_atol",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            sprintf(optionName, "-%sksp_divtol",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            sprintf(optionName, "-%sksp_max_it",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            ierr = KSPSetTolerances(lksp, PETSC_DEFAULT, PETSC_DEFAULT,
                PETSC_DEFAULT, 1); 
            CHKERRQ(ierr);

            ierr = KSPSetConvergenceTest(lksp, KSPSkipConverged, PETSC_NULL);
            CHKERRQ(ierr);

            ierr = KSPSetInitialGuessNonzero(lksp, PETSC_FALSE);
            CHKERRQ(ierr);

            ierr  = KSPGetPC(lksp, &pc); CHKERRQ(ierr);

            PCGetOptionsPrefix(pc, &clearOptionPrefix);
            sprintf(optionName, "-%spc_type",clearOptionPrefix);
            ierr = PetscOptionsClearValue(optionName); CHKERRQ(ierr);
            ierr  = PCSetType(pc, PCSHELL); CHKERRQ(ierr);
            ierr = PCShellSetName(pc, "PC_KSP_Shell"); CHKERRQ(ierr);

            PC_KSP_Shell* pcShellContext = new PC_KSP_Shell;
            pcShellContext->sol_private = NULL;
            pcShellContext->rhs_private = NULL;
            pcShellContext->ksp_private = NULL;
            pcShellContext->pc = pc;
            pcShellContext->iAmActive = damg[0]->da->iAmActive();
            pcShellContext->commActive = damg[0]->da->getCommActive();
            ierr = PCShellSetContext(pc, pcShellContext); CHKERRQ(ierr);
            ierr = PCShellSetSetUp(pc, PC_KSP_Shell_SetUp); CHKERRQ(ierr);
            ierr = PCShellSetApply(pc, PC_KSP_Shell_Apply); CHKERRQ(ierr);
            ierr = PCShellSetDestroy(pc, PC_KSP_Shell_Destroy); CHKERRQ(ierr);
          }//end if PC==MG
        }//end if all procs active on coarsest level

        damg[i]->solve = DAMGSolveKSP;
        damg[i]->rhs   = rhs;
      }//end for i
    }//end if ksp

#ifdef __DEBUG_MG__
    MPI_Barrier(damg[0]->comm);
    if(!rank) {
      std::cout<<"Finished setting up ksp for all levels."<<std::endl;
      fflush(stdout);
    }
    MPI_Barrier(damg[0]->comm);
#endif

    for (i = 0; i < nlevels; i++) {
      if(compJac) {
        ierr = (*compJac)(damg[i],damg[i]->J,damg[i]->B); CHKERRQ(ierr);
      }
      damg[i]->matricesset = PETSC_TRUE;
    }

#ifdef __DEBUG_MG__
    MPI_Barrier(damg[0]->comm);
    if(!rank) {
      std::cout<<"Finished building matrices for all levels."<<std::endl;
      fflush(stdout);
    }
    MPI_Barrier(damg[0]->comm);
#endif

    PROF_MG_SET_KSP_END
  }

  PetscErrorCode DAMGSetNullSpace(DAMG* damg, PetscTruth has_cnst, PetscInt n,
      PetscErrorCode (*func)(DAMG,Vec[]))
  {
    PetscErrorCode ierr;
    int       i, j, nlevels = damg[0]->nlevels;
    Vec            *nulls = 0;
    MatNullSpace   nullsp;
    KSP            iksp;
    PC             pc,ipc;
    PetscTruth     ismg,isred;

    PetscFunctionBegin;
    if (!damg) SETERRQ(PETSC_ERR_ARG_NULL,"Passing null as DAMG");
    if (!damg[0]->ksp) SETERRQ(PETSC_ERR_ORDER,"Must call AFTER DAMGSetKSP() or DAMGSetSNES()");
    if ((n && !func) || (!n && func)) SETERRQ(PETSC_ERR_ARG_INCOMP,"Both n and func() must be set together");
    if (n < 0) SETERRQ1(PETSC_ERR_ARG_OUTOFRANGE,"Cannot have negative number of vectors in null space n = %D",n)

      for (i = 0; i < nlevels; i++) {
        if (n) {
          ierr = VecDuplicateVecs(damg[i]->b,n,&nulls);CHKERRQ(ierr);
          ierr = (*func)(damg[i],nulls);CHKERRQ(ierr);
        }
        ierr = MatNullSpaceCreate(damg[i]->comm,has_cnst,n,nulls,&nullsp);CHKERRQ(ierr);
        ierr = KSPSetNullSpace(damg[i]->ksp,nullsp);CHKERRQ(ierr);
        for (j = i; j < nlevels; j++) {
          ierr = KSPGetPC(damg[j]->ksp,&pc);CHKERRQ(ierr);
          ierr = PetscTypeCompare((PetscObject)pc,PCMG,&ismg);CHKERRQ(ierr);
          if (ismg) {
            ierr = PCMGGetSmoother(pc,i,&iksp);CHKERRQ(ierr);
            ierr = KSPSetNullSpace(iksp, nullsp);CHKERRQ(ierr);
          }
        }
        ierr = MatNullSpaceDestroy(nullsp);CHKERRQ(ierr);
        if (n) {
          ierr = PetscFree(nulls);CHKERRQ(ierr);
        }
      }//end for i

    /* make all the coarse grid solvers have LU shift since they are singular */
    for (i = 0; i < nlevels; i++) {
      ierr = KSPGetPC(damg[i]->ksp,&pc);CHKERRQ(ierr);
      ierr = PetscTypeCompare((PetscObject)pc,PCMG,&ismg);CHKERRQ(ierr);
      if (ismg) {
        ierr = PCMGGetSmoother(pc,0,&iksp);CHKERRQ(ierr);
        ierr = KSPGetPC(iksp,&ipc);CHKERRQ(ierr);
        ierr = PetscTypeCompare((PetscObject)ipc,PCREDUNDANT,&isred);CHKERRQ(ierr);
        if (isred) {
          ierr = PCRedundantGetPC(ipc,&ipc);CHKERRQ(ierr);
        }
        ierr = PCFactorSetShiftPd(ipc,PETSC_TRUE);CHKERRQ(ierr); 
      }
    }//end for i

    PetscFunctionReturn(0);
  }

  PetscErrorCode DAMGSetInitialGuess(DAMG* damg, PetscErrorCode (*guess)(DAMG, Vec)) {
    int i, nlevels = damg[0]->nlevels;
    for(i = 0; i < nlevels; i++) {
      damg[i]->initialguess = guess;
    }
    return(0);
  }

  PetscErrorCode DAMGInitialGuessCurrent(DAMG damg, Vec vec) {
    return (0);
  }

  PetscErrorCode DAMGSolve(DAMG* damg)
  {
    PetscErrorCode ierr;
    int       i,nlevels = damg[0]->nlevels;
    PetscTruth     gridseq,vecmonitor;

    PetscFunctionBegin;
    ierr = PetscOptionsHasName(0,"-damg_grid_sequence",&gridseq);CHKERRQ(ierr);
    ierr = PetscOptionsHasName(0,"-damg_vecmonitor",&vecmonitor);CHKERRQ(ierr);
    if (gridseq) {    
      if(damg[0]->initialguess) {
        (*(damg[0]->initialguess))(damg[0],damg[0]->x);
        KSPSetInitialGuessNonzero(damg[0]->ksp, PETSC_TRUE);
      }
      for (i=0; i<nlevels-1; i++) {
        ierr = (*damg[i]->solve)(damg,i);CHKERRQ(ierr);
        if (vecmonitor) {
          ierr = VecView(damg[i]->x,PETSC_VIEWER_DRAW_(damg[i]->comm));CHKERRQ(ierr);
        }
        ierr = MatInterpolate(damg[i+1]->R,damg[i]->x,damg[i+1]->x);CHKERRQ(ierr);      
        KSPSetInitialGuessNonzero(damg[i+1]->ksp, PETSC_TRUE);
      }    
    }else {
      if(damg[nlevels-1]->initialguess) {
        (*(damg[nlevels-1]->initialguess))(damg[nlevels-1], damg[nlevels-1]->x);
        KSPSetInitialGuessNonzero(damg[nlevels-1]->ksp, PETSC_TRUE);
      }
    }
    ierr = (*DAMGGetFine(damg)->solve)(damg, nlevels-1);CHKERRQ(ierr);
    if (vecmonitor) {
      ierr = VecView(damg[nlevels-1]->x,PETSC_VIEWER_DRAW_(damg[nlevels-1]->comm));CHKERRQ(ierr);
    }  
    PetscFunctionReturn(0);
  }

  PetscErrorCode DAMGSolveKSP(DAMG *damg, int level)
  {
    PetscErrorCode ierr;

    PetscFunctionBegin;

    if (damg[level]->rhs) {
      ierr = (*damg[level]->rhs)(damg[level],damg[level]->b);CHKERRQ(ierr); 
    }

    if (damg[level]->matricesset) {
      ierr = KSPSetOperators(damg[level]->ksp, damg[level]->J,
          damg[level]->B, SAME_NONZERO_PATTERN);CHKERRQ(ierr);
      damg[level]->matricesset = PETSC_FALSE;
    }

    ierr = KSPSolve(damg[level]->ksp, damg[level]->b, damg[level]->x);CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

  PetscErrorCode DAMGSetUpLevel(DAMG* damg, KSP ksp, int nlevels)
  {
    PetscErrorCode ierr;
    PC             pc;
    PetscTruth     ismg,monitor,ismf,isshell,ismffd;
    PetscTruth     useRTLMG;
    KSP            lksp; /* solver internal to the multigrid preconditioner */
    MPI_Comm       *comms,comm;
    PetscViewer    ascii;

    PetscFunctionBegin;
    if (!damg) SETERRQ(PETSC_ERR_ARG_NULL,"Passing null as DAMG");

    ierr = PetscOptionsHasName(PETSC_NULL,"-damg_ksp_monitor",&monitor);CHKERRQ(ierr);
    ierr = PetscOptionsHasName(PETSC_NULL,"-damg_useRTLMG",&useRTLMG);CHKERRQ(ierr);

    if (monitor) {
      ierr = PetscObjectGetComm((PetscObject)ksp,&comm);CHKERRQ(ierr);
      ierr = PetscViewerASCIIOpen(comm,"stdout",&ascii);CHKERRQ(ierr);
      ierr = PetscViewerASCIISetTab(ascii,1+(damg[0]->nlevels)-nlevels);CHKERRQ(ierr);
      ierr = KSPMonitorSet(ksp, KSPMonitorDefault, ascii,
          (PetscErrorCode(*)(void*))PetscViewerDestroy);CHKERRQ(ierr);
    }

    /* use fgmres on outer iteration by default */
    ierr  = KSPSetType(ksp, KSPFGMRES); CHKERRQ(ierr);
    ierr  = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
    ierr  = PCSetType(pc,PCMG); CHKERRQ(ierr);

    const char *mainKSPprefix;
    KSPGetOptionsPrefix(ksp, &mainKSPprefix);

    if(useRTLMG) {
      char prefixName[256];
      sprintf(prefixName, "rtlmg_finest_");
      PCSetOptionsPrefix(pc, mainKSPprefix);
      PCAppendOptionsPrefix(pc, prefixName);
      //finest to coarsest
      for(int i = (nlevels-1); i >= 1; i--) {
        ierr  = PetscMalloc(2*sizeof(MPI_Comm),&comms);CHKERRQ(ierr);
        for (int j = 0; j < 2; j++) {
          comms[j] = damg[j]->comm;
        }
        PCMGSetLevels(pc,2,comms);
        PetscFree(comms);
        PCMGSetType(pc,PC_MG_FULL);
        PetscTypeCompare((PetscObject)pc,PCMG,&ismg);
        if (ismg) {
          /*Finer Level*/
          PCMGGetSmoother(pc,1,&lksp);
          KSPSetOperators(lksp,damg[i]->B,damg[i]->B,DIFFERENT_NONZERO_PATTERN);
          PCMGSetR(pc,1,damg[i]->r);
          PCMGSetResidual(pc,1,PCMGDefaultResidual,damg[i]->B);
          /*Coarser Level*/
          PCMGGetSmoother(pc,0,&lksp);
          KSPSetOperators(lksp,damg[i-1]->J,damg[i-1]->J,DIFFERENT_NONZERO_PATTERN);
          PCMGSetX(pc,0,damg[i-1]->x);
          PCMGSetRhs(pc,0,damg[i-1]->b);
          /* Set interpolation/restriction between levels */
          PCMGSetInterpolation(pc,1,damg[i]->R);
          PCMGSetRestriction(pc,1,damg[i]->R);
          if(i > 1) {
            /*Get PC for the next level*/
            KSPSetType(lksp,KSPFGMRES);
            KSPGetPC(lksp,&pc);
            PCSetType(pc,PCMG);
            sprintf(prefixName, "rtlmg_levels_%d_",((int)(i-1)));
            PCSetOptionsPrefix(pc, mainKSPprefix);
            PCAppendOptionsPrefix(pc, prefixName);
          }
        }//end if mg
      }//end for i
    } else {
      //Standard Scheme same as original DMMG 
      ierr  = PetscMalloc(nlevels*sizeof(MPI_Comm),&comms);CHKERRQ(ierr);
      for (int i = 0; i < nlevels; i++) {
        comms[i] = damg[i]->comm;
      }
      ierr  = PCMGSetLevels(pc,nlevels,comms);CHKERRQ(ierr);
      ierr =  PCMGSetType(pc,PC_MG_FULL);CHKERRQ(ierr);

      ierr  = PetscFree(comms);CHKERRQ(ierr); 

      ierr = PetscTypeCompare((PetscObject)pc,PCMG,&ismg);CHKERRQ(ierr);
      if (ismg) {
        /* set solvers for each level */
        for (int i = 0; i < nlevels; i++) {
          ierr = PCMGGetSmoother(pc, i, &lksp);CHKERRQ(ierr);

          ierr = KSPSetOperators(lksp, damg[i]->J, damg[i]->B,
              DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);

          if (i < nlevels-1) { /* don't set for finest level, they are set in PCApply_MG()*/
            ierr = PCMGSetX(pc,i,damg[i]->x);CHKERRQ(ierr); 
            ierr = PCMGSetRhs(pc,i,damg[i]->b);CHKERRQ(ierr); 
          }
          if (i > 0) {
            ierr = PCMGSetR(pc,i,damg[i]->r);CHKERRQ(ierr); 
            ierr = PCMGSetResidual(pc,i,PCMGDefaultResidual,damg[i]->J);CHKERRQ(ierr);
          }
          if (monitor) {
            ierr = PetscObjectGetComm((PetscObject)lksp,&comm);CHKERRQ(ierr);
            ierr = PetscViewerASCIIOpen(comm,"stdout",&ascii);CHKERRQ(ierr);
            ierr = PetscViewerASCIISetTab(ascii,1+damg[0]->nlevels-i);CHKERRQ(ierr);
            ierr = KSPMonitorSet(lksp,KSPMonitorDefault,ascii,
                (PetscErrorCode(*)(void*))PetscViewerDestroy);CHKERRQ(ierr);
          }
          /* If using a matrix free multiply and did not provide an explicit matrix to build
             the preconditioner then must use no preconditioner 
             */
          ierr = PetscTypeCompare((PetscObject)damg[i]->B,MATSHELL,&isshell);CHKERRQ(ierr);
          ierr = PetscTypeCompare((PetscObject)damg[i]->B,MATDAAD,&ismf);CHKERRQ(ierr);
          ierr = PetscTypeCompare((PetscObject)damg[i]->B,MATMFFD,&ismffd);CHKERRQ(ierr);
          if (isshell || ismf || ismffd) {
            PC  lpc;
            ierr = KSPGetPC(lksp,&lpc);CHKERRQ(ierr);
            //This is only the default value. It can be modified later
            ierr = PCSetType(lpc,PCNONE);CHKERRQ(ierr);
          }
        }//end for i

        /* Set interpolation/restriction between levels */
        for (int i=1; i<nlevels; i++) {
          ierr = PCMGSetInterpolation(pc,i,damg[i]->R);CHKERRQ(ierr); 
          ierr = PCMGSetRestriction(pc,i,damg[i]->R);CHKERRQ(ierr); 
        }//end for i
      }//end if mg
    }//end if TLMG

    PetscFunctionReturn(0);
  }//end fn.

  //level = 0 is the coarsest, level = (nlevels-1) is the finest.
  //nlevels for each level is the number of levels finer than this level.
  //New implementation. Written on April 24, 2008
  PetscErrorCode DAMGCreateAndSetDA(MPI_Comm comm, int & nlevels, 
      void* user, DAMG** damg, std::vector<ot::TreeNode>& finestOctree,
      unsigned int dof, double loadFac,
      bool compressLut, bool incCorner) {

    PROF_MG_SET_DA_BEGIN

      PROF_SET_DA_STAGE1_BEGIN

      PetscErrorCode ierr;
    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

#ifndef __SILENT_MODE__
    if(!rank) {
      std::cout<<"MG-Load Fac: "<<loadFac<<std::endl;
      fflush(stdout);
    }
#endif

    par::partitionW<ot::TreeNode>(finestOctree, NULL, comm);

    assert(nlevels > 0);
    if(finestOctree.empty()) {
      std::cout<<"Processor "<<rank<<
        " called DAMGCreateAndSetDA with an empty finest octree."<<std::endl;
      fflush(stdout);
      assert(false);
    }

    unsigned int dim = finestOctree[0].getDim();
    unsigned int maxDepth = finestOctree[0].getMaxDepth();

    PROF_SET_DA_STAGE1_END
#ifdef __PROF_WITH_BARRIER__
      MPI_Barrier(comm);
#endif
    PROF_SET_DA_STAGE2_BEGIN

      int idxOfCoarsestLev = -1; 

    std::vector<ot::TreeNode>* coarserOctrees = NULL;

    bool* activeStatesInCoarseBal = NULL;
    MPI_Comm* activeCommsInCoarseBal = NULL;
    int* activeNpesInCoarseBal = NULL;

    if(nlevels > 1) {
      coarserOctrees = new std::vector<ot::TreeNode> [nlevels-1];
      activeStatesInCoarseBal = new bool[nlevels - 1];
      activeCommsInCoarseBal = new MPI_Comm[nlevels - 1];
      activeNpesInCoarseBal = new int[nlevels - 1];

      //Default value is false so if a processor becomes inactive at some
      //level and exits the loop it will still have the correct value for the
      //remaining levels 
      for(int i = 0; i < (nlevels-1); i++) {
        activeStatesInCoarseBal[i] = false;
      }//end for i

      MPI_Comm tmpComm1 = comm;
      MPI_Comm tmpComm2;
      bool iAmActiveForCoarsening = true;
      bool repeatLoop = true;
      while( (idxOfCoarsestLev < (nlevels-2)) && (repeatLoop) ) {
        std::vector<ot::TreeNode> tmpOctree;      
        //First coarsen
        if (idxOfCoarsestLev == -1) {
          //We can skip Partition for the first call to coarsen since
          //finestOctree is explicitly partitioned above. 
          ot::coarsenOctree(finestOctree, tmpOctree, dim, maxDepth, 
              tmpComm1, true, &tmpComm2, &iAmActiveForCoarsening);
        } else {
          //We cannot skip Partition here since there is no explicit call to
          //partition at the end of balancing and after the intial block
          //partition, balancing could have introduced a load imbalance across
          //the processors
          ot::coarsenOctree(coarserOctrees[idxOfCoarsestLev], tmpOctree,
              dim, maxDepth, tmpComm1, false, &tmpComm2, &iAmActiveForCoarsening);
        }//end if finest level

        idxOfCoarsestLev++;

        //Then balance
        if(iAmActiveForCoarsening) {
          ot::balanceOctree(tmpOctree, coarserOctrees[idxOfCoarsestLev],
              dim, maxDepth, incCorner, tmpComm2, &tmpComm1, &iAmActiveForCoarsening);
        }
        tmpOctree.clear();

        //All active processors for this level will have the correct comm set
        //and that's all we care. We do not care about the comms for inactive
        //processors
        activeStatesInCoarseBal[idxOfCoarsestLev] = iAmActiveForCoarsening;
        activeCommsInCoarseBal[idxOfCoarsestLev] = tmpComm1;

        if(iAmActiveForCoarsening) {
          MPI_Comm_size(tmpComm1, (activeNpesInCoarseBal + idxOfCoarsestLev));
          if( (activeNpesInCoarseBal[idxOfCoarsestLev] == 1) &&
              (coarserOctrees[idxOfCoarsestLev].size() < 500) ) {
            repeatLoop = false;
          }
        } else {
          assert(coarserOctrees[idxOfCoarsestLev].empty());
          break;
        }

      }//end while
    }//end if initial nlevels > 1

    PROF_SET_DA_STAGE2_END
#ifdef __PROF_WITH_BARRIER__
      MPI_Barrier(comm);
#endif
    PROF_SET_DA_STAGE3_BEGIN

      //Only processor 0 is guaranteed to have the correct idxOfCoarsestLev
      par::Mpi_Bcast<int>(&idxOfCoarsestLev, 1, 0, comm);

    //Reset nlevels
    nlevels = (idxOfCoarsestLev + 2);

    //Need to synchronize activeNpesInCoarseBal across all processors.
    //Only processor 0 is guaranteed to be active at all levels. Hence, only it
    //will have the correct values at all levels. Inactive processors will have
    //junk values.  Similarly, only processor 0 will have the correct
    //globalSizes for all levels since only it is active on all levels.
    //Inactive processors will have junk values.
    if(nlevels > 1) {
      par::Mpi_Bcast<int>(activeNpesInCoarseBal, (nlevels - 1), 0, comm);
    }

#ifdef __DEBUG_MG__
    for(int lev = 0; lev < (nlevels - 1); lev++) {
      char fname[256];
      sprintf(fname, "coarseOctBeforeBdy_%d_%d_%d.ot", lev, rank, npes);
      ot::writeNodesToFile(fname, coarserOctrees[lev]);
    }
#endif

#ifdef __USE_PVT_DA_IN_MG__
    //All processors are active on the finest level
    std::vector<ot::TreeNode> positiveBoundaries;
    ot::addBoundaryNodesType1(finestOctree, positiveBoundaries, dim, maxDepth);

    //update maxdepth
    maxDepth = maxDepth + 1;

    //Most processors will not add any positive boundaries. So, there is no need
    //to unnecessarily distribute positive boundaries on all procs just for
    //sorting. So only the few processors touching the positive boundary will
    //participate in the parallel sort
    MPI_Comm bdyComm;
    par::splitComm2way((positiveBoundaries.empty()), &bdyComm, comm);

    if(!(positiveBoundaries.empty())) {
      //Call Sample Sort  
      std::vector<ot::TreeNode > tmpVecTN;
      par::sampleSort<ot::TreeNode>(positiveBoundaries, tmpVecTN, bdyComm);
      positiveBoundaries = tmpVecTN;
      tmpVecTN.clear();
    }

    par::concatenate<ot::TreeNode>(finestOctree, positiveBoundaries, comm);
    positiveBoundaries.clear();

    //Now for all the coarser levels
    for (int lev = 0; lev < (nlevels - 1); lev++) {
      if(activeStatesInCoarseBal[lev]) {        
        ot::addBoundaryNodesType1(coarserOctrees[lev], positiveBoundaries, dim, maxDepth-1);

        //Most processors will not add any positive boundaries. So, there is no need
        //to unnecessarily distribute positive boundaries on all procs just for
        //sorting. So only the few processors touching the positive boundary will
        //participate in the parallel sort
        par::splitComm2way((positiveBoundaries.empty()), &bdyComm, activeCommsInCoarseBal[lev]);

        if(!(positiveBoundaries.empty())) {
          //Call Sample Sort  
          std::vector<ot::TreeNode > tmpVecTN;
          par::sampleSort<ot::TreeNode>(positiveBoundaries, tmpVecTN, bdyComm);
          positiveBoundaries = tmpVecTN;
          tmpVecTN.clear();
        }

        par::concatenate<ot::TreeNode>(coarserOctrees[lev], positiveBoundaries, activeCommsInCoarseBal[lev]);
        positiveBoundaries.clear();
      }
    }//end for lev

#endif

#ifdef __DEBUG_MG__
    for(int lev = 0; lev < (nlevels - 1); lev++) {
      char fname[256];
      sprintf(fname, "coarseOctAfterBdy_%d_%d_%d.ot", lev, rank, npes);
      ot::writeNodesToFile(fname, coarserOctrees[lev]);
    }
#endif


    //All processors should know the global size at each level
    DendroIntL* localOctreeSizeForThisLevel = new DendroIntL[nlevels];
    DendroIntL* globalOctreeSizeForThisLevel = new DendroIntL[nlevels];
    localOctreeSizeForThisLevel[0] = finestOctree.size();
    for(int lev = 0; lev < (nlevels - 1); lev++) {
      localOctreeSizeForThisLevel[lev + 1] = coarserOctrees[lev].size();
    }//end for lev
    par::Mpi_Allreduce<DendroIntL>(localOctreeSizeForThisLevel,
        globalOctreeSizeForThisLevel, nlevels, MPI_SUM, comm);
    delete [] localOctreeSizeForThisLevel;

#ifdef __DEBUG_MG__
    MPI_Barrier(comm);
    for(int i = 0; i < (nlevels - 1); i++) {
      if(rank < activeNpesInCoarseBal[i]) {
        assert(activeStatesInCoarseBal[i]);
      } else {
        assert(!(activeStatesInCoarseBal[i]));
      }
    }//end for i
    MPI_Barrier(comm);
#endif

    //0 is the finest and (nlevels-1) is the coarsest
    int *maxProcsForThisLevel = new int [nlevels];

    for(int i = 0; i < nlevels; i++) {
      const DendroIntL THOUSAND = 1000;
      if(globalOctreeSizeForThisLevel[i] < (THOUSAND*npes)) {
        int maxProcsToUse = (globalOctreeSizeForThisLevel[i]/THOUSAND);
        if(maxProcsToUse == 0) {
          maxProcsToUse = 1;
        }
        maxProcsForThisLevel[i] = maxProcsToUse;
      } else {
        //use all the processors
        maxProcsForThisLevel[i] = npes;
      }
    }//end for i

    //maxProcsForThisLevel was computed above using the final balanced octree
    //size for that level. However, the octree itself may be distributed on
    //(a) fewer or (b) more processors. The former case occurs when too many
    //new octants were generated during the balancing step. The latter case
    //occurs when too many extra boundaries were removed. Here we measure the
    //discrepancy between the two and decide which one to pick. On one hand if
    //we change the number of active processors, then we must call scatter to 
    //redistribute the octants. On the other hand if the grain size is much
    //larger than what it should be we will unnecessarily increase the cost
    //for subsequent operations. We will do this only for the coarsest level,
    //this will automatically influence the finer levels in the subsequent
    //step. 
    if(nlevels > 1) {
      if(activeNpesInCoarseBal[nlevels - 2] < maxProcsForThisLevel[nlevels - 1]) {
        DendroIntL currentGrainSize = 
          (globalOctreeSizeForThisLevel[nlevels - 1]/activeNpesInCoarseBal[nlevels - 2]);

        DendroIntL expectedGrainSize = 
          (globalOctreeSizeForThisLevel[nlevels - 1]/maxProcsForThisLevel[nlevels - 1]);

        if(static_cast<double>(currentGrainSize) < 
            (1.5*static_cast<double>(expectedGrainSize))) {
          maxProcsForThisLevel[nlevels - 1] = activeNpesInCoarseBal[nlevels - 2];
        }
      }
    }

    //avoid repartitioning at every level.
    //Loop from coarsest to finest
    //If the fine level does not fit on all processors, then
    //the coarser levels will not fit on all p either.
    //The procs for the coarse level will be >= 1/8 that
    //of the finer level and <= the finer level procs. This is because the
    //number of elements on the coarse level can be anything from almost same
    //(just 8 octants coarsened) to a 8 times reduction in size.  
    //So the grain size for the fine octree computed on the coarse grid procs
    //will be >= the true grain size and <= 8*(the true grain size) 
    //If P_coarse = P_fine, then the cost of the matvec goes up, but there is
    //no additional meshing and no scatters during transfers
    //If P_coarse < P_fine, then there is an extra fine grid mesh using P_fine
    //procs and extra scatters, but the cost of the matvec goes down (due to
    //smaller grain size)

    for(int i = (nlevels-2); i >= 0; i--) {
      if(maxProcsForThisLevel[i] > maxProcsForThisLevel[i+1]) {
        DendroIntL avgGrainSize = (globalOctreeSizeForThisLevel[i]/maxProcsForThisLevel[i]);

        DendroIntL avgGrainSizeUsing1LevCoarserNpes = 
          (globalOctreeSizeForThisLevel[i]/maxProcsForThisLevel[i+1]);

        if(static_cast<double>(avgGrainSizeUsing1LevCoarserNpes) <
            (1.5*static_cast<double>(avgGrainSize))) {
          maxProcsForThisLevel[i] = maxProcsForThisLevel[i+1];
        }
      } else {
        assert(maxProcsForThisLevel[i] == maxProcsForThisLevel[i+1]);
      }
    }//end for i

#ifndef __SILENT_MODE__
    if(!rank) {
      std::cout<<" Using "<<nlevels<<" Multigrid Levels."<<std::endl; 
      fflush(stdout);
      for(int i = 0; i < nlevels; i++) {
        std::cout<<" maxNpes for lev "<< i <<": "<< maxProcsForThisLevel[i] <<std::endl;
        fflush(stdout);
      }
    }
#endif

    PROF_SET_DA_STAGE3_END
#ifdef __PROF_WITH_BARRIER__
      MPI_Barrier(comm);
#endif
    PROF_SET_DA_STAGE4_BEGIN

      //Create activeComms for each level.
      //0 is the finest and (nlevels-1) is the coarsest
      MPI_Comm* activeComms = new MPI_Comm[nlevels];

    if(maxProcsForThisLevel[0] == npes) {
      activeComms[0] = comm;
    } else {
#ifndef __SILENT_MODE__
      if(!rank) {
        std::cout<<"splitting Comm in SetDA (using comm) "<<npes<<" -> "<<maxProcsForThisLevel[0]<<std::endl; 
      }
#endif
      par::splitCommUsingSplittingRank(maxProcsForThisLevel[0], activeComms, comm);
    }
    for(int lev = 1; lev < nlevels; lev++) {
      if(maxProcsForThisLevel[lev] == maxProcsForThisLevel[lev - 1]) {
        //coarse and fine grid use the same comm
        activeComms[lev] = activeComms[lev - 1];
      } else {
        assert(maxProcsForThisLevel[lev] < maxProcsForThisLevel[lev - 1]);
        if(maxProcsForThisLevel[lev] == activeNpesInCoarseBal[lev - 1]) {
          //The correct comm for this level was already created during coarsening
          activeComms[lev] = activeCommsInCoarseBal[lev - 1];
        } else if(maxProcsForThisLevel[lev] < activeNpesInCoarseBal[lev - 1]) {
          //The correct comm is a subset of the comm that was already created for this level
          if(maxProcsForThisLevel[lev - 1] < activeNpesInCoarseBal[lev - 1]) {
#ifndef __SILENT_MODE__
            if(!rank) {
              std::cout<<"splitting Comm in SetDA (using fine activeComm) "<<
                maxProcsForThisLevel[lev - 1]<<" -> "<<maxProcsForThisLevel[lev]<<std::endl; 
            }
#endif
            if(rank < maxProcsForThisLevel[lev - 1]) {
              par::splitCommUsingSplittingRank(maxProcsForThisLevel[lev],
                  (activeComms + lev), activeComms[lev - 1]);
            }
          } else {
#ifndef __SILENT_MODE__
            if(!rank) {
              std::cout<<"splitting Comm in SetDA (using coarseBalComm) "<<
                activeNpesInCoarseBal[lev - 1]<<" -> "<<maxProcsForThisLevel[lev]<<std::endl; 
            }
#endif
            if(rank < activeNpesInCoarseBal[lev - 1]) {
              par::splitCommUsingSplittingRank(maxProcsForThisLevel[lev],
                  (activeComms + lev), activeCommsInCoarseBal[lev - 1]);
            }
          }
        } else {
#ifndef __SILENT_MODE__
          if(!rank) {
            std::cout<<"splitting Comm in SetDA (using fine activeComm) "<<
              maxProcsForThisLevel[lev - 1]<<" -> "<<maxProcsForThisLevel[lev]<<std::endl; 
          }
#endif
          if(rank < maxProcsForThisLevel[lev - 1]) {
            par::splitCommUsingSplittingRank(maxProcsForThisLevel[lev],
                (activeComms + lev), activeComms[lev - 1]);
          }
        }
      }
    }//end for lev

    PROF_SET_DA_STAGE4_END
#ifdef __PROF_WITH_BARRIER__
      MPI_Barrier(comm);
#endif
    PROF_SET_DA_STAGE5_BEGIN

      //Distribute the octrees on maxProcsForThisLevel processors
      //finestOctree is not empty on any processor and so all processors
      //participate in this scatter
      if(maxProcsForThisLevel[0] < npes) {
        std::vector<ot::TreeNode> tmpOctree;
        if(rank < maxProcsForThisLevel[0]) {
          DendroIntL avgGlobalOctSize = (globalOctreeSizeForThisLevel[0]/maxProcsForThisLevel[0]);
          int extraOcts = (globalOctreeSizeForThisLevel[0] % maxProcsForThisLevel[0]);
          if(rank < extraOcts) {
            par::scatterValues<ot::TreeNode>(finestOctree, tmpOctree,
                (avgGlobalOctSize+1), comm);
          }else {
            par::scatterValues<ot::TreeNode>(finestOctree, tmpOctree,
                avgGlobalOctSize, comm);
          }
        }else {
          par::scatterValues<ot::TreeNode>(finestOctree, tmpOctree,
              0, comm);
        }
        finestOctree = tmpOctree;
        tmpOctree.clear();
      }

    //coarseOctrees[i] is only distributed on activeCommsInCoarseBal[i]
    //It is empty on all other processors
    //So by using "activeCommsInCoarseBal[i]" instead of "comm", we can skip idle
    //processors from participating in the scatter. 
    //If "activeCommsInCoarseBal[i]" is different from "comm", it would have
    //been created using the function "splitCommUsingSplittingRank" inside
    //"CoarsenOctree" or "BalanceOctree". Moreover, this will happen only when
    //the the total input size for that function is < 1000*npes. When this
    //happens the input would be redistributed to the first few processors.
    //"splitCommUsingSplittingRank" preserves the relative order of ranks in
    //the old and new comms. Hence, on all active processors the "rank" in
    //"comm" and in "activeCommsInCoarseBal[i]" will be the same. Similarly,
    //"activeComms" is also created above using "splitCommUsingSplittingRank",
    //hence the ranks in "activeComms" and "comm" will also be the same. 

    for(int i = 1; i < nlevels; i++) {
      if(activeNpesInCoarseBal[i-1] < maxProcsForThisLevel[i]) {
        //This case happens when the grain size using activeNpes is more than
        //1.5 times the grain size using maxProcs 
        //We need to expand the communicator
        //Both active and inactive processors need to participate here
        //But not all inactive processors need to participate. We have already
        //created the final activeComms above. We can use that here
        std::vector<ot::TreeNode> tmpOctree;
        if(rank < maxProcsForThisLevel[i]) {
          DendroIntL avgGlobalOctSize = (globalOctreeSizeForThisLevel[i]/maxProcsForThisLevel[i]);
          int extraOcts = (globalOctreeSizeForThisLevel[i] % maxProcsForThisLevel[i]);
          if(rank < extraOcts) {
            par::scatterValues<ot::TreeNode>(coarserOctrees[i-1], tmpOctree,
                (avgGlobalOctSize+1), activeComms[i]);
          }else {
            par::scatterValues<ot::TreeNode>(coarserOctrees[i-1], tmpOctree,
                avgGlobalOctSize, activeComms[i]);
          }
        }
        coarserOctrees[i-1] = tmpOctree;
        tmpOctree.clear();
      } else if(activeNpesInCoarseBal[i-1] > maxProcsForThisLevel[i]) {
        //We need to shrink the communicator
        //This case happens when we decide to use the same number of processors
        //as the immediate coarser level or too many extra boundaries were
        //discarded resulting in new empty processors
        //Inactive processors can be skipped here
        if(activeStatesInCoarseBal[i-1]) {
          std::vector<ot::TreeNode> tmpOctree;
          if(rank < maxProcsForThisLevel[i]) {
            DendroIntL avgGlobalOctSize = (globalOctreeSizeForThisLevel[i]/maxProcsForThisLevel[i]);
            int extraOcts = (globalOctreeSizeForThisLevel[i] % maxProcsForThisLevel[i]);
            if(rank < extraOcts) {
              par::scatterValues<ot::TreeNode>(coarserOctrees[i-1], tmpOctree,
                  (avgGlobalOctSize+1), activeCommsInCoarseBal[i-1]);
            }else {
              par::scatterValues<ot::TreeNode>(coarserOctrees[i-1], tmpOctree,
                  avgGlobalOctSize, activeCommsInCoarseBal[i-1]);
            }
          }else {
            par::scatterValues<ot::TreeNode>(coarserOctrees[i-1], tmpOctree,
                0, activeCommsInCoarseBal[i-1]);
          }
          coarserOctrees[i-1] = tmpOctree;
          tmpOctree.clear();
        }//end if active
      }//end if expanding or shrinking
    }//end for i

#ifdef __USE_PVT_DA_IN_MG__
    //A simple implementation for now. This just calls flagNodesType3 for all
    //levels. A smarter implementation would use the fact that regular coarse
    //grid nodes remain regular on all finer octrees 
    ot::markHangingNodesAtAllLevels(finestOctree, nlevels, coarserOctrees,
        activeComms, dim, maxDepth);
#endif

    if(activeStatesInCoarseBal) {
      delete [] activeStatesInCoarseBal;
      activeStatesInCoarseBal = NULL;
    }

    if(activeCommsInCoarseBal) {
      delete [] activeCommsInCoarseBal;
      activeCommsInCoarseBal = NULL;
    }

    if(activeNpesInCoarseBal) {
      delete [] activeNpesInCoarseBal;
      activeNpesInCoarseBal = NULL;
    }

    //nlevels has been corrected if necessary. Create the DAMG now...
    DAMG     *tmpDAMG = new DAMG[nlevels];  
    for (int i = 0; i < nlevels; i++) {
      tmpDAMG[i] = new _p_DAMG;    
      tmpDAMG[i]->nlevels  = nlevels - i;
      tmpDAMG[i]->totalLevels = nlevels;
      tmpDAMG[i]->comm     = comm;
      tmpDAMG[i]->user     = user;    
      tmpDAMG[i]->initialguess = NULL;
      tmpDAMG[i]->suppressedDOF = NULL; 
      tmpDAMG[i]->suppressedDOFaux = NULL; 
      tmpDAMG[i]->dof = dof;
      tmpDAMG[i]->da = NULL;
      tmpDAMG[i]->da_aux = NULL;
      tmpDAMG[i]->x = NULL;
      tmpDAMG[i]->b = NULL;
      tmpDAMG[i]->r = NULL;
      tmpDAMG[i]->J = NULL;
      tmpDAMG[i]->B = NULL;
      tmpDAMG[i]->R = NULL;
      tmpDAMG[i]->solve = NULL;		
      // KSP only 
      tmpDAMG[i]->ksp = NULL;             
      tmpDAMG[i]->rhs = NULL;
      // User had called stsDMMGSetKSP() and the matrices have been computed 
      tmpDAMG[i]->matricesset = PETSC_FALSE;
    }//end for i  
    *damg = tmpDAMG;

    PROF_SET_DA_STAGE5_END
#ifdef __PROF_WITH_BARRIER__
      MPI_Barrier(comm);
#endif
    PROF_SET_DA_STAGE6_BEGIN

      if(nlevels == 1) {
        //Single level only

#ifndef __USE_PVT_DA_IN_MG__
        tmpDAMG[0]->da = new DA(finestOctree, comm, activeComms[0], compressLut, NULL, NULL);
#else
        tmpDAMG[0]->da = new DA(1, finestOctree, comm, activeComms[0], compressLut, NULL, NULL);
#endif

        if(!rank) {
          int activeNpes = tmpDAMG[0]->da->getNpesActive();
          if(activeNpes < npes) {
            std::cout<<"WARNING: Some processors will be"<<
              " idle even on the finest level."<<
              " You might want to use fewer processors."<<std::endl;
            fflush(stdout);
          }
        }

        finestOctree.clear();

        //Free memory
        delete [] maxProcsForThisLevel;
        maxProcsForThisLevel = NULL;

        delete [] globalOctreeSizeForThisLevel;
        globalOctreeSizeForThisLevel = NULL;

        delete [] activeComms;
        activeComms = NULL;

        //If initial nlevels was 1, then coarserOctree would never have been
        //created. If nlevels was reset to 1 later then coarserOctrees would
        //have been created.
        if(coarserOctrees != NULL) {
          delete [] coarserOctrees;
          coarserOctrees = NULL;
        }

        PROF_SET_DA_STAGE6_END

          ierr = DAMGSetUp(tmpDAMG);CHKERRQ(ierr); 

        PROF_MG_SET_DA_END
      }//special case: single grid  

    //Start with the coarsest octree and mesh. 
    std::vector<ot::TreeNode>* blocksPtr = NULL;
    ot::DA* newDa = NULL;
    while(idxOfCoarsestLev >= 0) {
#ifdef __DEBUG_MG__
      if(!rank) {
        std::cout<<"Meshing lev: "<<idxOfCoarsestLev<<std::endl;
      }
#endif

      if(newDa == NULL) {
#ifndef __USE_PVT_DA_IN_MG__
        newDa = new DA(coarserOctrees[idxOfCoarsestLev], comm, 
            activeComms[idxOfCoarsestLev + 1], compressLut,
            blocksPtr, NULL);
#else
        newDa = new DA(1, coarserOctrees[idxOfCoarsestLev], comm, 
            activeComms[idxOfCoarsestLev + 1], compressLut,
            blocksPtr, NULL);
#endif
      }

      //The constructor will modify the input octree. So it's useless
      //afterwards.
      coarserOctrees[idxOfCoarsestLev].clear();
      tmpDAMG[nlevels-(2+idxOfCoarsestLev)]->da = newDa;
      newDa = NULL;

      bool procsNotUsedWell;
      if(!rank) {  
        int activeNpes = tmpDAMG[nlevels-(2+idxOfCoarsestLev)]->da->getNpesActive();
        if( (activeNpes < maxProcsForThisLevel[idxOfCoarsestLev+1]) ||
            (maxProcsForThisLevel[idxOfCoarsestLev+1] <
             maxProcsForThisLevel[idxOfCoarsestLev]) ) {
          procsNotUsedWell = true;
        } else {
          procsNotUsedWell = false;
        }
      }

      par::Mpi_Bcast<bool>(&procsNotUsedWell, 1, 0, comm);

      //For the mat-vecs, 
      //Communication cost: O(dependent)
      //Computation cost: O(own + pre-Ghost)
      //Heuristic: Give 75% weight to computation and 25% to communication
      //load =  0.75*computation + 0.25*communication  
      double minMaxLoad[2]; //0 for minLoad, 1 for maxLoad

      if(tmpDAMG[nlevels-(2+idxOfCoarsestLev)]->da->iAmActive()) {
        ot::DA* currDa = tmpDAMG[nlevels-(2+idxOfCoarsestLev)]->da;

        unsigned int localPreGhostSize = currDa->getIdxElementBegin();

        unsigned int localIndependentSize = currDa->getIndependentSize();

        unsigned int localOwnSize = currDa->getElementSize();

        unsigned int localDependentSize = (localOwnSize - localIndependentSize);

#ifdef __NO_GHOST_LOOP__
        double computationLoad = localOwnSize;
#else
        double computationLoad = (localOwnSize + localPreGhostSize);
#endif
        double communicationLoad = localDependentSize;

        double localLoad = ((0.75*computationLoad) + (0.25*communicationLoad));

        par::Mpi_Reduce<double>(&localLoad, &(minMaxLoad[0]), 1, 
            MPI_MIN, 0, currDa->getCommActive());

        par::Mpi_Reduce<double>(&localLoad, &(minMaxLoad[1]), 1, 
            MPI_MAX, 0, currDa->getCommActive());
      }

      par::Mpi_Bcast<double>(minMaxLoad, 2, 0, comm);

      //Give the finer grids a chance to be well load balanced and
      //be distributed on as many procs as possible.
      bool needAuxGrid = (procsNotUsedWell ||
          (minMaxLoad[1] > (minMaxLoad[0]*loadFac)));

      //Get the coarse grid partition. This will be needed whether or not the
      //aux_mesh is required
      if(blocksPtr == NULL) {
        blocksPtr = new std::vector<ot::TreeNode>;
        (*blocksPtr) = tmpDAMG[nlevels-(2+idxOfCoarsestLev)]->da->getBlocks();
      }

      if(needAuxGrid) {
        //Need DA_AUX for the immediate finer octree
        //Create DA_AUX for the immediate finer octree using the current
        //partition
        //We don't want to alter the actual fine octree. So we work with a copy
        //instead.
        std::vector<ot::TreeNode> fineOctreeCopy;
        if(idxOfCoarsestLev > 0) {
          fineOctreeCopy = coarserOctrees[idxOfCoarsestLev-1];
        }else {
          fineOctreeCopy = finestOctree;
        }

        //Need to redistribute fineOctreeCopy such that the input is only
        //distributed on the processors that are active on the coarse grid.
        bool iAmActive = tmpDAMG[nlevels-(2+idxOfCoarsestLev)]->da->iAmActive();
        MPI_Comm newComm = tmpDAMG[nlevels-(2+idxOfCoarsestLev)]->da->getCommActive();

        DendroIntL mySize = fineOctreeCopy.size();
        DendroIntL totalSize = globalOctreeSizeForThisLevel[idxOfCoarsestLev];   

        std::vector<ot::TreeNode> fineOctAfterPart;

        if(rank < maxProcsForThisLevel[idxOfCoarsestLev]) {
          if(iAmActive) {
            int newRank;
            int newNpes;
            MPI_Comm_rank(newComm, &newRank);
            MPI_Comm_size(newComm, &newNpes);
            DendroIntL avgSize = (totalSize / newNpes);
            int leftOvers = (totalSize % newNpes);
            if(newRank < leftOvers) {
              par::scatterValues<ot::TreeNode>(fineOctreeCopy, fineOctAfterPart,
                  (avgSize+1), activeComms[idxOfCoarsestLev]);
            }else {
              par::scatterValues<ot::TreeNode>(fineOctreeCopy, fineOctAfterPart,
                  avgSize, activeComms[idxOfCoarsestLev]);
            }
          }else {
            par::scatterValues<ot::TreeNode>(fineOctreeCopy, fineOctAfterPart,
                0, activeComms[idxOfCoarsestLev]);
          }//end if active on coarse grid
        }//end if active on fine grid

        fineOctreeCopy.clear();

        //This DA is aligned with the coarser grid
#ifndef __USE_PVT_DA_IN_MG__
        newDa = new DA(fineOctAfterPart, comm, newComm, compressLut, blocksPtr, NULL);
#else
        newDa = new DA(1, fineOctAfterPart, comm, newComm, compressLut, blocksPtr, NULL);
#endif

        fineOctAfterPart.clear();

        //Re-calculate imbalance. Check if aux mesh was really necessary
        //a-priori estimates could be misleading at times...

        if(!procsNotUsedWell) {
          if(newDa->iAmActive()) {
            ot::DA* currDa = newDa;

            unsigned int localPreGhostSize = currDa->getIdxElementBegin();

            unsigned int localOwnSize = currDa->getElementSize();

            unsigned int localIndependentSize = currDa->getIndependentSize();

            unsigned int localDependentSize = (localOwnSize - localIndependentSize);

#ifdef __NO_GHOST_LOOP__
            double computationLoad = localOwnSize;
#else
            double computationLoad = (localOwnSize + localPreGhostSize);
#endif
            double communicationLoad = localDependentSize;

            double localLoad = ((0.75*computationLoad) + (0.25*communicationLoad));

            par::Mpi_Reduce<double>(&localLoad, &minMaxLoad[0], 1, 
                MPI_MIN, 0, currDa->getCommActive());

            par::Mpi_Reduce<double>(&localLoad, &minMaxLoad[1], 1, 
                MPI_MAX, 0, currDa->getCommActive());
          }

          par::Mpi_Bcast<double>(minMaxLoad, 2, 0, comm);

          needAuxGrid = (minMaxLoad[1] > (minMaxLoad[0]*loadFac));
        } else {
          needAuxGrid = true;
        }

        if(needAuxGrid) {
          tmpDAMG[nlevels-(1+idxOfCoarsestLev)]->da_aux = newDa;
          newDa = NULL;

          //We don't want to use the same partition for the next finer level.
          blocksPtr->clear();
          delete blocksPtr;
          blocksPtr = NULL;
        }
      }//end check for aux_mesh 

      idxOfCoarsestLev--;
    }//end loop to mesh coarser grids

    //Free the memory for the pointers to the coarser octrees.
    delete [] coarserOctrees;
    coarserOctrees = NULL;

    delete [] globalOctreeSizeForThisLevel;
    globalOctreeSizeForThisLevel = NULL;

    delete [] maxProcsForThisLevel;
    maxProcsForThisLevel = NULL;

    //mesh finest level here
    if(newDa == NULL) {
#ifndef __USE_PVT_DA_IN_MG__
      newDa = new DA(finestOctree, comm, activeComms[0], compressLut, blocksPtr, NULL);
#else
      newDa = new DA(1, finestOctree, comm, activeComms[0], compressLut, blocksPtr, NULL);
#endif
    }

    tmpDAMG[nlevels-1]->da = newDa;
    newDa = NULL;

    finestOctree.clear();

    delete [] activeComms;
    activeComms = NULL;

    if(!rank) {
      int activeNpes = tmpDAMG[nlevels-1]->da->getNpesActive();
      if(activeNpes < npes) {
        std::cout<<"WARNING: Some processors will be idle"
          <<" even on the finest level."
          <<" You might want to use fewer processors."<<std::endl;
        fflush(stdout);
      }
    }

    if(blocksPtr != NULL) {
      blocksPtr->clear();
      delete blocksPtr;
      blocksPtr = NULL;
    }

    PROF_SET_DA_STAGE6_END

      ierr = DAMGSetUp(tmpDAMG); CHKERRQ(ierr); 

    PROF_MG_SET_DA_END

  }//end fn.

  PetscErrorCode DAMGSetUp(DAMG* damg)
  {
    PROF_MG_SETUP_BEGIN

      int       i,nlevels = damg[0]->nlevels;


#ifdef __DEBUG_MG__
    int rank;
    MPI_Comm_rank(damg[0]->comm,&rank);
    if(!rank) {
      std::cout<<"In DAMGSetUp."<<std::endl;
      fflush(stdout);
    }	
#endif

    /* Create work vectors and matrix for each level */
    damg[0]->da->createVector(damg[0]->x, false, false, damg[0]->dof);
    VecDuplicate(damg[0]->x,&damg[0]->b); 
    VecDuplicate(damg[0]->x,&damg[0]->r); 

#ifdef __DEBUG_MG__
    assert(damg[0]->da_aux == NULL);
#endif

    for (i = 1; i < nlevels; i++) {
      damg[i]->da->createVector(damg[i]->x, false, false, damg[i]->dof);
      VecDuplicate(damg[i]->x,&damg[i]->b); 
      VecDuplicate(damg[i]->x,&damg[i]->r); 

      /* Create interpolation/restriction between levels. */  
      //nlevels must be atleast 2 for Restriction/Interpolation  
      //i is fine, i-1 is coarse
      if(damg[i]->da_aux == NULL) {
        //Type-1. Assume Coarse and Fine are aligned. Do not use aux.
#ifdef __DEBUG_MG__
        if(!rank) {
          std::cout<<"Lev: "<<i<<" R-type: 1"<<std::endl;
          fflush(stdout);
        }
#endif
        createInterpolationType1(damg[i-1],damg[i],&damg[i]->R); 
      }else {
        //Type-2. Use Aux. Coarse and Fine are not aligned.
#ifdef __DEBUG_MG__
        if(!rank) {
          std::cout<<"Lev: "<<i<<" R-type: 2"<<std::endl;
          fflush(stdout);
        }
#endif
        createInterpolationType2(damg[i-1],damg[i],&damg[i]->R); 
      }
    }//end for i

    PROF_MG_SETUP_END
  }

  /*Matrix-Free Intergrid Transfer Operators 
Type1: Assume Coarse and Fine are aligned. Do not use aux.
Type2: Use Aux. Coarse and Fine are not aligned.
*/
  PetscErrorCode createInterpolationType1(DAMG damgc, DAMG damgf, Mat *R)
  {
    PROF_MG_CREATE_RP1_BEGIN

      ot::DA*   dac = damgc->da; 
    ot::DA*   daf = damgf->da;	
    unsigned int dof = damgc->dof;
    MPI_Comm comm = damgc->comm;
    unsigned int  mc,mf;
    //The size this processor owns ( without ghosts).
    mc=dac->getNodeSize();
    mf=daf->getNodeSize();
    TransferOpData* data = new TransferOpData;

    unsigned int localCoarseIndSize =  dac->getIndependentSize();

    if(dac->iAmActive()) {
      par::Mpi_Allreduce<unsigned int>(&localCoarseIndSize, &(data->minIndependentSize), 1,
          MPI_MIN, dac->getCommActive());
    } else {
      data->minIndependentSize = 0;
    }

    data->dac = dac;
    data->daf = daf;
    data->suppressedDOFc = damgc->suppressedDOF;
    data->suppressedDOFf = damgf->suppressedDOF;
    data->comm = comm;
    data->dof = dof;
    data->fineTouchedFlags = new std::vector<ot::FineTouchedStatus>;

    data->tmp = NULL;
    data->addRtmp = NULL;
    data->addPtmp = NULL;

    data->sendSzP = NULL;
    data->sendOffP = NULL;
    data->recvSzP = NULL;
    data->recvOffP = NULL;
    data->sendSzR = NULL;
    data->sendOffR = NULL;
    data->recvSzR = NULL;
    data->recvOffR = NULL;

    iC(MatCreateShell(comm, mc*dof ,mf*dof,PETSC_DETERMINE,PETSC_DETERMINE, (void*)(data), R));  
    iC(MatShellSetOperation(*R,MATOP_MULT_TRANSPOSE,(void(*)(void))prolongMatVecType1));
    iC(MatShellSetOperation(*R,MATOP_MULT,(void(*)(void))restrictMatVecType1));
    iC(MatShellSetOperation(*R,MATOP_MULT_ADD,(void(*)(void))addRestrictMatVec));
    iC(MatShellSetOperation(*R,MATOP_MULT_TRANSPOSE_ADD,(void(*)(void))addProlongMatVec));
    iC(MatShellSetOperation(*R ,MATOP_DESTROY, (void(*)(void)) rpMatDestroy));

    //dummy call to restrict-type-1
    data->daf->createVector<ot::FineTouchedStatus >(*(data->fineTouchedFlags), false, false, 1);
    dummyRestrictMatVecType1(data);

    PROF_MG_CREATE_RP1_END
  }

  PetscErrorCode  createInterpolationType2(DAMG damgc, DAMG damgf, Mat *R)	
  {
    PROF_MG_CREATE_RP2_BEGIN

      ot::DA*   dac = damgc->da; 
    ot::DA*   daf = damgf->da;	
    ot::DA*  daf_aux = damgf->da_aux;
    unsigned int dof = damgc->dof;
    MPI_Comm comm = damgc->comm;
    unsigned int  mc,mf;
    //The size this processor owns ( without ghosts).
    mc=dac->getNodeSize();
    mf=daf->getNodeSize();
    TransferOpData* data = new TransferOpData;

    unsigned int localCoarseIndSize =  dac->getIndependentSize();

    if(dac->iAmActive()) {
      par::Mpi_Allreduce<unsigned int>(&localCoarseIndSize, &(data->minIndependentSize), 1,
          MPI_MIN, dac->getCommActive());
    } else {
      data->minIndependentSize = 0;
    }

    data->dac = dac;
    data->daf = daf_aux;	
    data->suppressedDOFc = damgc->suppressedDOF;
    data->suppressedDOFf = damgf->suppressedDOFaux;
    data->comm = comm;
    data->dof = dof;
    data->fineTouchedFlags = new std::vector<ot::FineTouchedStatus>;
    daf_aux->createVector(data->tmp, false, false, dof);
    data->addRtmp = NULL;
    data->addPtmp = NULL;

    data->sendSzP = NULL;
    data->sendOffP = NULL;
    data->recvSzP = NULL;
    data->recvOffP = NULL;
    data->sendSzR = NULL;
    data->sendOffR = NULL;
    data->recvSzR = NULL;
    data->recvOffR = NULL;

    iC(MatCreateShell(comm, mc*dof ,mf*dof,PETSC_DETERMINE,PETSC_DETERMINE, (void*)(data), R));  
    iC(MatShellSetOperation(*R,MATOP_MULT_TRANSPOSE,(void(*)(void))prolongMatVecType2));
    iC(MatShellSetOperation(*R,MATOP_MULT,(void(*)(void))restrictMatVecType2));
    iC(MatShellSetOperation(*R,MATOP_MULT_TRANSPOSE_ADD,(void(*)(void))addProlongMatVec));
    iC(MatShellSetOperation(*R,MATOP_MULT_ADD,(void(*)(void))addRestrictMatVec));
    iC(MatShellSetOperation(*R ,MATOP_DESTROY, (void(*)(void)) rpMatDestroy));

    //dummy call to restrict-type-1.
    //Note type2 is not called. Since we don't need to scatter values for the
    //dummy loop. In fact, it will be wrong to use scatter values since it will
    //create the scatter contexts the first time it is called.
    data->daf->createVector<ot::FineTouchedStatus >(*(data->fineTouchedFlags), false, false, 1);
    dummyRestrictMatVecType1(data);

    PROF_MG_CREATE_RP2_END
  }

  PetscErrorCode  rpMatDestroy(Mat R) {
    TransferOpData *data;
    PetscFunctionBegin;
    iC(MatShellGetContext( R, (void **)&data));
    if(data) {
      if(data->fineTouchedFlags) { 
        data->fineTouchedFlags->clear(); 
        delete data->fineTouchedFlags;
        data->fineTouchedFlags = NULL;
      }
      if(data->tmp) { 
        iC(VecDestroy(data->tmp));
        data->tmp = NULL;
      }  
      if(data->addRtmp) {
        VecDestroy(data->addRtmp);
        data->addRtmp = NULL;
      }
      if(data->addPtmp) {
        VecDestroy(data->addPtmp);
        data->addPtmp = NULL;
      }
      if(data->sendSzP) {
        delete [] data->sendSzP;
        data->sendSzP = NULL;
      }
      if(data->sendOffP) {
        delete [] data->sendOffP;
        data->sendOffP = NULL;
      }
      if(data->recvSzP) {
        delete [] data->recvSzP;
        data->recvSzP = NULL;
      }
      if(data->recvOffP) {
        delete [] data->recvOffP;
        data->recvOffP = NULL;
      }
      if(data->sendSzR) {
        delete [] data->sendSzR;
        data->sendSzR = NULL;
      }
      if(data->sendOffR) {
        delete [] data->sendOffR;
        data->sendOffR = NULL;
      }
      if(data->recvSzR) {
        delete [] data->recvSzR;
        data->recvSzR = NULL;
      }
      if(data->recvOffR) {
        delete [] data->recvOffR;
        data->recvOffR = NULL;
      }
      delete data;
      data = NULL;
    }

    PetscFunctionReturn(0);
  }

  void PrintDAMG(DAMG* damg) {
    int  nlevels = damg[0]->nlevels;
    int rank;
    MPI_Comm_rank(damg[0]->comm,&rank);
    for (int i = 0; i < nlevels; i++) {
      MPI_Comm activeComm = damg[i]->da->getCommActive();
      MPI_Comm activeAuxComm;

      if(damg[i]->da_aux) {
        activeAuxComm = damg[i]->da_aux->getCommActive();
      }

      DendroIntL globalBlocksSize;
      DendroIntL globalBlocksSizeAux;
      DendroIntL localBlocksSize = damg[i]->da->getNumBlocks();
      DendroIntL localBlocksSizeAux;

      if(damg[i]->da_aux) {
        localBlocksSizeAux = damg[i]->da_aux->getNumBlocks();
      }

      DendroIntL elemSizeLocal = damg[i]->da->getElementSize();
      DendroIntL ghostedElemSizeLocal = damg[i]->da->getGhostedElementSize();
      DendroIntL elemSizeLocalAux;	
      DendroIntL ghostedElemSizeLocalAux;	

      if(damg[i]->da_aux) {
        elemSizeLocalAux = damg[i]->da_aux->getElementSize();
        ghostedElemSizeLocalAux = damg[i]->da_aux->getGhostedElementSize();
      }	

      DendroIntL elemSizeGlobal;
      DendroIntL minElemSize;
      DendroIntL maxElemSize;
      DendroIntL minGhostedElemSize;
      DendroIntL maxGhostedElemSize;

      DendroIntL elemSizeGlobalAux;
      DendroIntL minElemSizeAux;
      DendroIntL maxElemSizeAux;
      DendroIntL minGhostedElemSizeAux;
      DendroIntL maxGhostedElemSizeAux;

      if(damg[i]->da->iAmActive()) {
        par::Mpi_Reduce<DendroIntL>(&localBlocksSize, &globalBlocksSize, 1, 
            MPI_SUM, 0, activeComm);

        par::Mpi_Reduce<DendroIntL>(&elemSizeLocal, &elemSizeGlobal, 1, 
            MPI_SUM, 0, activeComm);

        par::Mpi_Reduce<DendroIntL>(&elemSizeLocal, &minElemSize, 1, 
            MPI_MIN, 0, activeComm);

        par::Mpi_Reduce<DendroIntL>(&elemSizeLocal, &maxElemSize, 1, 
            MPI_MAX, 0, activeComm);

        par::Mpi_Reduce<DendroIntL>(&ghostedElemSizeLocal, &minGhostedElemSize, 1, 
            MPI_MIN, 0, activeComm);

        par::Mpi_Reduce<DendroIntL>(&ghostedElemSizeLocal, &maxGhostedElemSize, 1, 
            MPI_MAX, 0, activeComm);
      }

      if(damg[i]->da_aux) {
        if(damg[i]->da_aux->iAmActive()) {
          par::Mpi_Reduce<DendroIntL>(&localBlocksSizeAux, &globalBlocksSizeAux, 1, 
              MPI_SUM, 0, activeAuxComm);

          par::Mpi_Reduce<DendroIntL>(&elemSizeLocalAux, &elemSizeGlobalAux, 1, 
              MPI_SUM, 0, activeAuxComm);

          par::Mpi_Reduce<DendroIntL>(&elemSizeLocalAux, &minElemSizeAux, 1,  
              MPI_MIN, 0, activeAuxComm);

          par::Mpi_Reduce<DendroIntL>(&elemSizeLocalAux, &maxElemSizeAux, 1, 
              MPI_MAX, 0, activeAuxComm);

          par::Mpi_Reduce<DendroIntL>(&ghostedElemSizeLocalAux, &minGhostedElemSizeAux, 1, 
              MPI_MIN, 0, activeAuxComm);

          par::Mpi_Reduce<DendroIntL>(&ghostedElemSizeLocalAux, &maxGhostedElemSizeAux, 1, 
              MPI_MAX, 0, activeAuxComm);
        }
      }	

      DendroIntL nodeSizeLocal = damg[i]->da->getNodeSize();
      DendroIntL nodeSizeLocalAux;

      if(damg[i]->da_aux) {
        nodeSizeLocalAux = damg[i]->da_aux->getNodeSize();
      }

      DendroIntL nodeSizeGlobal;
      DendroIntL minNodeSize;
      DendroIntL maxNodeSize;

      if(damg[i]->da->iAmActive()) {
        par::Mpi_Reduce<DendroIntL>(&nodeSizeLocal, &nodeSizeGlobal, 1, 
            MPI_SUM, 0, activeComm);

        par::Mpi_Reduce<DendroIntL>(&nodeSizeLocal, &minNodeSize, 1, 
            MPI_MIN, 0, activeComm);

        par::Mpi_Reduce<DendroIntL>(&nodeSizeLocal, &maxNodeSize, 1, 
            MPI_MAX, 0, activeComm);
      }

      DendroIntL nodeSizeGlobalAux;
      DendroIntL minNodeSizeAux;
      DendroIntL maxNodeSizeAux;

      if(damg[i]->da_aux) {
        if(damg[i]->da_aux->iAmActive()) {
          par::Mpi_Reduce<DendroIntL>(&nodeSizeLocalAux, &nodeSizeGlobalAux, 1, 
              MPI_SUM, 0, activeAuxComm);

          par::Mpi_Reduce<DendroIntL>(&nodeSizeLocalAux, &minNodeSizeAux, 1, 
              MPI_MIN, 0, activeAuxComm);

          par::Mpi_Reduce<DendroIntL>(&nodeSizeLocalAux, &maxNodeSizeAux, 1, 
              MPI_MAX, 0, activeAuxComm);
        }
      }

      unsigned int globalMinLev = getGlobalMinLevel(damg[i]->da);
      unsigned int globalMaxLev = getGlobalMaxLevel(damg[i]->da);

      unsigned int globalMinLevAux;
      if(damg[i]->da_aux) {	
        globalMinLevAux = getGlobalMinLevel(damg[i]->da_aux);
      }
      unsigned int globalMaxLevAux;
      if(damg[i]->da_aux) {
        globalMaxLevAux = getGlobalMaxLevel(damg[i]->da_aux);
      }

      int activeNpes = damg[i]->da->getNpesActive();
      int activeNpesAux;
      if(damg[i]->da_aux) {
        activeNpesAux = damg[i]->da_aux->getNpesActive();
      }
      if(!rank) {
        std::cout<<"At Level "<<i<<" Elem size: "<<elemSizeGlobal
          <<", GlobalBlocksSize: "<<globalBlocksSize
          <<", node size: "<<nodeSizeGlobal
          <<", minElemSize: "<<minElemSize
          <<", maxElemSize: "<<maxElemSize
          <<", minGhostedElemSize: "<<minGhostedElemSize
          <<", maxGhostedElemSize: "<<maxGhostedElemSize
          <<", minNodeSize: "<<minNodeSize
          <<", maxNodeSize: "<<maxNodeSize
          <<", min Lev: "<<globalMinLev
          <<", max Lev: "<<globalMaxLev
          <<", Active Npes: "<<activeNpes<<std::endl; 
        fflush(stdout);
        if(damg[i]->da_aux) {
          std::cout<<"At AUX Level "<<i<<" Elem size: "<<elemSizeGlobalAux
            <<", GlobalBlocksSize: "<<globalBlocksSizeAux
            <<", node size: "<<nodeSizeGlobalAux
            <<", minElemSize: "<<minElemSizeAux
            <<", maxElemSize: "<<maxElemSizeAux
            <<", minGhostedElemSize: "<<minGhostedElemSizeAux
            <<", maxGhostedElemSize: "<<maxGhostedElemSizeAux
            <<", minNodeSize: "<<minNodeSizeAux
            <<", maxNodeSize: "<<maxNodeSizeAux
            <<", min Lev: "<<globalMinLevAux
            <<", max Lev: "<<globalMaxLevAux
            <<", Active Npes: "<<activeNpesAux<<std::endl; 
          fflush(stdout);
        }
      }//end if p0
    }//end for i
  }//end function

  //Private Functions...

  PetscErrorCode PC_KSP_Shell_SetUp(void* ctx) {

    PetscFunctionBegin;

    PC_KSP_Shell* data = static_cast<PC_KSP_Shell*>(ctx); 

    MPI_Comm commActive = data->commActive;
    bool iAmActive = data->iAmActive;

    //This points to the shell itself
    PC pc = data->pc;

    Mat Amat;
    PCGetOperators(pc, &Amat, NULL, NULL);      
    PetscTruth isshell;
    PetscTypeCompare((PetscObject)Amat, MATSHELL, &isshell);

    if(!isshell) {
      SETERRQ(PETSC_ERR_SUP, " Expected a MATSHELL.");
    }

    //Create ksp_private, rhs_private, sol_private,
    if(iAmActive) {
      Mat Amat_private;
      Mat Pmat_private;
      MatStructure pFlag;

      if(getPrivateMatricesForKSP_Shell) {
        (*getPrivateMatricesForKSP_Shell)(Amat, &Amat_private,
            &Pmat_private, &pFlag);
      } else {
        SETERRQ(PETSC_ERR_USER,
            " Expected function to be set:\
            getPrivateMatricesForKSP_Shell");
      }

      //Sanity Checks
      assert(Amat_private != NULL);

      PetscInt globalRowSize;
      PetscInt globalColSize;
      PetscInt localRowSize;
      PetscInt localColSize;

      PetscInt globalRowSizePrivate;
      PetscInt globalColSizePrivate;
      PetscInt localRowSizePrivate;
      PetscInt localColSizePrivate;

      MatGetSize(Amat, &globalRowSize, &globalColSize);
      MatGetSize(Amat_private, &globalRowSizePrivate, &globalColSizePrivate);

      MatGetLocalSize(Amat, &localRowSize, &localColSize);
      MatGetLocalSize(Amat_private, &localRowSizePrivate, &localColSizePrivate);

      assert(globalRowSize == globalRowSizePrivate);
      assert(globalColSize == globalColSizePrivate);

      assert(localRowSize == localRowSizePrivate);
      assert(localColSize == localColSizePrivate);

      MPI_Comm privateComm;
      PetscObjectGetComm((PetscObject)Amat_private, &privateComm);

      int commCompareResult;

      //PetscHeaderCreate_Private duplicates comm
      //instead of a simple copy. So, the comms will not exactly be
      //identical. But, they will be equivalent  
      MPI_Comm_compare(privateComm, commActive, &commCompareResult);

      assert( (commCompareResult == MPI_CONGRUENT) || (commCompareResult == MPI_IDENT));

      if(pc->setupcalled == 0) {
        assert(data->ksp_private == NULL);
        assert(data->rhs_private == NULL);
        assert(data->sol_private == NULL);
      } else {
        assert(data->ksp_private != NULL);
        assert(data->rhs_private != NULL);
        assert(data->sol_private != NULL);
      }

      if(pc->setupcalled == 0) {
        KSPCreate(commActive, &(data->ksp_private));

        const char *prefix;
        PCGetOptionsPrefix(pc, &prefix);

        //These functions also set the correct prefix for the inner pc 
        KSPSetOptionsPrefix(data->ksp_private, prefix);
        KSPAppendOptionsPrefix(data->ksp_private, "private_");

        //Default Types for KSP and PC
        KSPSetType(data->ksp_private, KSPPREONLY);

        PC privatePC;
        KSPGetPC(data->ksp_private, &privatePC);
        PCSetType(privatePC, PCLU);

        //The command line options get higher precedence.
        //This also calls PCSetFromOptions for the private pc internally
        KSPSetFromOptions(data->ksp_private);  

        MatGetVecs(Amat_private, &(data->sol_private), &(data->rhs_private));
      }

      KSPSetOperators(data->ksp_private, Amat_private, Pmat_private, pFlag);

    } else {
      data->sol_private = NULL;
      data->rhs_private = NULL;
      data->ksp_private = NULL;
    }//end if active

    PetscFunctionReturn(0);
  }

  PetscErrorCode PC_KSP_Shell_Destroy(void* ctx) {

    PetscFunctionBegin;

    PC_KSP_Shell* data = static_cast<PC_KSP_Shell*>(ctx); 

    if(data) {
      if(data->ksp_private) {
        KSPDestroy(data->ksp_private);
        data->ksp_private = NULL;
      }

      if(data->rhs_private) {
        VecDestroy(data->rhs_private);
        data->rhs_private = NULL;
      }

      if(data->sol_private) {
        VecDestroy(data->sol_private);
        data->sol_private = NULL;
      }

      delete data;
      data = NULL;
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode PC_KSP_Shell_Apply(void* ctx, Vec rhs, Vec sol) {

    PetscFunctionBegin;

    PC_KSP_Shell* data = static_cast<PC_KSP_Shell*>(ctx); 

    if(data->iAmActive) {      
      PetscScalar* rhsArray;
      PetscScalar* solArray;

      //There are no copies and no mallocs involved.

      VecGetArray(rhs, &rhsArray);
      VecGetArray(sol, &solArray);

      VecPlaceArray(data->rhs_private, rhsArray);
      VecPlaceArray(data->sol_private, solArray);

      KSPSolve(data->ksp_private, data->rhs_private, data->sol_private);

      VecResetArray(data->rhs_private);
      VecResetArray(data->sol_private);

      VecRestoreArray(rhs, &rhsArray);
      VecRestoreArray(sol, &solArray);
    }

    PetscFunctionReturn(0);
  }

  int scatterValues(Vec in, Vec out, PetscInt inSz, PetscInt outSz,
      int *& sendSz, int *& sendOff, int *& recvSz, int *& recvOff, MPI_Comm comm ) {

    PROF_MG_SCATTER_BEGIN 

      bool isNew = (sendSz == NULL);

    if(isNew) {
      int rank, npes;

      MPI_Comm_size(comm, &npes);
      MPI_Comm_rank(comm, &rank);

      MPI_Request request;
      MPI_Status status;

      PetscInt inSize, outSize; 

      VecGetLocalSize(in, &inSize);
      VecGetLocalSize(out, &outSize);

      assert(inSize == inSz);
      assert(outSize == outSz);

      PetscInt  off1 = 0, off2 = 0;
      PetscInt* scnIn = NULL;
      if(inSz) {
        scnIn = new PetscInt [inSz]; 
      }

      // perform a local scan first ...
      PetscInt zero=0;
      if(inSz) {
        scnIn[0] = 1;
        for (PetscInt i = 1; i < inSz; i++) {
          scnIn[i] = scnIn[i-1] + 1;
        }//end for
        // now scan with the final members of 
        par::Mpi_Scan<PetscInt>(scnIn+inSz-1, &off1, 1, MPI_SUM, comm ); 
      } else{
        par::Mpi_Scan<PetscInt>(&zero, &off1, 1, MPI_SUM, comm ); 
      }

      // communicate the offsets ...
      if (rank < (npes-1)){
        par::Mpi_Issend<PetscInt>( &off1, 1, rank+1, 0, comm, &request );
      }
      if (rank){
        par::Mpi_Recv<PetscInt>( &off2, 1, rank-1, 0, comm, &status );
      } else{
        off2 = 0; 
      }

      // add offset to local array
      for (PetscInt i = 0; i < inSz; i++) {
        scnIn[i] = scnIn[i] + off2;  // This has the global scan results now ...
      }//end for

      //Gather Scan of outCnts
      PetscInt* outCnts;
      outCnts = new PetscInt[npes];

      if(rank < (npes-1)) {
        MPI_Status statusWait;
        MPI_Wait(&request, &statusWait);
      }

      if( outSz ) {
        par::Mpi_Scan<PetscInt>( &outSz, &off1, 1, MPI_SUM, comm ); 
      }else {
        par::Mpi_Scan<PetscInt>( &zero, &off1, 1, MPI_SUM, comm ); 
      }

      par::Mpi_Allgather<PetscInt>( &off1, outCnts, 1, comm);

#ifdef __DEBUG_MG__
      assert(sendSz == NULL);
      assert(recvSz == NULL);
      assert(sendOff == NULL);
      assert(recvOff == NULL);
#endif

      sendSz = new int [npes];
      recvSz = new int [npes];
      sendOff = new int [npes];
      recvOff = new int [npes];

      // compute the partition offsets and sizes so that All2Allv can be performed.
      // initialize ...
      for (int i = 0; i < npes; i++) {
        sendSz[i] = 0;
      }

      //The Heart of the algorithm....
      //scnIn and outCnts are both sorted 
      PetscInt inCnt = 0;
      int pCnt = 0;
      while( (inCnt < inSz) && (pCnt < npes) ) {
        if( scnIn[inCnt] <= outCnts[pCnt]  ) {
          sendSz[pCnt]++;
          inCnt++;
        }else {
          pCnt++;
        }
      }

      // communicate with other procs how many you shall be sending and get how
      // many to recieve from whom.
      par::Mpi_Alltoall<int>(sendSz, recvSz, 1, comm);

      PetscInt nn=0; // new value of nlSize, ie the local nodes.
      for (int i = 0; i < npes; i++) {
        nn += recvSz[i];
      }

      // compute offsets ...
      sendOff[0] = 0;
      recvOff[0] = 0;
      for (int i=1; i<npes; i++) {
        sendOff[i] = sendOff[i-1] + sendSz[i-1];
        recvOff[i] = recvOff[i-1] + recvSz[i-1];
      }

      assert(nn == outSz);
      // clean up...
      if(scnIn) {
        delete [] scnIn;
      }
      delete [] outCnts;
    }//end if isNew

    // perform All2All  ... 
    PetscScalar * inArr;
    PetscScalar * outArr;

    VecGetArray(in,&inArr);
    VecGetArray(out,&outArr);

    par::Mpi_Alltoallv_sparse<PetscScalar>(inArr, sendSz, sendOff, 
        outArr, recvSz, recvOff, comm);

    VecRestoreArray(in,&inArr);
    VecRestoreArray(out,&outArr);

    PROF_MG_SCATTER_END 
  }//end function

}//end namespace



