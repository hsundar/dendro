
/**
  @file VtxMaps.C
  @brief Load stencils for Restriction matvec
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  */

#include <cstdio>
#include <iostream>
#include <cassert>

namespace ot {

  int readVtxMaps(unsigned short ****&map1, unsigned short *****&map2,
      unsigned short *****&map3, unsigned short ******&map4) {

    FILE* infile;
    char fname[100];
    sprintf(fname,"vtxMap.inp");
    infile = fopen(fname,"r");
    if(!infile) {
      std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
      assert(false);
    }

    typedef unsigned short***** us5Ptr;
    typedef unsigned short**** us4Ptr;
    typedef unsigned short*** us3Ptr;
    typedef unsigned short** us2Ptr;
    typedef unsigned short* usPtr;

    map1 = new us3Ptr[8];
    map2 = new us4Ptr[8];
    for(int fineElemNum = 0; fineElemNum < 8;fineElemNum++) {
      map1[fineElemNum] = new us2Ptr[8];
      for(int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
        map1[fineElemNum][cNumCoarse] = new usPtr[18];
        for(int eTypeCoarse = 0; eTypeCoarse < 18; eTypeCoarse++) {
          map1[fineElemNum][cNumCoarse][eTypeCoarse] = new unsigned short[8];
          for(int vCtr = 0; vCtr < 8; vCtr++) {
            fscanf(infile,"%hu",
                &(map1[fineElemNum][cNumCoarse][eTypeCoarse][vCtr]));
          }
        }
      }
      map2[fineElemNum] = new us3Ptr[8];
      for(int cNumFine = 0; cNumFine < 8; cNumFine++) {
        map2[fineElemNum][cNumFine] = new us2Ptr[8];
        for(int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
          map2[fineElemNum][cNumFine][cNumCoarse] = new usPtr[18];
          for(int eTypeCoarse = 0; eTypeCoarse < 18; eTypeCoarse++) {
            map2[fineElemNum][cNumFine][cNumCoarse][eTypeCoarse] = new unsigned short[8];
            for(int vCtr = 0; vCtr < 8; vCtr++) {
              fscanf(infile,"%hu",
                  &(map2[fineElemNum][cNumFine][cNumCoarse][eTypeCoarse][vCtr]));
            }
          }
        }
      }
    }

    map3 = new us4Ptr[7];
    map4 = new us5Ptr[7];
    for(int fineElemNum = 0; fineElemNum < 7; fineElemNum++) { 
      map3[fineElemNum] = new us3Ptr[2];
      map4[fineElemNum] = new us4Ptr[2];
      for(int scalingCtr = 0; scalingCtr < 2; scalingCtr++) {
        map3[fineElemNum][scalingCtr] = new us2Ptr[8];
        for(int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
          map3[fineElemNum][scalingCtr][cNumCoarse] = new usPtr[18];
          for(int eTypeCoarse = 0; eTypeCoarse < 18; eTypeCoarse++) {
            map3[fineElemNum][scalingCtr][cNumCoarse][eTypeCoarse] =new unsigned short[8];
            for(int vCtr = 0; vCtr < 8; vCtr++) {
              fscanf(infile,"%hu",
                  &(map3[fineElemNum][scalingCtr][cNumCoarse][eTypeCoarse][vCtr]));
            }
          }
        }
        map4[fineElemNum][scalingCtr] = new us3Ptr[8]; 
        for(int cNumFine = 0; cNumFine < 8; cNumFine++) {
          map4[fineElemNum][scalingCtr][cNumFine] = new us2Ptr[8];
          for(int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
            map4[fineElemNum][scalingCtr][cNumFine][cNumCoarse] = new usPtr[18];
            for(int eTypeCoarse = 0; eTypeCoarse < 18; eTypeCoarse++) {
              map4[fineElemNum][scalingCtr][cNumFine][cNumCoarse][eTypeCoarse] = new unsigned short[8];
              for(int vCtr = 0; vCtr < 8; vCtr++) {
                fscanf(infile,"%hu",
                    &(map4[fineElemNum][scalingCtr][cNumFine][cNumCoarse][eTypeCoarse][vCtr]));
              }
            }
          }
        }
      }
    }

    fclose(infile);
    return 1;
  }//end of function

  int IreadVtxMaps(unsigned short ****&map1, unsigned short *****&map2,
      unsigned short *****&map3, unsigned short ******&map4, int rank) {

    FILE* infile;
    char fname[250];
    sprintf(fname,"vtxMap_%d.inp", rank);
    infile = fopen(fname,"r");
    if(!infile) {
      std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
      assert(false);
    }

    typedef unsigned short***** us5Ptr;
    typedef unsigned short**** us4Ptr;
    typedef unsigned short*** us3Ptr;
    typedef unsigned short** us2Ptr;
    typedef unsigned short* usPtr;

    map1 = new us3Ptr[8];
    map2 = new us4Ptr[8];
    for(int fineElemNum = 0; fineElemNum < 8;fineElemNum++) {
      map1[fineElemNum] = new us2Ptr[8];
      for(int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
        map1[fineElemNum][cNumCoarse] = new usPtr[18];
        for(int eTypeCoarse = 0; eTypeCoarse < 18; eTypeCoarse++) {
          map1[fineElemNum][cNumCoarse][eTypeCoarse] = new unsigned short[8];
          for(int vCtr = 0; vCtr < 8; vCtr++) {
            fscanf(infile,"%hu",
                &(map1[fineElemNum][cNumCoarse][eTypeCoarse][vCtr]));
          }
        }
      }
      map2[fineElemNum] = new us3Ptr[8];
      for(int cNumFine = 0; cNumFine < 8; cNumFine++) {
        map2[fineElemNum][cNumFine] = new us2Ptr[8];
        for(int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
          map2[fineElemNum][cNumFine][cNumCoarse] = new usPtr[18];
          for(int eTypeCoarse = 0; eTypeCoarse < 18; eTypeCoarse++) {
            map2[fineElemNum][cNumFine][cNumCoarse][eTypeCoarse] = new unsigned short[8];
            for(int vCtr = 0; vCtr < 8; vCtr++) {
              fscanf(infile,"%hu",
                  &(map2[fineElemNum][cNumFine][cNumCoarse][eTypeCoarse][vCtr]));
            }
          }
        }
      }
    }

    map3 = new us4Ptr[7];
    map4 = new us5Ptr[7];
    for(int fineElemNum = 0; fineElemNum < 7; fineElemNum++) { 
      map3[fineElemNum] = new us3Ptr[2];
      map4[fineElemNum] = new us4Ptr[2];
      for(int scalingCtr = 0; scalingCtr < 2; scalingCtr++) {
        map3[fineElemNum][scalingCtr] = new us2Ptr[8];
        for(int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
          map3[fineElemNum][scalingCtr][cNumCoarse] = new usPtr[18];
          for(int eTypeCoarse = 0; eTypeCoarse < 18; eTypeCoarse++) {
            map3[fineElemNum][scalingCtr][cNumCoarse][eTypeCoarse] =new unsigned short[8];
            for(int vCtr = 0; vCtr < 8; vCtr++) {
              fscanf(infile,"%hu",
                  &(map3[fineElemNum][scalingCtr][cNumCoarse][eTypeCoarse][vCtr]));
            }
          }
        }
        map4[fineElemNum][scalingCtr] = new us3Ptr[8]; 
        for(int cNumFine = 0; cNumFine < 8; cNumFine++) {
          map4[fineElemNum][scalingCtr][cNumFine] = new us2Ptr[8];
          for(int cNumCoarse = 0; cNumCoarse < 8; cNumCoarse++) {
            map4[fineElemNum][scalingCtr][cNumFine][cNumCoarse] = new usPtr[18];
            for(int eTypeCoarse = 0; eTypeCoarse < 18; eTypeCoarse++) {
              map4[fineElemNum][scalingCtr][cNumFine][cNumCoarse][eTypeCoarse] = new unsigned short[8];
              for(int vCtr = 0; vCtr < 8; vCtr++) {
                fscanf(infile,"%hu",
                    &(map4[fineElemNum][scalingCtr][cNumFine][cNumCoarse][eTypeCoarse][vCtr]));
              }
            }
          }
        }
      }
    }

    fclose(infile);
    return 1;
  }//end of function

  int destroyVtxMaps(unsigned short ****&map1, unsigned short *****&map2,
      unsigned short *****&map3, unsigned short ******&map4) {

    for(int i = 0; i < 8; i++) {
      for(int j = 0; j < 8; j++) {
        for(int k = 0; k < 18; k++) {
          delete [] map1[i][j][k];
          map1[i][j][k] = NULL;
        }
        delete [] map1[i][j];
        map1[i][j] = NULL;
      }
      delete [] map1[i];
      map1[i] = NULL;
    }
    delete [] map1;
    map1 = NULL;

    for(int i = 0; i < 8; i++) {
      for(int j = 0; j < 8; j++) {
        for(int k = 0; k < 8; k++) {
          for(int l = 0;l < 18; l++) {
            delete [] map2[i][j][k][l];
            map2[i][j][k][l] = NULL;
          }
          delete [] map2[i][j][k];
          map2[i][j][k] = NULL;
        }
        delete [] map2[i][j];
        map2[i][j] = NULL;
      }
      delete [] map2[i];
      map2[i] = NULL;
    }
    delete [] map2;
    map2 = NULL;

    for(int i = 0; i < 7; i++) {
      for(int j = 0; j < 2; j++) {
        for(int k = 0; k < 8; k++) {
          for(int l = 0; l < 18; l++) {
            delete [] map3[i][j][k][l];
            map3[i][j][k][l] = NULL;
          }
          delete [] map3[i][j][k];
          map3[i][j][k] = NULL;
        }
        delete [] map3[i][j];
        map3[i][j] = NULL;
      }
      delete [] map3[i];
      map3[i] = NULL;
    }
    delete [] map3;
    map3 = NULL;

    for(int i = 0; i < 7; i++) {
      for(int j = 0; j < 2; j++) {
        for(int k = 0; k < 8; k++) {
          for(int l = 0; l < 8; l++) {
            for(int m = 0; m < 18; m++) {
              delete [] map4[i][j][k][l][m];
              map4[i][j][k][l][m] = NULL;
            }
            delete [] map4[i][j][k][l];
            map4[i][j][k][l] = NULL;
          }
          delete [] map4[i][j][k];
          map4[i][j][k] = NULL;
        }
        delete [] map4[i][j];
        map4[i][j] = NULL;
      }
      delete [] map4[i];
      map4[i] = NULL;
    }
    delete [] map4;
    map4 = NULL;

    return 1;
  }//end of function

}//end namespace



