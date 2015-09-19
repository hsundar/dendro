

/**
 @file RmatType2StencilsAux.C
 @brief Load stencils for Restriction
 @author Rahul S. Sampath, rahul.sampath@gmail.com
  */

#include <cstdio>
#include <iostream>
#include <cassert>

namespace ot {

  int readRmatType2Stencil(double ****&lut) {
    typedef double*** double3Ptr;
    typedef double** double2Ptr;
    typedef double* doublePtr;
    FILE* infile;
    int res;
    char fname[100];
    sprintf(fname,"RmatType2Stencils.inp");
    infile = fopen(fname,"r");
    if(!infile) {
      std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
      assert(false);
    }
    lut = new double3Ptr[8];
    for(int j=0;j<8;j++) {
      lut[j] = new double2Ptr[18];
      for(int k=0;k<18;k++) {
        lut[j][k] = new doublePtr[8];
        for(int m=0;m<8;m++) {
          lut[j][k][m] = new double[8];
          for(int n=0;n<8;n++) {
            res = fscanf(infile,"%lf",&(lut[j][k][m][n]));
          }//end for n
        }//end for m
      }//end for k
    }//end for j

    fclose(infile);
    return 1;
  }//end of function

  int IreadRmatType2Stencil(double ****&lut, int rank) {
    typedef double*** double3Ptr;
    typedef double** double2Ptr;
    typedef double* doublePtr;
    FILE* infile;
    int res;
    char fname[250];
    sprintf(fname,"RmatType2Stencils_%d.inp", rank);
    infile = fopen(fname,"r");
    if(!infile) {
      std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
      assert(false);
    }
    lut = new double3Ptr[8];
    for(int j=0;j<8;j++) {
      lut[j] = new double2Ptr[18];
      for(int k=0;k<18;k++) {
        lut[j][k] = new doublePtr[8];
        for(int m=0;m<8;m++) {
          lut[j][k][m] = new double[8];
          for(int n=0;n<8;n++) {
            res = fscanf(infile,"%lf",&(lut[j][k][m][n]));
          }//end for n
        }//end for m
      }//end for k
    }//end for j

    fclose(infile);
    return 1;
  }//end of function

  int destroyRmatType2Stencil(double ****&lut) {
    for(int j=0;j<8;j++) {
      for(int k=0;k<18;k++) {
        for(int m=0;m<8;m++) {
          delete [] lut[j][k][m];
          lut[j][k][m] = NULL;
        }//end for m
        delete [] lut[j][k];
        lut[j][k] = NULL;
      }//end for k
      delete [] lut[j];
      lut[j] = NULL;
    }//end for j
    delete [] lut;
    lut = NULL;
    return 1;
  }//end of function

}//end namespace



