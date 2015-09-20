
/**
  @file RmatType1StencilsAux.C
  @brief Load Stencils for Restriction
  @author Rahul S. Sampath,rahul.sampath@gmail.com
  */

#include <cstdio>
#include <iostream>
#include <cassert>

namespace ot {

  int readRmatType1Stencil(double *****&lut) {
    typedef double**** double4Ptr;
    typedef double*** double3Ptr;
    typedef double** double2Ptr;
    typedef double* doublePtr;
    FILE* infile;
    int res;
    char fname[100];
    sprintf(fname,"RmatType1Stencils.inp");
    infile = fopen(fname,"r");
    if(!infile) {
      std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
      assert(false);
    }
    lut = new double4Ptr[8];
    for(int i=0;i<8;i++) {
      lut[i] = new double3Ptr[8];
      for(int j=0;j<8;j++) {
        lut[i][j] = new double2Ptr[18];
        for(int k=0;k<18;k++) {
          lut[i][j][k] = new doublePtr[8];
          for(int l=0;l<8;l++) {
            lut[i][j][k][l] = new double[8];
            for(int m=0;m<8;m++) {
              res = fscanf(infile,"%lf",&(lut[i][j][k][l][m]));
            }//end for m
          }//end for l
        }//end for k
      }//end for j
    }//end for i

    fclose(infile);
    return 1;
  }//end of function

  int IreadRmatType1Stencil(double *****&lut, int rank) {
    typedef double**** double4Ptr;
    typedef double*** double3Ptr;
    typedef double** double2Ptr;
    typedef double* doublePtr;
    FILE* infile;
    int res;
    char fname[250];
    sprintf(fname,"RmatType1Stencils_%d.inp",rank);
    infile = fopen(fname,"r");
    if(!infile) {
      std::cout<<"The file "<<fname<<" is not good for reading."<<std::endl;
      assert(false);
    }
    lut = new double4Ptr[8];
    for(int i=0;i<8;i++) {
      lut[i] = new double3Ptr[8];
      for(int j=0;j<8;j++) {
        lut[i][j] = new double2Ptr[18];
        for(int k=0;k<18;k++) {
          lut[i][j][k] = new doublePtr[8];
          for(int l=0;l<8;l++) {
            lut[i][j][k][l] = new double[8];
            for(int m=0;m<8;m++) {
              res = fscanf(infile,"%lf",&(lut[i][j][k][l][m]));
            }//end for m
          }//end for l
        }//end for k
      }//end for j
    }//end for i

    fclose(infile);
    return 1;
  }//end of function


  int destroyRmatType1Stencil(double *****&lut) {
    for(int i=0;i<8;i++) {
      for(int j=0;j<8;j++) {
        for(int k=0;k<18;k++) {
          for(int l=0;l<8;l++) {
            delete [] lut[i][j][k][l];
            lut[i][j][k][l] = NULL;
          }//end for l
          delete [] lut[i][j][k];
          lut[i][j][k] = NULL;
        }//end for k
        delete [] lut[i][j];
        lut[i][j] = NULL;
      }//end for j
      delete [] lut[i]; 
      lut[i] = NULL;
    }//end for i
    delete [] lut;
    lut = NULL;
    return 1;
  }//end of function

}//end namespace


