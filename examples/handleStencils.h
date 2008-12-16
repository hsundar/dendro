
/**
  @file handleStencils.h
  @author Rahul Sampath, rahul.sampath@gmail.com
  */

#ifndef __HANDLE_STENCILS_H
#define __HANDLE_STENCILS_H

int createLmatType2(double ****& Lmat2);
int createMmatType2(double ****& Mmat2);
int createGDmatType2(double ****& GDmat2);

int createLmatType2_Type1(double ****& Lmat2);
int createMmatType2_Type1(double ****& Mmat2);
int createGDmatType2_Type1(double ****& GDmat2);

int createLmatType2_Type2(double ****& Lmat2);
int createMmatType2_Type2(double ****& Mmat2);
int createGDmatType2_Type2(double ****& GDmat2);

int createLmatType2_Type3(double ****& Lmat2);
int createMmatType2_Type3(double ****& Mmat2);
int createGDmatType2_Type3(double ****& GDmat2);

int createLmatType1(double *****& Lmat1);
int createMmatType1(double *****& Mmat1);

int createLmatType1_Type1(double *****& Lmat1);
int createMmatType1_Type1(double *****& Mmat1);

int createLmatType1_Type2(double *****& Lmat1);
int createMmatType1_Type2(double *****& Mmat1);

int createLmatType1_Type3(double *****& Lmat1);
int createMmatType1_Type3(double *****& Mmat1);

int createShFnMat(double******& shFnMat);
int createRHSType2(double***& RHS_data);

int createShFnMat_Type1(double******& shFnMat);
int createRHSType2_Type1(double***& RHS_data);

int createShFnMat_Type2(double******& shFnMat);
int createRHSType2_Type2(double***& RHS_data);

int createShFnMat_Type3(double******& shFnMat);
int createRHSType2_Type3(double***& RHS_data);

int destroyLmatType1(double *****& Lmat1);
int destroyMmatType1(double *****& Mmat1);

int destroyShFnMat(double******& shFnMat);
int destroyRHSType2(double***& RHS_data);

int destroyLmatType2(double ****& Lmat2);
int destroyMmatType2(double ****& Mmat2);
int destroyGDmatType2(double ****& GDmat2);

#endif

