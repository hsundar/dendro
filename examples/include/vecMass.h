
/**
  @file vecMass.h
  @author Rahul Sampath, rahul.sampath@gmail.com
  */

#ifndef __VEC_MASS_H
#define __VEC_MASS_H

PetscErrorCode CreateConstVecMass(ot::DAMG damg, Mat *B);
PetscErrorCode ComputeConstVecMass(ot::DAMG damg, Mat J, Mat B);

PetscErrorCode ConstVecMassMatMult(Mat, Vec, Vec);
PetscErrorCode ConstVecMassMatGetDiagonal(Mat, Vec);
PetscErrorCode ConstVecMassMatDestroy(Mat);

#endif

