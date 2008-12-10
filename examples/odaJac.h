
/**
  @file odaJac.h
  @author Rahul Sampath, rahul.sampath@gmail.com
  */

#ifndef __ODA_JAC_H
#define __ODA_JAC_H

#include "oda.h"
#include "petscmat.h"
#include "petscvec.h"

extern int _internal_ierr;

#ifndef iC
#define iC(fun) {CHKERRQ(fun);}
#endif

struct Jac1MFreeData {
  ot::DA *da;
  bool isFinestLevel;
  Mat Jmat_private;
  Vec inTmp;
  Vec outTmp;
};

PetscErrorCode CreateJacobian1(ot::DA*,Mat *);
PetscErrorCode CreateActiveJacobian1(ot::DA*,Mat *);
PetscErrorCode ComputeJacobian1(ot::DA*,Mat);
PetscErrorCode CreateAndComputeFullJacobian1(ot::DA*,Mat *);
PetscErrorCode CreateAndComputeFullActiveJacobian1(ot::DA*,Mat *);
PetscErrorCode CreateAndComputeMassMatrix(ot::DA*,Mat *);

PetscErrorCode Jacobian1MatMult(Mat, Vec, Vec);
PetscErrorCode Jacobian1ShellMatMult(Mat, Vec, Vec);

PetscErrorCode MassMatMult(Mat, Vec, Vec);
PetscErrorCode Jacobian1MatGetDiagonal(Mat, Vec);
PetscErrorCode MassMatGetDiagonal(Mat, Vec);
PetscErrorCode Jacobian1MatDestroy(Mat);
PetscErrorCode MassMatDestroy(Mat);

#endif

