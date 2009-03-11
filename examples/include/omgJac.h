
/**
  @file omgJac.h
  @author Rahul Sampath, rahul.sampath@gmail.com
  */

#ifndef __OMG_JAC_H
#define __OMG_JAC_H

//Stuff for the case where the material properties on the coarser grids are
//constructed by averaging the material properties of the immediately finer
//grid 
struct Jac2MFreeData {
  std::vector<double>* matProp;
  bool isFinestLevel;
  Mat Jmat_private;
  Vec inTmp;
  Vec outTmp;
};

struct Jac3MFreeData { 
  ot::DA* daf;
  bool isFinestLevel;
  bool isCoarsestLevel;
  std::vector<double>* matProp;
  std::vector<double>* matPropFine;
  bool changedPartition; 
  Mat JmatThisLevel;
  Mat BmatThisLevel;
  Mat Jmat_private;
  Vec inTmp;
  Vec outTmp;
};

struct DirichletJacData {
  unsigned char* bdyArr;
  Mat Jmat_private;
  Vec inTmp;
  Vec outTmp;
};

void SetDirichletJacContexts(ot::DAMG* damg);
void DestroyDirichletJacContexts(ot::DAMG* damg);

void SetUserContexts(ot::DAMG* damg);
void SetUserContextsCoarsestToFinest(ot::DAMG* damg);

void SetCoarseToFineFromPts(ot::DAMG* damg,
    const std::vector<double>& pts,
    const std::vector<double> & lapJump);

void SetUserContextsFromPts(ot::DAMG* damg, const std::vector<double>& pts,
    const std::vector<double> & lapJump);

void DestroyUserContexts(ot::DAMG* damg);

PetscErrorCode CreateTmpDirichletLaplacian(ot::DAMG damg,Mat *B);
PetscErrorCode CreateDirichletLaplacian(ot::DAMG damg,Mat *B);
PetscErrorCode ComputeDirichletLaplacian(ot::DAMG damg,Mat J, Mat B);

PetscErrorCode TmpDirichletLaplacianMatMult(Mat, Vec, Vec);
PetscErrorCode DirichletLaplacianMatMult(Mat, Vec, Vec);
PetscErrorCode DirichletLaplacianShellMatMult(Mat, Vec, Vec);
PetscErrorCode DirichletLaplacianMatGetDiagonal(Mat, Vec);
PetscErrorCode DirichletLaplacianMatDestroy(Mat);

PetscErrorCode CreateTmpDirichletJacobian(ot::DAMG damg,Mat *B);
PetscErrorCode CreateDirichletJacobian(ot::DAMG damg,Mat *B);
PetscErrorCode ComputeDirichletJacobian(ot::DAMG damg,Mat J, Mat B);

PetscErrorCode TmpDirichletJacobianMatMult(Mat, Vec, Vec);
PetscErrorCode DirichletJacobianMatMult(Mat, Vec, Vec);
PetscErrorCode DirichletJacobianShellMatMult(Mat, Vec, Vec);
PetscErrorCode DirichletJacobianMatGetDiagonal(Mat, Vec);
PetscErrorCode DirichletJacobianMatDestroy(Mat);

PetscErrorCode CreateJacobian1(ot::DAMG damg,Mat *B);
PetscErrorCode ComputeJacobian1(ot::DAMG damg,Mat J, Mat B);
PetscErrorCode CreateAndComputeFullJacobian1(ot::DAMG damg,Mat *B);
PetscErrorCode CreateAndComputeFullActiveJacobian1(ot::DAMG damg,Mat *B);

PetscErrorCode CreateAndComputeMassMatrix(ot::DAMG damg, Mat* B);

PetscErrorCode CreateJacobian2(ot::DAMG damg,Mat *B);
PetscErrorCode ComputeJacobian2(ot::DAMG damg,Mat J, Mat B);

PetscErrorCode Jacobian2MatMult(Mat, Vec, Vec);
PetscErrorCode Jacobian2ShellMatMult(Mat, Vec, Vec);
PetscErrorCode Jacobian2MatGetDiagonal(Mat, Vec);
PetscErrorCode Jacobian2MatDestroy(Mat);

PetscErrorCode ComputeSol(ot::DAMG damg, Vec expectedSol); 
PetscErrorCode ComputeRandomRHS(ot::DAMG damg, Vec in);
PetscErrorCode ComputeConsistentRandomRHS(ot::DAMG damg, Vec in);
PetscErrorCode ComputeRHS0(ot::DAMG damg, Vec rhs);
PetscErrorCode ComputeRHS1(ot::DAMG damg, Vec rhs);
PetscErrorCode ComputeRHS2(ot::DAMG damg, Vec in);
PetscErrorCode ComputeRHS3(ot::DAMG damg, Vec in);
PetscErrorCode ComputeRHS4(ot::DAMG damg, Vec in);
PetscErrorCode ComputeRHS5(ot::DAMG damg, Vec in);
PetscErrorCode ComputeRHS6(ot::DAMG damg, Vec in);
PetscErrorCode ComputeRHS7(ot::DAMG damg, Vec in);
PetscErrorCode ComputeRHS8(ot::DAMG damg, Vec in);
PetscErrorCode ComputeRHS9(ot::DAMG damg, Vec in);
PetscErrorCode ComputeTestFBM_RHS(ot::DAMG damg, Vec in);
PetscErrorCode ComputeFBM2_RHS(ot::DAMG damg, Vec in);
PetscErrorCode ComputeFBM2_RHS_Part1(ot::DAMG damg, Vec in);
PetscErrorCode ComputeFBM2_RHS_Part2(ot::DAMG damg, Vec in);
PetscErrorCode SetSolutionFBM2(ot::DAMG damg, Vec in);
PetscErrorCode EnforceZeroFBM2(ot::DAMG damg, Vec in);
PetscErrorCode ComputeFBM_RHS(ot::DAMG damg, Vec in);
PetscErrorCode ComputeFBM_RHS_Part1(ot::DAMG damg, Vec in);
PetscErrorCode ComputeFBM_RHS_Part2(ot::DAMG damg, Vec in);
PetscErrorCode ComputeFBM_RHS_Part3(ot::DAMG damg, Vec in);
PetscErrorCode SetSolutionFBM(ot::DAMG damg, Vec in);
PetscErrorCode EnforceZeroFBM(ot::DAMG damg, Vec in);
PetscErrorCode SetSolution5(ot::DA* da, Vec in);

double ComputeError5(ot::DA* da, Vec in);
double TestError5(ot::DA* da);

PetscErrorCode CreateJacobian3(ot::DAMG damg, Mat *B);
PetscErrorCode ComputeJacobian3(ot::DAMG damg, Mat J, Mat B);

PetscErrorCode Jacobian3MatMult(Mat, Vec, Vec);
PetscErrorCode Jacobian3ShellMatMult(Mat, Vec, Vec);
PetscErrorCode Jacobian3MatGetDiagonal(Mat, Vec);
PetscErrorCode Jacobian3MatDestroy(Mat);

void getActiveStateAndActiveCommForKSP_Shell_DirichletJac(Mat mat,
    bool & activeState, MPI_Comm & activeComm);

void getActiveStateAndActiveCommForKSP_Shell_Jac1(Mat mat,
    bool & activeState, MPI_Comm & activeComm);

void getActiveStateAndActiveCommForKSP_Shell_Jac2or3(Mat mat,
    bool & activeState, MPI_Comm & activeComm);

void getPrivateMatricesForKSP_Shell_DirichletJac(Mat mat,
    Mat *AmatPrivate, Mat *PmatPrivate, MatStructure* pFlag);

void getPrivateMatricesForKSP_Shell_Jac1(Mat mat,
    Mat *AmatPrivate, Mat *PmatPrivate, MatStructure* pFlag);

void getPrivateMatricesForKSP_Shell_Jac2(Mat mat,
    Mat *AmatPrivate, Mat *PmatPrivate, MatStructure* pFlag);

void getPrivateMatricesForKSP_Shell_Jac3(Mat mat,
    Mat *AmatPrivate, Mat *PmatPrivate, MatStructure* pFlag);

#endif

