
/**
  @file elasticityJac.h
  @author Rahul Sampath, rahul.sampath@gmail.com
  */

#ifndef __ELASTICITY_JAC_H
#define __ELASTICITY_JAC_H

//Stuff for the case where the material properties on the coarser grids are
//constructed by averaging the material properties of the immediately finer
//grid 
struct ElasticityData {
  unsigned char* bdyArr;
  PetscReal mu;
  PetscReal lambda;
  Mat Jmat_private;
  Vec inTmp;
  Vec outTmp;
};

void SetElasticityContexts(ot::DAMG* damg);
void DestroyElasticityContexts(ot::DAMG* damg);

PetscErrorCode CreateElasticityMat(ot::DAMG damg,Mat *B);
PetscErrorCode ComputeElasticityMat(ot::DAMG damg,Mat J, Mat B);

PetscErrorCode ElasticityMatMult(Mat, Vec, Vec);
PetscErrorCode ElasticityMatGetDiagonal(Mat, Vec);
PetscErrorCode ElasticityShellMatMult(Mat, Vec, Vec);
PetscErrorCode ElasticityMatDestroy(Mat);

PetscErrorCode ComputeElasticityRHS(ot::DAMG damg,Vec rhs);

//Functions required for PC_BlockDiag.
void getDofAndNodeSizeForElasticityMat(Mat J, unsigned int & dof, unsigned int & nodeSize);

void computeInvBlockDiagEntriesForElasticityMat(Mat J, double **invBlockDiagEntries);

//Function for using KSP_Shell (will be used at the coarsest grid if not all
//processors are active on the coarsest grid)
void getPrivateMatricesForKSP_Shell_Elas(Mat mat,
    Mat *AmatPrivate, Mat *PmatPrivate, MatStructure* pFlag);

#endif

