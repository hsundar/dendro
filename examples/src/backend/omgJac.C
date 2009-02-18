
/**
  @file omgJac.C
  @author Rahul Sampath, rahul.sampath@gmail.com
  */

#include "petsc.h"
#include "petscmat.h"
#include "petscvec.h"
#include "omg.h"
#include "oda.h"
#include "omgJac.h"
#include "odaJac.h"


#ifdef PETSC_USE_LOG
extern int Jac2MultEvent;
extern int Jac2DiagEvent;
extern int Jac2FinestMultEvent;
extern int Jac2FinestDiagEvent;

extern int Jac3MultEvent;
extern int Jac3DiagEvent;
extern int Jac3FinestMultEvent;
extern int Jac3FinestDiagEvent;
#endif


extern double***** LaplacianType1Stencil; 
extern double***** MassType1Stencil; 

extern double**** LaplacianType2Stencil; 
extern double**** MassType2Stencil; 

PetscErrorCode Jacobian2ShellMatMult(Mat J, Vec in, Vec out) {
  PetscFunctionBegin;

  ot::DAMG damg;

  MatShellGetContext(J, (void**)(&damg));

  Jac2MFreeData* ctx = static_cast<Jac2MFreeData*>(damg->user);

  if(damg->da->iAmActive()) {      
    PetscScalar* inArray;
    PetscScalar* outArray;

    VecGetArray(in, &inArray);
    VecGetArray(out, &outArray);

    VecPlaceArray(ctx->inTmp, inArray);
    VecPlaceArray(ctx->outTmp, outArray);

    MatMult(ctx->Jmat_private, ctx->inTmp, ctx->outTmp);

    VecResetArray(ctx->inTmp);
    VecResetArray(ctx->outTmp);

    VecRestoreArray(in, &inArray);
    VecRestoreArray(out, &outArray);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode Jacobian3ShellMatMult(Mat J, Vec in, Vec out) {
  PetscFunctionBegin;

  ot::DAMG damg;

  MatShellGetContext(J, (void**)(&damg));

  Jac3MFreeData* ctx = static_cast<Jac3MFreeData*>(damg->user);

  if(damg->da->iAmActive()) {      
    PetscScalar* inArray;
    PetscScalar* outArray;

    VecGetArray(in, &inArray);
    VecGetArray(out, &outArray);

    VecPlaceArray(ctx->inTmp, inArray);
    VecPlaceArray(ctx->outTmp, outArray);

    MatMult(ctx->Jmat_private, ctx->inTmp, ctx->outTmp);

    VecResetArray(ctx->inTmp);
    VecResetArray(ctx->outTmp);

    VecRestoreArray(in, &inArray);
    VecRestoreArray(out, &outArray);
  }

  PetscFunctionReturn(0);
}

void getActiveStateAndActiveCommForKSP_Shell_Jac1(Mat mat,
    bool & activeState, MPI_Comm & activeComm) {
  PetscTruth isshell;
  PetscTypeCompare((PetscObject)mat, MATSHELL, &isshell);
  assert(isshell);
  Jac1MFreeData *data;
  MatShellGetContext(mat, (void**)(&data));
  ot::DA* da = data->da;
  activeState = da->iAmActive();
  activeComm = da->getCommActive();
}

void getActiveStateAndActiveCommForKSP_Shell_Jac2or3(Mat mat,
    bool & activeState, MPI_Comm & activeComm) {
  PetscTruth isshell;
  PetscTypeCompare((PetscObject)mat, MATSHELL, &isshell);
  assert(isshell);
  ot::DAMG damg;
  MatShellGetContext(mat, (void**)(&damg));
  ot::DA* da = damg->da;
  activeState = da->iAmActive();
  activeComm = da->getCommActive();
}

void getPrivateMatricesForKSP_Shell_Jac1(Mat mat,
    Mat *AmatPrivate, Mat *PmatPrivate, MatStructure* pFlag) {
  PetscTruth isshell;
  PetscTypeCompare((PetscObject)mat, MATSHELL, &isshell);
  assert(isshell);
  Jac1MFreeData *data;
  MatShellGetContext(mat, (void**)(&data));
  *AmatPrivate = data->Jmat_private;
  *PmatPrivate = data->Jmat_private;
  *pFlag = DIFFERENT_NONZERO_PATTERN;
}

void getPrivateMatricesForKSP_Shell_Jac2(Mat mat,
    Mat *AmatPrivate, Mat *PmatPrivate, MatStructure* pFlag) {
  PetscTruth isshell;
  PetscTypeCompare((PetscObject)mat, MATSHELL, &isshell);
  assert(isshell);
  ot::DAMG damg;
  MatShellGetContext(mat, (void**)(&damg));
  Jac2MFreeData *data = (static_cast<Jac2MFreeData*>(damg->user));
  *AmatPrivate = data->Jmat_private;
  *PmatPrivate = data->Jmat_private;
  *pFlag = DIFFERENT_NONZERO_PATTERN;
}

void getPrivateMatricesForKSP_Shell_Jac3(Mat mat,
    Mat *AmatPrivate, Mat *PmatPrivate, MatStructure* pFlag) {
  PetscTruth isshell;
  PetscTypeCompare((PetscObject)mat, MATSHELL, &isshell);
  assert(isshell);
  ot::DAMG damg;
  MatShellGetContext(mat, (void**)(&damg));
  Jac3MFreeData *data = (static_cast<Jac3MFreeData*>(damg->user));
  *AmatPrivate = data->Jmat_private;
  *PmatPrivate = data->Jmat_private;
  *pFlag = DIFFERENT_NONZERO_PATTERN;
}

PetscErrorCode CreateAndComputeMassMatrix(ot::DAMG damg, Mat* jac) {
  unsigned int  m,n;
  PetscFunctionBegin;
  ot::DA* da = damg->da;
  //The size this processor owns ( without ghosts).
  m=n=da->getNodeSize();
  Jac1MFreeData* data = new Jac1MFreeData;
  data->da = da;
  data->inTmp = NULL;
  data->outTmp = NULL;
  data->Jmat_private = NULL;
  iC(MatCreateShell(da->getComm(), m ,n,PETSC_DETERMINE,PETSC_DETERMINE, (void*)(data), jac));
  iC(MatShellSetOperation(*jac ,MATOP_MULT, (void(*)(void)) MassMatMult));
  iC(MatShellSetOperation(*jac ,MATOP_GET_DIAGONAL, (void(*)(void)) MassMatGetDiagonal));
  iC(MatShellSetOperation(*jac ,MATOP_DESTROY, (void(*)(void)) MassMatDestroy));
  PetscFunctionReturn(0);
}

/* WRAPPERS FOR JACOBIAN TYPE-1 MATRICES */
PetscErrorCode CreateJacobian1(ot::DAMG damg,Mat *jac) {
  PetscFunctionBegin;
  PetscInt buildFullCoarseMat;
  PetscInt buildFullMatAll;
  int totalLevels;
  PetscTruth flg;
  PetscOptionsGetInt(PETSC_NULL, "-buildFullMatAll", &buildFullMatAll, &flg);
  PetscOptionsGetInt(PETSC_NULL, "-buildFullCoarseMat", &buildFullCoarseMat, &flg);
  if(buildFullMatAll) {
    buildFullCoarseMat = 1;
  }
  totalLevels = damg->totalLevels;
  ot::DA* da = damg->da;
  int myRank;
  MPI_Comm_rank(da->getComm(), &myRank);
  //The size this processor owns ( without ghosts).
  unsigned int  m,n;
  m=n=da->getNodeSize();
  if(totalLevels == damg->nlevels) {
    //This is the coarsest.
    if( (!myRank) && buildFullCoarseMat ) {
      std::cout<<"Building Full Coarse Mat."<<std::endl;
    }
    if(buildFullCoarseMat) {
      if(!(da->computedLocalToGlobal())) {
        da->computeLocalToGlobalMappings();
      }
    }
    bool requirePrivateMats = (da->getNpesActive() != da->getNpesAll());    
    if(requirePrivateMats) {
      Jac1MFreeData* data = new Jac1MFreeData;
      data->da = da;
      if(damg->nlevels == 1) {
        data->isFinestLevel = true;
      }else {
        data->isFinestLevel = false;
      }
      if(da->iAmActive()) {
        if(buildFullCoarseMat) {
          CreateAndComputeFullActiveJacobian1(damg, &(data->Jmat_private));
        } else {
          iC(MatCreateShell(damg->da->getCommActive(), m ,n,PETSC_DETERMINE,PETSC_DETERMINE,
                (void*)(data), &(data->Jmat_private)));

          iC(MatShellSetOperation( data->Jmat_private, MATOP_MULT,
                (void(*)(void)) Jacobian1MatMult));

          iC(MatShellSetOperation( data->Jmat_private, MATOP_GET_DIAGONAL,
                (void(*)(void)) Jacobian1MatGetDiagonal));

          iC(MatShellSetOperation(data->Jmat_private ,MATOP_DESTROY,
                (void(*)(void)) Jacobian1MatDestroy));
        }
        MatGetVecs(data->Jmat_private, &(data->inTmp), &(data->outTmp));
      } else {
        data->Jmat_private = NULL;
        data->inTmp = NULL;
        data->outTmp = NULL;
      }
      MatCreateShell(damg->comm, m ,n, PETSC_DETERMINE,PETSC_DETERMINE,
          (void*)(data), jac);
      MatShellSetOperation(*jac ,MATOP_DESTROY, (void(*)(void)) Jacobian1MatDestroy);
      MatShellSetOperation(*jac ,MATOP_MULT, (void(*)(void)) Jacobian1ShellMatMult);
    } else {
      if(buildFullCoarseMat) {
        CreateAndComputeFullJacobian1(damg, jac);
      } else {
        Jac1MFreeData* data = new Jac1MFreeData;
        data->da = da;
        if(damg->nlevels == 1) {
          data->isFinestLevel = true;
        }else {
          data->isFinestLevel = false;
        }
        data->Jmat_private = NULL;
        data->inTmp = NULL;
        data->outTmp = NULL;
        MatCreateShell(damg->comm, m ,n,PETSC_DETERMINE,PETSC_DETERMINE,
            (void*)(data), jac);

        MatShellSetOperation(*jac ,MATOP_MULT,
            (void(*)(void)) Jacobian1MatMult);

        MatShellSetOperation(*jac ,MATOP_GET_DIAGONAL,
            (void(*)(void)) Jacobian1MatGetDiagonal);

        MatShellSetOperation(*jac ,MATOP_DESTROY,
            (void(*)(void)) Jacobian1MatDestroy);
      }
    }
    if( (!myRank) && buildFullCoarseMat ) {
      std::cout<<"Finished Building Full Coarse Mat."<<std::endl;
    }
  } else {
    //This is some finer level. So no need for KSP_Shell
    if(buildFullMatAll) {
      if(!(da->computedLocalToGlobal())) {
        da->computeLocalToGlobalMappings();
      }
      if(!myRank) {
        std::cout<<"Building Full Mat for level: "<<(damg->nlevels)<<std::endl;
      }
      CreateAndComputeFullJacobian1(damg,jac);
      if(!myRank) {
        std::cout<<"Finished Building Full Mat for level: "<<(damg->nlevels)<<std::endl;
      }
    } else {
      Jac1MFreeData* data = new Jac1MFreeData;
      data->da = da;
      if(damg->nlevels == 1) {
        data->isFinestLevel = true;
      }else {
        data->isFinestLevel = false;
      }
      data->Jmat_private = NULL;
      data->inTmp = NULL;
      data->outTmp = NULL;
      MatCreateShell(damg->comm, m ,n,PETSC_DETERMINE,PETSC_DETERMINE,
          (void*)(data), jac);
      MatShellSetOperation(*jac ,MATOP_MULT, (void(*)(void)) Jacobian1MatMult);
      MatShellSetOperation(*jac ,MATOP_GET_DIAGONAL, (void(*)(void)) Jacobian1MatGetDiagonal);
      MatShellSetOperation(*jac ,MATOP_DESTROY, (void(*)(void)) Jacobian1MatDestroy);
    }    
  }

  PetscFunctionReturn(0);
}

PetscErrorCode CreateAndComputeFullJacobian1(ot::DAMG damg,Mat *mat)
{
  PetscFunctionBegin;
  CreateAndComputeFullJacobian1(damg->da, mat);
  PetscFunctionReturn(0);
}

PetscErrorCode CreateAndComputeFullActiveJacobian1(ot::DAMG damg,Mat *mat)
{
  PetscFunctionBegin;
  CreateAndComputeFullActiveJacobian1(damg->da, mat);
  PetscFunctionReturn(0);
}

PetscErrorCode ComputeJacobian1(ot::DAMG damg,Mat J, Mat B) {
  //damg not used here. might be used if user defined props exist for each level.
  //Since J and B are the same, only B is built.
  PetscFunctionBegin;
  //Nothing to do here
  PetscFunctionReturn(0);
}

PetscErrorCode CreateJacobian2(ot::DAMG damg, Mat *jac) {
  PetscFunctionBegin;
  int totalLevels;
  PetscTruth flg;
  PetscInt buildFullCoarseMat;
  PetscInt buildFullMatAll;
  PetscOptionsGetInt(PETSC_NULL,"-buildFullMatAll",&buildFullMatAll,&flg);
  PetscOptionsGetInt(PETSC_NULL,"-buildFullCoarseMat",&buildFullCoarseMat,&flg);
  if(buildFullMatAll) {
    buildFullCoarseMat = 1;
  }
  totalLevels = damg->totalLevels;
  ot::DA* da = damg->da;
  int myRank;
  MPI_Comm_rank(da->getComm(),&myRank);
  //The size this processor owns ( without ghosts).
  unsigned int  m,n;
  m=n=da->getNodeSize();
  if(totalLevels == damg->nlevels) {
    //This is the coarsest.
    if( (!myRank) && buildFullCoarseMat ) {
      std::cout<<"Building Full Coarse Mat."<<std::endl;
    }
    char matType[30];
    if(buildFullCoarseMat) {
      if(!(da->computedLocalToGlobal())) {
        da->computeLocalToGlobalMappings();
      }
      PetscTruth typeFound;
      PetscOptionsGetString(PETSC_NULL,"-fullJacMatType",matType,30,&typeFound);
      if(!typeFound) {
        std::cout<<"I need a MatType for the full Jacobian matrix!"<<std::endl;
        assert(false);
      } 
    }
    bool requirePrivateMats = (da->getNpesActive() != da->getNpesAll());    
    if(requirePrivateMats ) {
      Jac2MFreeData *data = (static_cast<Jac2MFreeData*>(damg->user));
      if(da->iAmActive()) {
        if(buildFullCoarseMat) {
          da->createActiveMatrix(data->Jmat_private, matType, 1);
        } else {
          MatCreateShell(da->getCommActive(), m, n, PETSC_DETERMINE,
              PETSC_DETERMINE, damg, &(data->Jmat_private));
          MatShellSetOperation(data->Jmat_private, MATOP_MULT,
              (void (*)(void)) Jacobian2MatMult);
          MatShellSetOperation(data->Jmat_private, MATOP_GET_DIAGONAL,
              (void (*)(void)) Jacobian2MatGetDiagonal);
          MatShellSetOperation(data->Jmat_private, MATOP_DESTROY,
              (void (*)(void)) Jacobian2MatDestroy);
        }
        MatGetVecs(data->Jmat_private, &(data->inTmp), &(data->outTmp));
      } else {
        data->Jmat_private = NULL;
        data->inTmp = NULL;
        data->outTmp = NULL;
      }
      MatCreateShell(damg->comm, m ,n, PETSC_DETERMINE, PETSC_DETERMINE, damg, jac);
      MatShellSetOperation(*jac ,MATOP_DESTROY, (void (*)(void)) Jacobian2MatDestroy);
      MatShellSetOperation(*jac ,MATOP_MULT, (void(*)(void)) Jacobian2ShellMatMult);
    } else {
      if(buildFullCoarseMat) {
        da->createMatrix(*jac, matType, 1);
      } else {
        MatCreateShell(damg->comm, m ,n, PETSC_DETERMINE, PETSC_DETERMINE, damg, jac);
        MatShellSetOperation(*jac ,MATOP_MULT, (void (*)(void)) Jacobian2MatMult);
        MatShellSetOperation(*jac ,MATOP_GET_DIAGONAL, (void (*)(void)) Jacobian2MatGetDiagonal);
        MatShellSetOperation(*jac ,MATOP_DESTROY, (void (*)(void)) Jacobian2MatDestroy);
      }
    }
    if((!myRank) && buildFullCoarseMat) {
      std::cout<<"Finished Building Full Coarse Mat."<<std::endl;
    }
  } else {
    //This is some finer level.
    if(buildFullMatAll) {
      if(!myRank) {
        std::cout<<"Building Full Mat for level: "<<(damg->nlevels)<<std::endl;
      }
      if(!(da->computedLocalToGlobal())) {
        da->computeLocalToGlobalMappings();
      }
      char matType[30];
      PetscTruth typeFound;
      PetscOptionsGetString(PETSC_NULL,"-fullJacMatType",matType,30,&typeFound);
      if(!typeFound) {
        std::cout<<"I need a MatType for the full Jacobian matrix!"<<std::endl;
        assert(false);
      } 
      da->createMatrix(*jac, matType, 1);
      if(!myRank) {
        std::cout<<"Finished Building Full Mat for level: "<<(damg->nlevels)<<std::endl;
      }
    } else {
      MatCreateShell(damg->comm, m ,n, PETSC_DETERMINE, PETSC_DETERMINE, damg, jac);
      MatShellSetOperation(*jac ,MATOP_MULT, (void (*)(void)) Jacobian2MatMult);
      MatShellSetOperation(*jac ,MATOP_GET_DIAGONAL, (void (*)(void)) Jacobian2MatGetDiagonal);
      MatShellSetOperation(*jac ,MATOP_DESTROY, (void (*)(void)) Jacobian2MatDestroy);
    }
  }

  PetscFunctionReturn(0);
}//end fn.

PetscErrorCode Jacobian2MatDestroy(Mat J) {
  PetscFunctionBegin;
  //Nothing to be done here. No new pointers were created during creation. 
  //The pointer were created in setUSerContext. So they will be destroyed in
  //DestroyUserContexts.
  PetscFunctionReturn(0);
}

#define JAC_TYPE2_ELEM_DIAG_BLOCK {\
  unsigned int idx = da->curr();\
  unsigned int lev = da->getLevel(idx);\
  double h = hFac*(1u << (maxD - lev));\
  double matP1 = matPropArr[2*idx];\
  double matP2 = matPropArr[2*idx+1];\
  double fac1 = matP1*h/2.0;\
  double fac2 = matP2*h*h*h/8.0;\
  unsigned int indices[8];\
  da->getNodeIndices(indices);\
  unsigned char childNum = da->getChildNumber();\
  unsigned char hnMask = da->getHangingNodeIndex(idx);\
  unsigned char elemType = 0;\
  GET_ETYPE_BLOCK(elemType,hnMask,childNum)\
  for(int k = 0; k < 8; k++) {\
    diagArr[indices[k]] +=  ((fac1*(LaplacianType2Stencil[childNum][elemType][k][k])) +\
        (fac2*(MassType2Stencil[childNum][elemType][k][k])));\
  } /*end k*/\
}

#define JAC_TYPE2_DIAG_BLOCK {\
  ot::DA* da = damg->da;\
  iC(VecZeroEntries(diag));\
  double *matPropArr;\
  /*Elem,Ghost,Read-only,1 dof*/\
  da->vecGetBuffer<double>(*(data->matProp),matPropArr,true,true,true,2);\
  PetscScalar *diagArr;\
  /*Nodal,Non-Ghosted,Write,1 dof*/\
  da->vecGetBuffer(diag,diagArr,false,false,false,1);\
  unsigned int maxD;\
  double hFac;\
  if(da->iAmActive()) {\
    maxD = (da->getMaxDepth());\
    hFac = 1.0/((double)(1u << (maxD-1)));\
    /*Loop through All Elements including ghosted*/\
    for(da->init<ot::DA_FLAGS::ALL>();\
        da->curr() < da->end<ot::DA_FLAGS::ALL>();\
        da->next<ot::DA_FLAGS::ALL>()) {\
      JAC_TYPE2_ELEM_DIAG_BLOCK\
    } /*end i*/\
  } /*end if active*/\
  da->vecRestoreBuffer<double>(*(data->matProp),matPropArr,true,true,true,2);\
  da->vecRestoreBuffer(diag,diagArr,false,false,false,1);\
  PetscLogFlops(44*(da->getGhostedElementSize()));\
}

PetscErrorCode Jacobian2MatGetDiagonal(Mat J, Vec diag) {
  PetscFunctionBegin;

  PetscLogEventBegin(Jac2DiagEvent,diag,0,0,0);

  ot::DAMG damg;
  iC(MatShellGetContext(J, (void**)(&damg)));

  Jac2MFreeData *data = (static_cast<Jac2MFreeData*>(damg->user));

  if(data->isFinestLevel) {
    PetscLogEventBegin(Jac2FinestDiagEvent,diag,0,0,0);
  }

  JAC_TYPE2_DIAG_BLOCK 

    if(data->isFinestLevel) {
      PetscLogEventEnd(Jac2FinestDiagEvent,diag,0,0,0);
    }

  PetscLogEventEnd(Jac2DiagEvent,diag,0,0,0);

  PetscFunctionReturn(0);
}

#define JAC_TYPE2_ELEM_MULT_BLOCK {\
  unsigned int idx = da->curr();\
  unsigned int lev = da->getLevel(idx);\
  double h = hFac*(1u << (maxD - lev));\
  double matP1 = matPropArr[2*idx];\
  double matP2 = matPropArr[2*idx+1];\
  double fac1 = matP1*h/2.0;\
  double fac2 = matP2*h*h*h/8.0;\
  unsigned int indices[8];\
  da->getNodeIndices(indices);\
  unsigned char childNum = da->getChildNumber();\
  unsigned char hnMask = da->getHangingNodeIndex(idx);\
  unsigned char elemType = 0;\
  GET_ETYPE_BLOCK(elemType,hnMask,childNum)\
  for(int k = 0;k < 8;k++) {\
    for(int j=0;j<8;j++) {\
      outArr[indices[k]] +=  (((fac1*(LaplacianType2Stencil[childNum][elemType][k][j])) +\
            (fac2*(MassType2Stencil[childNum][elemType][k][j])))*inArr[indices[j]]);\
    }/*end for j*/\
  }/*end for k*/\
}

#define JAC_TYPE2_MULT_BLOCK {\
  ot::DA* da = damg->da;\
  iC(VecZeroEntries(out));\
  unsigned int maxD;\
  double hFac;\
  if(da->iAmActive()) {\
    maxD = da->getMaxDepth();\
    hFac = 1.0/((double)(1u << (maxD-1)));\
  }\
  double *matPropArr = NULL;\
  /*Elem,Ghost,Read-only,2 dof*/\
  da->vecGetBuffer<double>(*(data->matProp),matPropArr,true,true,true,2);\
  PetscScalar *outArr=NULL;\
  PetscScalar *inArr=NULL;\
  /*Nodal,Non-Ghosted,Read,1 dof*/\
  da->vecGetBuffer(in,inArr,false,false,true,1);\
  if(da->iAmActive()) {\
    da->ReadFromGhostsBegin<PetscScalar>(inArr,1);\
  }\
  /*Nodal,Non-Ghosted,Write,1 dof*/\
  da->vecGetBuffer(out,outArr,false,false,false,1);\
  /*Loop through All Independent Elements*/\
  if(da->iAmActive()) {\
    for(da->init<ot::DA_FLAGS::INDEPENDENT>();\
        da->curr() < da->end<ot::DA_FLAGS::INDEPENDENT>();\
        da->next<ot::DA_FLAGS::INDEPENDENT>() ) {\
      JAC_TYPE2_ELEM_MULT_BLOCK\
    } /*end independent*/\
  } /*end if active*/\
  if(da->iAmActive()) {\
    da->ReadFromGhostsEnd<PetscScalar>(inArr);\
  }\
  /*Loop through All Dependent Elements,*/\
  /*i.e. Elements which have atleast one*/\
  /*vertex owned by this processor and at least one*/\
  /*vertex not owned by this processor.*/\
  if(da->iAmActive()) {\
    for(da->init<ot::DA_FLAGS::DEPENDENT>();\
        da->curr() < da->end<ot::DA_FLAGS::DEPENDENT>();\
        da->next<ot::DA_FLAGS::DEPENDENT>() ) {\
      JAC_TYPE2_ELEM_MULT_BLOCK\
    } /*end loop for dependent elems*/\
  } /*end if active*/\
  da->vecRestoreBuffer<double>(*(data->matProp),matPropArr,true,true,true,2);\
  da->vecRestoreBuffer(in,inArr,false,false,true,1);\
  da->vecRestoreBuffer(out,outArr,false,false,false,1);\
  PetscLogFlops(332*(da->getGhostedElementSize()));\
}

#define JAC_TYPE2_MULT_DEBUG_BLOCK {\
  ot::DA* da = damg->da;\
  iC(VecZeroEntries(out));\
  unsigned int maxD;\
  double hFac;\
  if(da->iAmActive()) {\
    maxD = da->getMaxDepth();\
    hFac = 1.0/((double)(1u << (maxD-1)));\
  }\
  double *matPropArr = NULL;\
  /*Elem,Ghost,Read-only,2 dof*/\
  da->vecGetBuffer<double>(*(data->matProp),matPropArr,true,true,true,2);\
  PetscScalar *outArr=NULL;\
  PetscScalar *inArr=NULL;\
  /*Nodal,Non-Ghosted,Read,1 dof*/\
  da->vecGetBuffer(in,inArr,false,false,true,1);\
  if(da->iAmActive()) {\
    da->ReadFromGhostsBegin<PetscScalar>(inArr,1);\
    da->ReadFromGhostsEnd<PetscScalar>(inArr);\
  }\
  /*Nodal,Non-Ghosted,Write,1 dof*/\
  da->vecGetBuffer(out,outArr,false,false,false,1);\
  if(da->iAmActive()) {\
    for(da->init<ot::DA_FLAGS::WRITABLE>();\
        da->curr() < da->end<ot::DA_FLAGS::WRITABLE>();\
        da->next<ot::DA_FLAGS::WRITABLE>() ) {\
      JAC_TYPE2_ELEM_MULT_BLOCK\
    } /*end loop */\
  } /*end if active*/\
  if(da->iAmActive()) {\
    da->WriteToGhostsBegin<PetscScalar>(outArr,1);\
    da->WriteToGhostsEnd<PetscScalar>(outArr,1);\
  }\
  da->vecRestoreBuffer<double>(*(data->matProp),matPropArr,true,true,true,2);\
  da->vecRestoreBuffer(in,inArr,false,false,true,1);\
  da->vecRestoreBuffer(out,outArr,false,false,false,1);\
  PetscLogFlops(332*(da->getGhostedElementSize()));\
}

PetscErrorCode Jacobian2MatMult(Mat J, Vec in, Vec out)
{
  PetscFunctionBegin;

  PetscLogEventBegin(Jac2MultEvent,in,out,0,0);

  ot::DAMG damg;
  iC(MatShellGetContext(J, (void**)(&damg)));

  Jac2MFreeData *data = (static_cast<Jac2MFreeData*>(damg->user));

  if(data->isFinestLevel) {
    PetscLogEventBegin(Jac2FinestMultEvent,in,out,0,0);
  }

  JAC_TYPE2_MULT_BLOCK 

    if(data->isFinestLevel) {
      PetscLogEventEnd(Jac2FinestMultEvent,in,out,0,0);
    }

  PetscLogEventEnd(Jac2MultEvent,in,out,0,0);

  PetscFunctionReturn(0);
}

#undef JAC_TYPE2_MULT_DEBUG_BLOCK 

#define BUILD_FULL_JAC_TYPE2_BLOCK(B) {\
  ot::DA* da = damg->da;\
  MatZeroEntries(B);\
  std::vector<ot::MatRecord> records;\
  unsigned int maxD;\
  double hFac;\
  if(da->iAmActive()) {\
    maxD = da->getMaxDepth();\
    hFac = 1.0/((double)(1u << (maxD-1)));\
  }\
  double *matPropArr = NULL;\
  /*Elem,Ghost,Read-only,2 dof*/\
  /*Note: All the flags are used for describing*/\
  /*the type of input vector, not*/\
  /*the type of the buffer.*/\
  da->vecGetBuffer<double>(*(data->matProp),matPropArr,true,true,true,2);\
  if(da->iAmActive()) {\
    for(da->init<ot::DA_FLAGS::WRITABLE>();\
        da->curr() < da->end<ot::DA_FLAGS::WRITABLE>();\
        da->next<ot::DA_FLAGS::WRITABLE>()) {\
      unsigned int idx = da->curr();\
      unsigned int lev = da->getLevel(idx);\
      double h = hFac*(1u << (maxD - lev));\
      double matP1 = matPropArr[2*idx];\
      double matP2 = matPropArr[2*idx+1];\
      double fac1 = matP1*h/2.0;\
      double fac2 = matP2*h*h*h/8.0;\
      unsigned int indices[8];\
      da->getNodeIndices(indices);\
      unsigned char childNum = da->getChildNumber();\
      unsigned char hnMask = da->getHangingNodeIndex(idx);\
      unsigned char elemType = 0;\
      GET_ETYPE_BLOCK(elemType,hnMask,childNum)\
      for(int k = 0;k < 8;k++) {\
        for(int j=0;j<8;j++) {\
          ot::MatRecord currRec;\
          currRec.rowIdx = indices[k];\
          currRec.colIdx = indices[j];\
          currRec.rowDim = 0;\
          currRec.colDim = 0;\
          currRec.val = ((fac1*(LaplacianType2Stencil[childNum][elemType][k][j])) +\
              (fac2*(MassType2Stencil[childNum][elemType][k][j])));\
          records.push_back(currRec);\
        } /*end for j*/\
      } /*end for k*/\
      if(records.size() > 1000) {\
        /*records will be cleared inside the function*/\
        da->setValuesInMatrix(B, records, 1, ADD_VALUES);\
      }\
    } /*end writable*/\
  } /*end if active*/\
  da->setValuesInMatrix(B, records, 1, ADD_VALUES);\
  MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);\
  da->vecRestoreBuffer<double>(*(data->matProp),matPropArr,true,true,true,2);\
  MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);\
}

PetscErrorCode ComputeJacobian2(ot::DAMG damg, Mat J, Mat B) {
  //For matShells nothing to be done here.
  PetscFunctionBegin;

  PetscTruth isshell;
  PetscTypeCompare((PetscObject)B, MATSHELL, &isshell);

  Jac2MFreeData *data = (static_cast<Jac2MFreeData*>(damg->user));

  assert(B == J);

  if(isshell) {
    if( data->Jmat_private == NULL ) {
      PetscFunctionReturn(0);
    } else {
      J = data->Jmat_private;
      B = data->Jmat_private;
    }
  }

  //B and J are the same.

  PetscTypeCompare((PetscObject)B, MATSHELL, &isshell);
  if(isshell) {
    //Nothing to be done here.
    PetscFunctionReturn(0);
  }

  BUILD_FULL_JAC_TYPE2_BLOCK(B) 

    PetscFunctionReturn(0);
}

/*Functions for Jacobian Type 3: Exact Fine grid integration using coarse grid
** shape functions  */
PetscErrorCode CreateJacobian3(ot::DAMG damg,Mat *jac) {
  PetscFunctionBegin;
  int totalLevels;
  PetscTruth flg;
  PetscInt buildFullMatAll;
  PetscInt buildFullCoarseMat;
  PetscOptionsGetInt(PETSC_NULL,"-buildFullMatAll",&buildFullMatAll,&flg);
  PetscOptionsGetInt(PETSC_NULL,"-buildFullCoarseMat",&buildFullCoarseMat,&flg);
  if(buildFullMatAll) {
    buildFullCoarseMat = 1;
  }
  totalLevels = damg->nlevels;
  ot::DA* da = damg->da;
  int myRank;
  MPI_Comm_rank(da->getComm(), &myRank);
  //The size this processor owns ( without ghosts).
  unsigned int  m,n;
  m=n=da->getNodeSize();
  if(totalLevels == damg->nlevels) {
    //This is the coarsest level
    if( (!myRank) && buildFullCoarseMat ) {
      std::cout<<"Building Full Mat for the coarsest level."<<std::endl;
    }
    char matType[30];
    if(buildFullCoarseMat) {
      if(!(da->computedLocalToGlobal())) {
        da->computeLocalToGlobalMappings();
      }
      PetscTruth typeFound;
      PetscOptionsGetString(PETSC_NULL,"-fullJacMatType",matType,30,&typeFound);
      if(!typeFound) {
        std::cout<<"I need a MatType for the full Jacobian matrix!"<<std::endl;
        assert(false);
      } 
    }
    bool requirePrivateMats = (da->getNpesActive() != da->getNpesAll());    
    if(requirePrivateMats) {
      Jac3MFreeData *data = (static_cast<Jac3MFreeData*>(damg->user));
      if(da->iAmActive()) {
        if(buildFullCoarseMat) {
          da->createActiveMatrix(data->Jmat_private, matType, 1);
        } else {
          MatCreateShell(damg->da->getCommActive(), m ,n, PETSC_DETERMINE,
              PETSC_DETERMINE, damg, &(data->Jmat_private) );

          MatShellSetOperation(data->Jmat_private ,MATOP_MULT,
              (void (*)(void)) Jacobian3MatMult);

          MatShellSetOperation(data->Jmat_private ,MATOP_GET_DIAGONAL,
              (void (*)(void)) Jacobian3MatGetDiagonal);

          MatShellSetOperation(data->Jmat_private ,MATOP_DESTROY,
              (void (*)(void)) Jacobian3MatDestroy);
        }
        MatGetVecs(data->Jmat_private, &(data->inTmp), &(data->outTmp));
      } else {
        data->Jmat_private = NULL;
        data->inTmp = NULL;
        data->outTmp = NULL;
      }
      MatCreateShell(damg->comm, m ,n, PETSC_DETERMINE, PETSC_DETERMINE, damg, jac);
      MatShellSetOperation(*jac ,MATOP_DESTROY, (void (*)(void)) Jacobian3MatDestroy);
      MatShellSetOperation(*jac ,MATOP_MULT, (void(*)(void)) Jacobian3ShellMatMult);
    } else {
      if(buildFullCoarseMat) {
        da->createMatrix(*jac, matType, 1);
      } else {
        MatCreateShell(damg->comm, m ,n, PETSC_DETERMINE, PETSC_DETERMINE, damg, jac);
        MatShellSetOperation(*jac ,MATOP_MULT, (void (*)(void)) Jacobian3MatMult);
        MatShellSetOperation(*jac ,MATOP_GET_DIAGONAL, (void (*)(void)) Jacobian3MatGetDiagonal);
        MatShellSetOperation(*jac ,MATOP_DESTROY, (void (*)(void)) Jacobian3MatDestroy);
      }
    }
    if( (!myRank) && buildFullCoarseMat ) {
      std::cout<<"Finished Building Full Mat for the coarsest level."<<std::endl;
    }
  } else {
    //This is some finer level
    if(buildFullMatAll) {
      if(!myRank) {
        std::cout<<"Building Full Mat for level: "<<(damg->nlevels)<<std::endl;
      }
      if(!(da->computedLocalToGlobal())) {
        da->computeLocalToGlobalMappings();
      }
      char matType[30];
      PetscTruth typeFound;
      PetscOptionsGetString(PETSC_NULL,"-fullJacMatType",matType,30,&typeFound);
      if(!typeFound) {
        std::cout<<"I need a MatType for the full Jacobian matrix!"<<std::endl;
        assert(false);
      } 
      da->createMatrix(*jac, matType, 1);
      if(!myRank) {
        std::cout<<"Finished Building Full Mat for level: "<<(damg->nlevels)<<std::endl;
      }
    } else {
      MatCreateShell(damg->comm, m ,n, PETSC_DETERMINE, PETSC_DETERMINE, damg, jac);
      MatShellSetOperation(*jac ,MATOP_MULT, (void (*)(void)) Jacobian3MatMult);
      MatShellSetOperation(*jac ,MATOP_GET_DIAGONAL, (void (*)(void)) Jacobian3MatGetDiagonal);
      MatShellSetOperation(*jac ,MATOP_DESTROY, (void (*)(void)) Jacobian3MatDestroy);
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode Jacobian3MatDestroy(Mat J) {
  PetscFunctionBegin;
  //Nothing to be done here. No new pointers were created during creation. 
  //The pointer were created in setUSerContext. So they will be destroyed in
  //DestroyUserContexts.
  PetscFunctionReturn(0);
}

#define JAC_TYPE3_DIAG_BLOCK {\
  /*The fine loop is always Writable, but the coarse loop*/\
  /*could be Independent or W_Dependent. Hence the fine counter must*/\
  /*be incremented properly to align with the coarse.*/\
  Point Cpt = da->getCurrentOffset();\
  while(daf->getCurrentOffset() != Cpt) {\
    daf->next<ot::DA_FLAGS::WRITABLE>();\
  }\
  unsigned int idxC= da->curr();\
  unsigned int lev = da->getLevel(idxC);\
  double h = hFac*(1u << (maxD - lev));\
  unsigned int indices[8];\
  da->getNodeIndices(indices);\
  unsigned char childNum = da->getChildNumber();\
  unsigned char hnMask = da->getHangingNodeIndex(idxC);\
  unsigned char elemType = 0;\
  GET_ETYPE_BLOCK(elemType,hnMask,childNum)\
  if(daf->getLevel(daf->curr()) == lev) {\
    /*The coarse and fine elements are the same,*/\
    type2Cnt++;\
    unsigned int idxF = daf->curr();\
    double matP1 = matPropArr[2*idxF];\
    double matP2 = matPropArr[2*idxF+1];\
    double fac1 = matP1*h/2.0;\
    double fac2 = matP2*h*h*h/8.0;\
    for(int k = 0; k < 8; k++) {\
      diagArr[indices[k]] +=  ((fac1*(LaplacianType2Stencil[childNum][elemType][k][k])) +\
          (fac2*(MassType2Stencil[childNum][elemType][k][k])));\
    } /*end k*/\
    daf->next<ot::DA_FLAGS::WRITABLE>();\
  }else {\
    type1Cnt++;\
    for(unsigned int cNumFine = 0; cNumFine < 8; cNumFine++) {\
      /*The coarse and fine elements are NOT the same. */\
      /*Loop over each of the 8 children of the coarse element.*/\
      /*These are the underlying fine elements.*/\
      unsigned int idxF = daf->curr();\
      double matP1 = matPropArr[2*idxF];\
      double matP2 = matPropArr[2*idxF+1];\
      double fac1 = matP1*h/2.0;\
      double fac2 = matP2*h*h*h/8.0;\
      for(int k = 0; k < 8; k++) {\
        diagArr[indices[k]] += \
        ((fac1*(LaplacianType1Stencil[childNum][elemType][cNumFine][k][k])) +\
         (fac2*(MassType1Stencil[childNum][elemType][cNumFine][k][k])));\
      } /*end k*/\
      daf->next<ot::DA_FLAGS::WRITABLE>();\
    } /*end loop over the 8 fine elements*/\
  }\
}

PetscErrorCode Jacobian3MatGetDiagonal(Mat J, Vec diag) {

  PetscFunctionBegin;

  PetscLogEventBegin(Jac3DiagEvent,0,0,0,0);

  ot::DAMG damg;
  iC(MatShellGetContext(J, (void**)(&damg)));

  Jac3MFreeData *data = (static_cast<Jac3MFreeData*>(damg->user));

  if(data->isFinestLevel) {
    PetscLogEventBegin(Jac3FinestDiagEvent,0,0,0,0);
  }

  if( (data->isFinestLevel) ||
      ( (data->JmatThisLevel != data->BmatThisLevel)
        && (J == data->BmatThisLevel) ) ) {
    JAC_TYPE2_DIAG_BLOCK 
  }else {
    //Use the fine grid material properties
    //Loop over the coarse and fine meshes simultaneously
    ot::DA* da = damg->da;
    ot::DA* daf = data->daf;
    assert(da->iAmActive() == daf->iAmActive());
    iC(VecZeroEntries(diag));
    double *matPropArr = NULL;
    if(data->changedPartition) {
      /*Elem,Non-Ghost,Read-only,1 dof*/
      daf->vecGetBuffer<double>(*(data->matPropFine),matPropArr,true,false,true,2);
    }else {
      /*Elem,Ghost,Read-only,1 dof*/
      daf->vecGetBuffer<double>(*(data->matPropFine),matPropArr,true,true,true,2);
    }
    PetscScalar *diagArr = NULL;
    /*Nodal,Non-Ghosted,Write,1 dof*/
    da->vecGetBuffer(diag,diagArr,false,false,false,1);
    unsigned int maxD;
    double hFac;
    unsigned int type1Cnt = 0; /*Coarse and Fine are not the same */
    unsigned int type2Cnt = 0; /*Coarse and Fine are the same */
    if(da->iAmActive()) {
      maxD = (da->getMaxDepth());
      hFac = 1.0/((double)(1u << (maxD-1)));
      for(da->init<ot::DA_FLAGS::W_DEPENDENT>(), daf->init<ot::DA_FLAGS::WRITABLE>();
          da->curr() < da->end<ot::DA_FLAGS::W_DEPENDENT>(); da->next<ot::DA_FLAGS::W_DEPENDENT>()) {
        JAC_TYPE3_DIAG_BLOCK 
      } /*end dependent loop*/
      da->WriteToGhostsBegin<PetscScalar>(diagArr, 1);
      for(da->init<ot::DA_FLAGS::INDEPENDENT>(), daf->init<ot::DA_FLAGS::WRITABLE>();
          da->curr() < da->end<ot::DA_FLAGS::INDEPENDENT>();  da->next<ot::DA_FLAGS::INDEPENDENT>()) {
        JAC_TYPE3_DIAG_BLOCK 
      } /*end Independent loop (overlapping with write to ghosts)*/
      da->WriteToGhostsEnd<PetscScalar>(diagArr, 1);
    } /*end if active*/

    if(data->changedPartition) {
      daf->vecRestoreBuffer<double>(*(data->matPropFine),matPropArr,true,false,true,2);
    }else {
      daf->vecRestoreBuffer<double>(*(data->matPropFine),matPropArr,true,true,true,2);
    }

    da->vecRestoreBuffer(diag,diagArr,false,false,false,1);

    PetscLogFlops( (type1Cnt*336) + (type2Cnt*44) );
  }

  if(data->isFinestLevel) {
    PetscLogEventEnd(Jac3FinestDiagEvent,0,0,0,0);
  }

  PetscLogEventEnd(Jac3DiagEvent,0,0,0,0);

  PetscFunctionReturn(0);
}

#undef JAC_TYPE3_DIAG_BLOCK 

#define JAC_TYPE3_MULT_BLOCK {\
  /*The fine loop is always Writable, but the coarse loop*/\
  /*could be Independent or W_Dependent. Hence the fine counter must*/\
  /*be incremented properly to align with the coarse.*/\
  Point Cpt = da->getCurrentOffset();\
  while(daf->getCurrentOffset() != Cpt) {\
    daf->next<ot::DA_FLAGS::WRITABLE>();\
  }\
  unsigned int idxC= da->curr();\
  unsigned int lev = da->getLevel(idxC);\
  double h = hFac*(1u << (maxD - lev));\
  unsigned int indices[8];\
  da->getNodeIndices(indices);\
  unsigned char childNum = da->getChildNumber();\
  unsigned char hnMask = da->getHangingNodeIndex(idxC);\
  unsigned char elemType = 0;\
  GET_ETYPE_BLOCK(elemType,hnMask,childNum)\
  if(daf->getLevel(daf->curr()) == lev) {\
    /*The coarse and fine elements are the same,*/\
    type2Cnt++;\
    unsigned int idxF = daf->curr();\
    double matP1 = matPropArr[2*idxF];\
    double matP2 = matPropArr[2*idxF+1];\
    double fac1 = matP1*h/2.0;\
    double fac2 = matP2*h*h*h/8.0;\
    for(int k = 0; k < 8; k++) {\
      for(int j = 0; j < 8; j++) {\
        outArr[indices[k]] += (((fac1*(LaplacianType2Stencil[childNum][elemType][k][j])) +\
              (fac2*(MassType2Stencil[childNum][elemType][k][j])))*inArr[indices[j]]);\
      } /*end j*/\
    } /*end k*/\
    daf->next<ot::DA_FLAGS::WRITABLE>();\
  }else {\
    double matP1[8];\
    double matP2[8];\
    for(unsigned int cNumFine = 0; cNumFine < 8; cNumFine++) {\
      /*The coarse and fine elements are NOT the same. */\
      /*Loop over each of the 8 children of the coarse element.*/\
      /*These are the underlying fine elements.*/\
      unsigned int idxF = daf->curr();\
      matP1[cNumFine] = matPropArr[2*idxF];\
      matP2[cNumFine] = matPropArr[2*idxF+1];\
      daf->next<ot::DA_FLAGS::WRITABLE>();\
    } /*end loop over the 8 fine elements*/\
    double maxP1 = matP1[0];\
    double minP1 = matP1[0];\
    double maxP2 = matP2[0];\
    double minP2 = matP2[0];\
    double totalP1 = 0.0;\
    double totalP2 = 0.0;\
    for(unsigned int cNumFine = 0; cNumFine < 8; cNumFine++) {\
      if(matP1[cNumFine] > maxP1) {\
        maxP1 = matP1[cNumFine];\
      }\
      if(matP2[cNumFine] > maxP2) {\
        maxP2 = matP2[cNumFine];\
      }\
      if(matP1[cNumFine] < minP1) {\
        minP1 = matP1[cNumFine];\
      }\
      if(matP2[cNumFine] < minP2) {\
        minP2 = matP2[cNumFine];\
      }\
      totalP1 += matP1[cNumFine];\
      totalP2 += matP2[cNumFine];\
    }\
    /*If the difference in the material properties is small*/\
    /*use the average value instead*/\
    if( ( (maxP1-minP1) <= (tolToCoarsenMatProp*maxP1) )\
        && ( (maxP2-minP2) <= (tolToCoarsenMatProp*maxP2) ) ) {\
      type2Cnt++;\
      double fac1 = totalP1*h/16.0;\
      double fac2 = totalP2*h*h*h/64.0;\
      for(int k = 0; k < 8; k++) {\
        for(int j = 0; j < 8; j++) {\
          outArr[indices[k]] += (((fac1*(LaplacianType2Stencil[childNum][elemType][k][j])) +\
                (fac2*(MassType2Stencil[childNum][elemType][k][j])))*inArr[indices[j]]);\
        } /*end j*/\
      } /*end k*/\
    }else {\
      type1Cnt++;\
      for(unsigned int cNumFine = 0; cNumFine < 8; cNumFine++) {\
        double fac1 = (matP1[cNumFine])*h/2.0;\
        double fac2 = (matP2[cNumFine])*h*h*h/8.0;\
        for(int k = 0; k < 8; k++) {\
          for(int j = 0; j < 8; j++) {\
            outArr[indices[k]] +=\
            (((fac1*(LaplacianType1Stencil[childNum][elemType][cNumFine][k][j])) +\
              (fac2*(MassType1Stencil[childNum][elemType][cNumFine][k][j])))*inArr[indices[j]]);\
          } /*end j*/\
        } /*end k*/\
      } /*end loop over the 8 fine elements*/\
    }\
  }\
}

PetscErrorCode Jacobian3MatMult(Mat J, Vec in, Vec out)
{
  PetscFunctionBegin;

  PetscLogEventBegin(Jac3MultEvent,in,out,0,0);

  ot::DAMG damg;
  iC(MatShellGetContext(J, (void**)(&damg)));

  Jac3MFreeData *data = (static_cast<Jac3MFreeData*>(damg->user));

  if(data->isFinestLevel) {
    PetscLogEventBegin(Jac3FinestMultEvent,in,out,0,0);
  }

  //The finest level and the pcMat for other levels use matProp directly
  //The Jmat for other levels use matPropFine
  //Only the finest and the coarsest grids have Jmat==Bmat
  //For the coarsest grid matPropFine must be used
  if( (data->isFinestLevel) ||
      ( (data->BmatThisLevel != data->JmatThisLevel) 
        && (J == data->BmatThisLevel) ) ) {
    JAC_TYPE2_MULT_BLOCK 
  }else {
    //Use the fine grid material properties
    //Loop over the coarse and fine meshes simultaneously
    PetscReal tolToCoarsenMatProp = 1.0e-12;
    PetscTruth optFound;
    PetscOptionsGetReal(0,"-tolToCoarsenMatProp",&tolToCoarsenMatProp,&optFound);
    assert(tolToCoarsenMatProp >= 0.0);

    ot::DA* da = damg->da;
    ot::DA* daf = data->daf;
    assert(da->iAmActive() == daf->iAmActive());
    double *matPropArr = NULL;
    if(data->changedPartition) {
      /*Elem,Non-Ghost,Read-only,1 dof*/
      daf->vecGetBuffer<double>(*(data->matPropFine),matPropArr,true,false,true,2);
    }else {
      /*Elem,Ghost,Read-only,1 dof*/
      daf->vecGetBuffer<double>(*(data->matPropFine),matPropArr,true,true,true,2);
    }
    PetscScalar *inArr = NULL;
    PetscScalar *outArr = NULL;
    /*Nodal,Non-Ghosted,ReadOnly,1 dof*/
    da->vecGetBuffer(in, inArr, false, false, true, 1);
    if(da->iAmActive()) {
      da->ReadFromGhostsBegin<PetscScalar>(inArr, 1);
    }
    iC(VecZeroEntries(out));
    /*Nodal,Non-Ghosted,Write,1 dof*/
    da->vecGetBuffer(out, outArr, false, false, false, 1);
    unsigned int maxD;
    double hFac;
    unsigned int type1Cnt = 0; /*Coarse and Fine are not the same */
    unsigned int type2Cnt = 0; /*Coarse and Fine are the same */
    if(da->iAmActive()) {
      maxD = (da->getMaxDepth());
      hFac = 1.0/((double)(1u << (maxD-1)));
      unsigned int loopCtr = 0;
      unsigned int numItersFirstLoop = static_cast<unsigned int>(0.3*static_cast<double>(da->getElementSize()));
      for(da->init<ot::DA_FLAGS::INDEPENDENT>(), daf->init<ot::DA_FLAGS::WRITABLE>();
          ( (daf->currWithInfo() == daf->currWithInfo()) && 
            (da->currWithInfo() < da->end<ot::DA_FLAGS::INDEPENDENT>()) && (loopCtr < numItersFirstLoop) );
          da->next<ot::DA_FLAGS::INDEPENDENT>(), loopCtr++) {
        JAC_TYPE3_MULT_BLOCK 
      } /*end Independent loop (overlapping with read from ghosts)*/
    } /* end if active */
    if(da->iAmActive()) {
      da->ReadFromGhostsEnd<PetscScalar>(inArr);
    }
    if(da->iAmActive()) {
      for(da->init<ot::DA_FLAGS::W_DEPENDENT>(), daf->init<ot::DA_FLAGS::WRITABLE>();
          da->curr() < da->end<ot::DA_FLAGS::W_DEPENDENT>(); da->next<ot::DA_FLAGS::W_DEPENDENT>()) {
        JAC_TYPE3_MULT_BLOCK 
      } /*end dependent loop*/
      da->WriteToGhostsBegin<PetscScalar>(outArr, 1);
      for(da->init<ot::DA_FLAGS::FROM_STORED>(), daf->init<ot::DA_FLAGS::FROM_STORED>();
          da->curr() < da->end<ot::DA_FLAGS::INDEPENDENT>(); da->next<ot::DA_FLAGS::INDEPENDENT>()) {
        JAC_TYPE3_MULT_BLOCK 
      } /*end Independent loop (overlapping with write to ghosts)*/
      da->WriteToGhostsEnd<PetscScalar>(outArr, 1);
    } /*end if active*/

    if(data->changedPartition) {
      daf->vecRestoreBuffer<double>(*(data->matPropFine),matPropArr,true,false,true,2);
    }else {
      daf->vecRestoreBuffer<double>(*(data->matPropFine),matPropArr,true,true,true,2);
    }

    da->vecRestoreBuffer(out, outArr, false, false, false, 1);
    da->vecRestoreBuffer(in, inArr, false, false, true, 1);

    PetscLogFlops( (type1Cnt*2643) + (type2Cnt*333) );
  }

  if(data->isFinestLevel) {
    PetscLogEventEnd(Jac3FinestMultEvent,in,out,0,0);
  }
  PetscLogEventEnd(Jac3MultEvent,in,out,0,0);

  PetscFunctionReturn(0);
}

#undef JAC_TYPE3_MULT_BLOCK 

#define BUILD_FULL_JAC_TYPE3_BLOCK(B) {\
  ot::DA* da = damg->da;\
  ot::DA* daf = data->daf;\
  MatZeroEntries(B);\
  std::vector<ot::MatRecord> records;\
  unsigned int maxD;\
  double hFac;\
  if(da->iAmActive()) {\
    maxD = da->getMaxDepth();\
    hFac = 1.0/((double)(1u << (maxD-1)));\
  }\
  double *matPropArr = NULL;\
  if(data->changedPartition) {\
    /*Elem,Non-Ghosted,Read-only,2 dof*/\
    daf->vecGetBuffer<double>(*(data->matPropFine),matPropArr,\
        true,false,true,2);\
  }else {\
    /*Elem,Ghosted,Read-only,2 dof*/\
    daf->vecGetBuffer<double>(*(data->matPropFine),matPropArr,\
        true,true,true,2);\
  }\
  if(da->iAmActive()) {\
    assert(daf->iAmActive());\
    for(da->init<ot::DA_FLAGS::WRITABLE>(),\
        daf->init<ot::DA_FLAGS::WRITABLE>();\
        da->curr() < da->end<ot::DA_FLAGS::WRITABLE>();\
        da->next<ot::DA_FLAGS::WRITABLE>()) {\
      Point Cpt = da->getCurrentOffset();\
      while(daf->getCurrentOffset() != Cpt) {\
        daf->next<ot::DA_FLAGS::WRITABLE>();\
      }\
      unsigned int idxC = da->curr();\
      unsigned int lev = da->getLevel(idxC);\
      double h = hFac*(1u << (maxD - lev));\
      unsigned int indices[8];\
      da->getNodeIndices(indices);\
      unsigned char childNum = da->getChildNumber();\
      unsigned char hnMask = da->getHangingNodeIndex(idxC);\
      unsigned char elemType = 0;\
      GET_ETYPE_BLOCK(elemType,hnMask,childNum)\
      if(daf->getLevel(daf->curr()) == lev) {\
        /*The coarse and fine elements are the same,*/\
        unsigned int idxF = daf->curr();\
        double matP1 = matPropArr[2*idxF];\
        double matP2 = matPropArr[2*idxF+1];\
        double fac1 = matP1*h/2.0;\
        double fac2 = matP2*h*h*h/8.0;\
        for(int k = 0; k < 8; k++) {\
          for(int j = 0; j < 8; j++) {\
            ot::MatRecord currRec;\
            currRec.rowIdx = indices[k];\
            currRec.colIdx = indices[j];\
            currRec.rowDim = 0;\
            currRec.colDim = 0;\
            currRec.val = ((fac1*\
                  (LaplacianType2Stencil[childNum][elemType][k][j])) +\
                (fac2*\
                 (MassType2Stencil[childNum][elemType][k][j])));\
            records.push_back(currRec);\
          } /*end for j*/\
        } /*end for k*/\
        daf->next<ot::DA_FLAGS::WRITABLE>();\
      }else {\
        for(unsigned char cNumFine = 0; cNumFine < 8; cNumFine++) {\
          /*The coarse and fine elements are NOT the same. */\
          /*Loop over each of the 8 children of the coarse element.*/\
          /*These are the underlying fine elements.*/\
          unsigned int idxF = daf->curr();\
          double matP1 = matPropArr[2*idxF];\
          double matP2 = matPropArr[2*idxF+1];\
          double fac1 = matP1*h/2.0;\
          double fac2 = matP2*h*h*h/8.0;\
          for(int k = 0; k < 8; k++) {\
            for(int j = 0; j < 8; j++) {\
              ot::MatRecord currRec;\
              currRec.rowIdx = indices[k];\
              currRec.colIdx = indices[j];\
              currRec.rowDim = 0;\
              currRec.colDim = 0;\
              currRec.val = ((fac1*(LaplacianType1Stencil\
                      [childNum][elemType][cNumFine][k][j]))+\
                  (fac2*(MassType1Stencil\
                         [childNum][elemType][cNumFine][k][j])));\
              records.push_back(currRec);\
            } /*end for j*/\
          } /*end for k*/\
          daf->next<ot::DA_FLAGS::WRITABLE>();\
        } /*end loop over the 8 fine elements*/\
      }\
      if(records.size() > 1000) {\
        /*records will be cleared inside the function*/\
        da->setValuesInMatrix(B, records, 1, ADD_VALUES);\
      }\
    } /*end writable*/\
  } else {\
    assert(!(daf->iAmActive()));\
  } /*end if active*/\
  da->setValuesInMatrix(B, records, 1, ADD_VALUES);\
  MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);\
  if(data->changedPartition) {\
    daf->vecRestoreBuffer<double>(*(data->matPropFine),matPropArr,\
        true,false,true,2);\
  }else {\
    daf->vecRestoreBuffer<double>(*(data->matPropFine),matPropArr,\
        true,true,true,2);\
  }\
  MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);\
}

PetscErrorCode ComputeJacobian3(ot::DAMG damg, Mat J, Mat B) {
  PetscFunctionBegin;

  Jac3MFreeData *data = (static_cast<Jac3MFreeData*>(damg->user));

  PetscTruth isShellB, isShellJ;
  PetscTypeCompare((PetscObject)B, MATSHELL, &isShellB);
  PetscTypeCompare((PetscObject)J, MATSHELL, &isShellJ);

  assert(isShellB == isShellJ);

  //For matShells nothing to be done here.
  if(isShellB) {
    if( data->Jmat_private == NULL ) {
      //This will be used to determine the type of matrix in the MatVecs
      data->JmatThisLevel = J;
      data->BmatThisLevel = B;

      PetscFunctionReturn(0);
    } else {
      //This must be the coarsest level: J and B will be the same
      J = data->Jmat_private;
      B = data->Jmat_private;
    }
  }

  //This will be used to determine the type of matrix in the MatVecs
  data->JmatThisLevel = J;
  data->BmatThisLevel = B;

  PetscTypeCompare((PetscObject)B, MATSHELL, &isShellB);
  PetscTypeCompare((PetscObject)J, MATSHELL, &isShellJ);

  if(J != B) {
    //Build both B ond J
    //Use matProp for J
    if(!isShellJ) {
      BUILD_FULL_JAC_TYPE2_BLOCK(J) 
    }
    //Use matPropFine for B
    //Loop over the coarse and fine meshes simultaneously
    if(!isShellB) {
      BUILD_FULL_JAC_TYPE3_BLOCK(B) 
    }
  } else {
    //Build B only
    if(!isShellB) {
      if(data->isFinestLevel) {
        //Use matProp
        BUILD_FULL_JAC_TYPE2_BLOCK(B) 
      } else {
        //This must be the coarsest grid
        //Use matPropFine
        //Loop over the coarse and fine meshes simultaneously
        BUILD_FULL_JAC_TYPE3_BLOCK(B) 
      }//end if finest
    }
  }//end if J == B

  PetscFunctionReturn(0);
}


#undef BUILD_FULL_JAC_TYPE2_BLOCK 
#undef BUILD_FULL_JAC_TYPE3_BLOCK 
#undef JAC_TYPE2_MULT_BLOCK 
#undef JAC_TYPE2_DIAG_BLOCK 
#undef JAC_TYPE2_ELEM_DIAG_BLOCK 
#undef JAC_TYPE2_ELEM_MULT_BLOCK 



