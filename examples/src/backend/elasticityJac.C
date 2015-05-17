
/**
  @file elasticityJac.C
  @author Rahul Sampath, rahul.sampath@gmail.com
  */

#include "petscmat.h"
#include "omg.h"
#include "oda.h"
#include "odaUtils.h"
#include "elasticityJac.h"

#ifdef PETSC_USE_LOG
extern int elasticityDiagEvent;
extern int elasticityMultEvent;
extern int elasticityFinestDiagEvent;
extern int elasticityFinestMultEvent;
#endif

extern double**** LaplacianType2Stencil; 
extern double**** GradDivType2Stencil; 

void getActiveStateAndActiveCommForKSP_Shell_Elas(Mat mat,
    bool & activeState, MPI_Comm & activeComm) {
  PetscBool isshell;
  PetscObjectTypeCompare((PetscObject)mat, MATSHELL, &isshell);
  assert(isshell);
  ot::DAMG damg;
  MatShellGetContext(mat, (void**)(&damg));
  ot::DA* da = damg->da;
  activeState = da->iAmActive();
  activeComm = da->getCommActive();
}

void getPrivateMatricesForKSP_Shell_Elas(Mat mat,
    Mat *AmatPrivate, Mat *PmatPrivate, MatStructure* pFlag) {
  PetscBool isshell;
  PetscObjectTypeCompare((PetscObject)mat, MATSHELL, &isshell);
  assert(isshell);
  ot::DAMG damg;
  MatShellGetContext(mat, (void**)(&damg));
  ElasticityData* data = (static_cast<ElasticityData*>(damg->user));
  *AmatPrivate = data->Jmat_private;
  *PmatPrivate = data->Jmat_private;
  *pFlag = DIFFERENT_NONZERO_PATTERN;
}

#define ELASTICITY_ELEM_DIAG_BLOCK {\
  unsigned int idx = da->curr();\
  unsigned int lev = da->getLevel(idx);\
  double h = hFac*(1u << (maxD - lev));\
  double fac = h/2.0;\
  unsigned int indices[8];\
  da->getNodeIndices(indices);\
  unsigned char childNum = da->getChildNumber();\
  unsigned char hnMask = da->getHangingNodeIndex(idx);\
  unsigned char elemType = 0;\
  GET_ETYPE_BLOCK(elemType,hnMask,childNum)\
  for(int k = 0; k < 8; k++) {\
    if(bdyArr[indices[k]]) {\
      /*Dirichlet Node*/\
      for(int dof = 0; dof < 3; dof++) {\
        diagArr[(3*indices[k])+dof] = 1.0;\
      } /*end dof*/\
    } else { \
      for(int dof = 0; dof < 3; dof++) {\
        diagArr[(3*indices[k])+dof] += (fac*(\
              (mu*LaplacianType2Stencil[childNum][elemType][k][k])\
              + ((mu+lambda)*GradDivType2Stencil[childNum][elemType][(3*k) + dof][(3*k) + dof])));\
      } /*end dof*/\
    }\
  } /*end k*/\
}

#define ELASTICITY_DIAG_BLOCK {\
  ot::DA* da = damg->da;\
  ElasticityData* data = (static_cast<ElasticityData*>(damg->user));\
  iC(VecZeroEntries(diag));\
  PetscScalar *diagArr = NULL;\
  unsigned char* bdyArr = data->bdyArr;\
  double mu = data->mu;\
  double lambda = data->lambda;\
  unsigned int maxD;\
  double hFac;\
  /*Nodal,Non-Ghosted,Write,3 dof*/\
  da->vecGetBuffer(diag,diagArr,false,false,false,3);\
  if(da->iAmActive()) {\
    maxD = (da->getMaxDepth());\
    hFac = 1.0/((double)(1u << (maxD-1)));\
    /*Loop through All Elements including ghosted*/\
    for(da->init<ot::DA_FLAGS::ALL>();\
        da->curr() < da->end<ot::DA_FLAGS::ALL>();\
        da->next<ot::DA_FLAGS::ALL>()) {\
 ELASTICITY_ELEM_DIAG_BLOCK \
    } /*end i*/\
  } /*end if active*/\
  da->vecRestoreBuffer(diag,diagArr,false,false,false,3);\
  /*2 IOP = 1 FLOP. Loop counters are included too.*/\
  PetscLogFlops(235*(da->getGhostedElementSize()));\
}

PetscErrorCode ElasticityMatGetDiagonal(Mat J, Vec diag) {
  PetscFunctionBegin;

  PetscLogEventBegin(elasticityDiagEvent,diag,0,0,0);

  ot::DAMG damg;
  iC(MatShellGetContext(J, (void**)(&damg)));

  bool isFinestLevel = (damg->nlevels == 1);

  if(isFinestLevel) {
    PetscLogEventBegin(elasticityFinestDiagEvent,diag,0,0,0);
  }

  ELASTICITY_DIAG_BLOCK  

    if(isFinestLevel) {
      PetscLogEventEnd(elasticityFinestDiagEvent,diag,0,0,0);
    }

  PetscLogEventEnd(elasticityDiagEvent,diag,0,0,0);

  PetscFunctionReturn(0);
}

#undef ELASTICITY_DIAG_BLOCK 
#undef ELASTICITY_ELEM_DIAG_BLOCK 

#define ELASTICITY_ELEM_MULT_BLOCK {\
  unsigned int idx = da->curr();\
  unsigned int lev = da->getLevel(idx);\
  double h = hFac*(1u << (maxD - lev));\
  double fac = h/2.0;\
  unsigned int indices[8];\
  da->getNodeIndices(indices);\
  unsigned char childNum = da->getChildNumber();\
  unsigned char hnMask = da->getHangingNodeIndex(idx);\
  unsigned char elemType = 0;\
  GET_ETYPE_BLOCK(elemType,hnMask,childNum)\
  for(int k = 0;k < 8;k++) {\
    if(bdyArr[indices[k]]) {\
      /*Dirichlet Node Row*/\
      for(int dof = 0; dof < 3; dof++) {\
        outArr[(3*indices[k]) + dof] =  inArr[(3*indices[k]) + dof];\
      }/*end for dof*/\
    } else {\
      for(int j=0;j<8;j++) {\
        /*Avoid Dirichlet Node Columns*/\
        if(!(bdyArr[indices[j]])) {\
          for(int dof = 0; dof < 3; dof++) {\
            outArr[(3*indices[k]) + dof] += (mu*fac*LaplacianType2Stencil[childNum][elemType][k][j]\
                *inArr[(3*indices[j]) + dof]);\
          }/*end for dof*/\
          for(int dofOut = 0; dofOut < 3; dofOut++) {\
            for(int dofIn = 0; dofIn < 3; dofIn++) {\
              outArr[(3*indices[k]) + dofOut] += ((mu+lambda)*fac*\
                  (GradDivType2Stencil[childNum][elemType][(3*k) + dofOut][(3*j) + dofIn])\
                  *inArr[(3*indices[j]) + dofIn]);\
            }/*end for dofIn*/\
          }/*end for dofOut*/\
        }\
      }/*end for j*/\
    }\
  }/*end for k*/\
}

#define ELASTICITY_MULT_BLOCK {\
  ot::DA* da = damg->da;\
  ElasticityData* data = (static_cast<ElasticityData*>(damg->user));\
  iC(VecZeroEntries(out));\
  unsigned int maxD;\
  double hFac;\
  if(da->iAmActive()) {\
    maxD = da->getMaxDepth();\
    hFac = 1.0/((double)(1u << (maxD-1)));\
  }\
  PetscScalar *outArr=NULL;\
  PetscScalar *inArr=NULL;\
  unsigned char* bdyArr = data->bdyArr;\
  double mu = data->mu;\
  double lambda = data->lambda;\
  /*Nodal,Non-Ghosted,Read,3 dof*/\
  da->vecGetBuffer(in,inArr,false,false,true,3);\
  /*Nodal,Non-Ghosted,Write,3 dof*/\
  da->vecGetBuffer(out,outArr,false,false,false,3);\
  if(da->iAmActive()) {\
    da->ReadFromGhostsBegin<PetscScalar>(inArr,3);\
    /*Loop through All Independent Elements*/\
    for(da->init<ot::DA_FLAGS::INDEPENDENT>();\
        da->curr() < da->end<ot::DA_FLAGS::INDEPENDENT>();\
        da->next<ot::DA_FLAGS::INDEPENDENT>() ) {\
 ELASTICITY_ELEM_MULT_BLOCK \
    } /*end independent*/\
    da->ReadFromGhostsEnd<PetscScalar>(inArr);\
    /*Loop through All Dependent Elements,*/\
    /*i.e. Elements which have atleast one*/\
    /*vertex owned by this processor and at least one*/\
    /*vertex not owned by this processor.*/\
    for(da->init<ot::DA_FLAGS::DEPENDENT>();\
        da->curr() < da->end<ot::DA_FLAGS::DEPENDENT>();\
        da->next<ot::DA_FLAGS::DEPENDENT>() ) {\
 ELASTICITY_ELEM_MULT_BLOCK \
    } /*end loop for dependent elems*/\
  } /*end if active*/\
  da->vecRestoreBuffer(in,inArr,false,false,true,3);\
  da->vecRestoreBuffer(out,outArr,false,false,false,3);\
  /*2 IOP = 1 FLOP. Loop counters are included too.*/\
  PetscLogFlops(6855*(da->getGhostedElementSize()));\
}

PetscErrorCode ElasticityMatMult(Mat J, Vec in, Vec out)
{
  PetscFunctionBegin;

  PetscLogEventBegin(elasticityMultEvent,in,out,0,0);

  ot::DAMG damg;
  iC(MatShellGetContext(J, (void**)(&damg)));

  bool isFinestLevel = (damg->nlevels == 1);

  if(isFinestLevel) {
    PetscLogEventBegin(elasticityFinestMultEvent,in,out,0,0);
  }

  ELASTICITY_MULT_BLOCK 

    if(isFinestLevel) {
      PetscLogEventEnd(elasticityFinestMultEvent,in,out,0,0);
    }

  PetscLogEventEnd(elasticityMultEvent,in,out,0,0);

  PetscFunctionReturn(0);
}

#undef ELASTICITY_ELEM_MULT_BLOCK 
#undef ELASTICITY_MULT_BLOCK 

#define BUILD_FULL_ELASTICITY_ELEM_ADD_BLOCK {\
  unsigned int idx = da->curr();\
  unsigned int lev = da->getLevel(idx);\
  double h = hFac*(1u << (maxD - lev));\
  double fac = h/2.0;\
  unsigned int indices[8];\
  da->getNodeIndices(indices);\
  unsigned char childNum = da->getChildNumber();\
  unsigned char hnMask = da->getHangingNodeIndex(idx);\
  unsigned char elemType = 0;\
  GET_ETYPE_BLOCK(elemType,hnMask,childNum)\
  for(int k = 0;k < 8;k++) {\
    /*Avoid Dirichlet Node Rows during ADD_VALUES loop.*/\
    /*Need a separate INSERT_VALUES loop for those*/\
    if(!(bdyArr[indices[k]])) {\
      for(int j=0;j<8;j++) {\
        if(!(bdyArr[indices[j]])) {\
          for(int dof = 0; dof < 3; dof++) {\
            ot::MatRecord currRec;\
            currRec.rowIdx = indices[k];\
            currRec.colIdx = indices[j];\
            currRec.rowDim = dof;\
            currRec.colDim = dof;\
            currRec.val = (mu*fac*\
                LaplacianType2Stencil[childNum][elemType][k][j]);\
            records.push_back(currRec);\
          } /*end for dof*/\
          for(int dofOut = 0; dofOut < 3; dofOut++) {\
            for(int dofIn = 0; dofIn < 3; dofIn++) {\
              ot::MatRecord currRec;\
              currRec.rowIdx = indices[k];\
              currRec.colIdx = indices[j];\
              currRec.rowDim = dofOut;\
              currRec.colDim = dofIn;\
              currRec.val = ((mu+lambda)*fac*\
                  GradDivType2Stencil[childNum][elemType][(3*k)+dofOut][(3*j)+dofIn]);\
              records.push_back(currRec);\
            } /*end for dofIn*/\
          } /*end for dofOut*/\
        }\
      } /*end for j*/\
    }\
  } /*end for k*/\
  if(records.size() > 1000) {\
    /*records will be cleared inside the function*/\
    da->setValuesInMatrix(B, records, 3, ADD_VALUES);\
  }\
}

#define BUILD_FULL_ELASTICITY_ELEM_INSERT_BLOCK {\
  unsigned int idx = da->curr();\
  unsigned int lev = da->getLevel(idx);\
  double h = hFac*(1u << (maxD - lev));\
  double fac = h/2.0;\
  unsigned int indices[8];\
  da->getNodeIndices(indices);\
  unsigned char childNum = da->getChildNumber();\
  unsigned char hnMask = da->getHangingNodeIndex(idx);\
  unsigned char elemType = 0;\
  GET_ETYPE_BLOCK(elemType,hnMask,childNum)\
  for(int k = 0;k < 8;k++) {\
    /*Insert values for Dirichlet Node Rows only.*/\
    if(bdyArr[indices[k]]) {\
      for(int dof = 0; dof < 3; dof++) {\
        ot::MatRecord currRec;\
        currRec.rowIdx = indices[k];\
        currRec.colIdx = indices[k];\
        currRec.rowDim = dof;\
        currRec.colDim = dof;\
        currRec.val = 1.0;\
        records.push_back(currRec);\
      } /*end for dof*/\
    }\
  } /*end for k*/\
  if(records.size() > 1000) {\
    /*records will be cleared inside the function*/\
    da->setValuesInMatrix(B, records, 3, INSERT_VALUES);\
  }\
}

#define BUILD_FULL_ELASTICITY_BLOCK {\
  ot::DA* da = damg->da;\
  MatZeroEntries(B);\
  std::vector<ot::MatRecord> records;\
  unsigned int maxD;\
  double hFac;\
  if(da->iAmActive()) {\
    maxD = da->getMaxDepth();\
    hFac = 1.0/((double)(1u << (maxD-1)));\
  }\
  unsigned char* bdyArr = data->bdyArr;\
  double mu = data->mu;\
  double lambda = data->lambda;\
  if(da->iAmActive()) {\
    for(da->init<ot::DA_FLAGS::WRITABLE>();\
        da->curr() < da->end<ot::DA_FLAGS::WRITABLE>();\
        da->next<ot::DA_FLAGS::WRITABLE>()) {\
 BUILD_FULL_ELASTICITY_ELEM_ADD_BLOCK \
    } /*end writable*/\
  } /*end if active*/\
  da->setValuesInMatrix(B, records, 3, ADD_VALUES);\
  MatAssemblyBegin(B, MAT_FLUSH_ASSEMBLY);\
  MatAssemblyEnd(B, MAT_FLUSH_ASSEMBLY);\
  if(da->iAmActive()) {\
    for(da->init<ot::DA_FLAGS::WRITABLE>();\
        da->curr() < da->end<ot::DA_FLAGS::WRITABLE>();\
        da->next<ot::DA_FLAGS::WRITABLE>()) {\
 BUILD_FULL_ELASTICITY_ELEM_INSERT_BLOCK \
    } /*end writable*/\
  } /*end if active*/\
  da->setValuesInMatrix(B, records, 3, INSERT_VALUES);\
  MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);\
  MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);\
  if(data->bdyArr) {\
    delete [] (data->bdyArr);\
    data->bdyArr = NULL;\
  }\
}

PetscErrorCode ComputeElasticityMat(ot::DAMG damg, Mat J, Mat B) {
  //For matShells nothing to be done here.
  PetscFunctionBegin;

  ElasticityData* data = (static_cast<ElasticityData*>(damg->user));

  PetscBool isshell;
  PetscObjectTypeCompare((PetscObject)B, MATSHELL, &isshell);

  assert(J == B);

  if(isshell) {
    if( data->Jmat_private == NULL ) {
      //inactive processors will return
      PetscFunctionReturn(0);
    } else {
      J = data->Jmat_private;
      B = data->Jmat_private;
    }
  }

  PetscObjectTypeCompare((PetscObject)B, MATSHELL, &isshell);

  if(isshell) {
    PetscFunctionReturn(0);
  }

  BUILD_FULL_ELASTICITY_BLOCK

    PetscFunctionReturn(0);
}

#undef BUILD_FULL_ELASTICITY_ELEM_ADD_BLOCK 
#undef BUILD_FULL_ELASTICITY_ELEM_INSERT_BLOCK 
#undef BUILD_FULL_ELASTICITY_BLOCK 

void SetElasticityContexts(ot::DAMG* damg) {
  int       nlevels = damg[0]->nlevels; //number of multigrid levels
  PetscReal muVal = 1.0;
  PetscReal lambdaVal = 1.0;
  PetscBool optFound;
  PetscOptionsGetReal("elasticity","-_mu", &muVal, &optFound);
  PetscOptionsGetReal("elasticity","-_lambda", &lambdaVal, &optFound);
  for(int i = 0; i < nlevels; i++) {
    ElasticityData* ctx= new ElasticityData;
    ctx->mu = muVal;
    ctx->lambda = lambdaVal;
    ctx->bdyArr = NULL;
    ctx->Jmat_private = NULL;
    ctx->inTmp = NULL;
    ctx->outTmp = NULL;

    std::vector<unsigned char> tmpBdyFlags;
    std::vector<unsigned char> tmpBdyFlagsAux;
    unsigned char* bdyArrAux = NULL;

    //This will create a nodal, non-ghosted, 1 dof array
    assignBoundaryFlags(damg[i]->da, tmpBdyFlags);
    ((damg[i])->da)->vecGetBuffer<unsigned char>(tmpBdyFlags, ctx->bdyArr, false, false, true, 1);
    if(damg[i]->da->iAmActive()) {
      ((damg[i])->da)->ReadFromGhostsBegin<unsigned char>(ctx->bdyArr, 1);
    }

    if(damg[i]->da_aux) {
      assignBoundaryFlags( damg[i]->da_aux, tmpBdyFlagsAux);
      ((damg[i])->da_aux)->vecGetBuffer<unsigned char>(tmpBdyFlagsAux, bdyArrAux,
        false, false, true, 1);
      if(damg[i]->da_aux->iAmActive()) {
        ((damg[i])->da_aux)->ReadFromGhostsBegin<unsigned char>(bdyArrAux, 1);
      }
    }

    tmpBdyFlags.clear();
    tmpBdyFlagsAux.clear();

    if(damg[i]->da->iAmActive()) {
      ((damg[i])->da)->ReadFromGhostsEnd<unsigned char>(ctx->bdyArr);
    }

    if((damg[i])->da_aux) {
      if(damg[i]->da_aux->iAmActive()) {
        ((damg[i])->da_aux)->ReadFromGhostsEnd<unsigned char>(bdyArrAux);
      }
    }

    for(int loopCtr = 0; loopCtr < 2; loopCtr++) {
      ot::DA* da = NULL;
      unsigned char* suppressedDOFptr = NULL;
      unsigned char* bdyArrPtr = NULL;
      if(loopCtr == 0) {
        da = damg[i]->da;
        suppressedDOFptr = damg[i]->suppressedDOF;
        bdyArrPtr = ctx->bdyArr;
      } else {
        da = damg[i]->da_aux;
        suppressedDOFptr = damg[i]->suppressedDOFaux;
        bdyArrPtr = bdyArrAux;
      }
      if(da) {
        if(da->iAmActive()) {
          for(da->init<ot::DA_FLAGS::ALL>(); 
              da->curr() < da->end<ot::DA_FLAGS::ALL>();
              da->next<ot::DA_FLAGS::ALL>()) {
            unsigned int indices[8];
            da->getNodeIndices(indices);
            for(unsigned int k = 0; k < 8; k++) {
              for(unsigned int d = 0; d < 3; d++) {
                suppressedDOFptr[(3*indices[k]) + d] = bdyArrPtr[indices[k]];
              }
            }
          }
        }
      }
    }

    if(bdyArrAux) {
      delete [] bdyArrAux;
      bdyArrAux = NULL;
    }

    (damg[i])->user = ctx;
  }//end for i

}//end fn.

void DestroyElasticityContexts(ot::DAMG* damg) {
  int       nlevels = damg[0]->nlevels; //number of multigrid levels
  for(int i = 0; i < nlevels; i++) {
    ElasticityData* ctx = (static_cast<ElasticityData*>(damg[i]->user));
    if(ctx->bdyArr) {
      delete [] (ctx->bdyArr);
      ctx->bdyArr = NULL;
    }
    if(ctx->Jmat_private) {
      MatDestroy(&(ctx->Jmat_private));
      ctx->Jmat_private = NULL;
    }
    if(ctx->inTmp) {
      VecDestroy(&(ctx->inTmp));
      ctx->inTmp = NULL;
    }
    if(ctx->outTmp) {
      VecDestroy(&(ctx->outTmp));
      ctx->outTmp = NULL;
    }
    delete ctx;
    ctx = NULL;
  }
}//end fn.

PetscErrorCode ElasticityShellMatMult(Mat J, Vec in, Vec out)
{
  PetscFunctionBegin;

  ot::DAMG damg;
  iC(MatShellGetContext(J, (void**)(&damg)));

  ElasticityData* ctx = (static_cast<ElasticityData*>(damg->user));

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

PetscErrorCode CreateElasticityMat(ot::DAMG damg, Mat *jac) {
  PetscFunctionBegin;
  PetscInt buildFullCoarseMat;
  PetscInt buildFullMatAll;
  int totalLevels;
  PetscBool flg;
  PetscOptionsGetInt(PETSC_NULL,"-buildFullCoarseMat",&buildFullCoarseMat,&flg);
  PetscOptionsGetInt(PETSC_NULL,"-buildFullMatAll",&buildFullMatAll,&flg);
  if(buildFullMatAll) {
    buildFullCoarseMat = 1;
  }
  totalLevels = damg->totalLevels;
  ot::DA* da = damg->da;
  int myRank;
  MPI_Comm_rank(da->getComm(), &myRank);
  if( totalLevels == damg->nlevels ) {
    //This is the coarsest.
    if( (!myRank) && (buildFullCoarseMat) ) {
      std::cout<<"Building Full elasticity Mat at the coarsest level."<<std::endl;
    }
    char matType[30];
    if(buildFullCoarseMat) {
      if(!(da->computedLocalToGlobal())) {
        da->computeLocalToGlobalMappings();
      }
      PetscBool typeFound;
      PetscOptionsGetString(PETSC_NULL,"-fullJacMatType",matType,30,&typeFound);
      if(!typeFound) {
        std::cout<<"I need a MatType for the full matrix!"<<std::endl;
        assert(false);
      }
    }
    bool requirePrivateMats = (da->getNpesActive() != da->getNpesAll());    
    if(requirePrivateMats) {
      unsigned int  m,n;
      m=n=(3*(da->getNodeSize()));
      ElasticityData* ctx = (static_cast<ElasticityData*>(damg->user));
      if(da->iAmActive()) {
        if(buildFullCoarseMat) {
          da->createActiveMatrix(ctx->Jmat_private, matType, 3);
        } else {
          MatCreateShell(da->getCommActive(), m ,n, PETSC_DETERMINE,
              PETSC_DETERMINE, damg, (&(ctx->Jmat_private)));

          MatShellSetOperation(ctx->Jmat_private, MATOP_MULT,
              (void (*)(void)) ElasticityMatMult);

          MatShellSetOperation(ctx->Jmat_private, MATOP_GET_DIAGONAL,
              (void (*)(void)) ElasticityMatGetDiagonal);

          MatShellSetOperation(ctx->Jmat_private, MATOP_DESTROY,
              (void (*)(void)) ElasticityMatDestroy);
        }
        MatGetVecs(ctx->Jmat_private, &(ctx->inTmp), &(ctx->outTmp));
      } else {
        ctx->Jmat_private = NULL;
      }
      //Need a MATShell wrapper anyway. But, the matvecs are not implemented for
      //this matrix. However, a matmult function is required for compute
      //residuals
      MatCreateShell(damg->comm, m ,n, PETSC_DETERMINE, PETSC_DETERMINE, damg, jac);
      MatShellSetOperation(*jac ,MATOP_DESTROY, (void (*)(void)) ElasticityMatDestroy);
      MatShellSetOperation(*jac, MATOP_MULT, (void (*)(void)) ElasticityShellMatMult);
    } else {
      if(buildFullCoarseMat) {
        da->createMatrix(*jac, matType, 3);
      } else {
        unsigned int  m,n;
        m=n=(3*(da->getNodeSize()));
        MatCreateShell(damg->comm, m ,n, PETSC_DETERMINE, PETSC_DETERMINE, damg, jac);
        MatShellSetOperation(*jac ,MATOP_MULT, (void (*)(void)) ElasticityMatMult);
        MatShellSetOperation(*jac ,MATOP_GET_DIAGONAL, (void (*)(void)) ElasticityMatGetDiagonal);
        MatShellSetOperation(*jac ,MATOP_DESTROY, (void (*)(void)) ElasticityMatDestroy);
      }
    }
    if( (!myRank) && (buildFullCoarseMat) ) {
      std::cout<<"Finished Building Full elasticity Mat at the coarsest level."<<std::endl;
    }
  }else {  
    //This is not the coarsest level. No need to bother with KSP_Shell
    if(buildFullMatAll) {
      if( !myRank ) {
        std::cout<<"Building Full elasticity Mat at level: "<<(damg->nlevels)<<std::endl;
      }
      if(!(da->computedLocalToGlobal())) {
        da->computeLocalToGlobalMappings();
      }
      char matType[30];
      PetscBool typeFound;
      PetscOptionsGetString(PETSC_NULL,"-fullJacMatType",matType,30,&typeFound);
      if(!typeFound) {
        std::cout<<"I need a MatType for the full matrix!"<<std::endl;
        assert(false);
      }
      da->createMatrix(*jac, matType, 3);
      if(!myRank) {
        std::cout<<"Finished Building Full elasticity Mat at level: "<<(damg->nlevels)<<std::endl;
      }
    } else {
      //Create a MATShell
      //The size this processor owns ( without ghosts).
      unsigned int  m,n;
      m=n=(3*(da->getNodeSize()));
      MatCreateShell(damg->comm, m ,n, PETSC_DETERMINE, PETSC_DETERMINE, damg, jac);
      MatShellSetOperation(*jac ,MATOP_MULT, (void (*)(void)) ElasticityMatMult);
      MatShellSetOperation(*jac ,MATOP_GET_DIAGONAL, (void (*)(void)) ElasticityMatGetDiagonal);
      MatShellSetOperation(*jac ,MATOP_DESTROY, (void (*)(void)) ElasticityMatDestroy);
    }
  }
  PetscFunctionReturn(0);
}//end fn.

PetscErrorCode ElasticityMatDestroy(Mat J) {
  PetscFunctionBegin;
  //Nothing to be done here. No new pointers were created during creation. 
  PetscFunctionReturn(0);
}

//Functions required in order to use BlockDiag PC

void computeInvBlockDiagEntriesForElasticityMat(Mat J, double **invBlockDiagEntries) {
  ot::DAMG damg;
  MatShellGetContext(J, (void**)(&damg));
  ElasticityData* data = (static_cast<ElasticityData*>(damg->user));
  ot::DA* da = damg->da;
  unsigned int dof = 3;
  unsigned int nodeSize = damg->da->getNodeSize();

  //Initialize
  for(int i = 0; i < (dof*nodeSize); i++ ) {
    for(int j = 0; j < dof; j++) {
      invBlockDiagEntries[i][j] = 0.0;
    }//end for j
  }//end for i

  unsigned char* bdyArr = data->bdyArr;
  double mu = data->mu;
  double lambda = data->lambda;
  unsigned int maxD;
  double hFac;

  std::vector<double> blockDiagVec;
  da->createVector<double>(blockDiagVec, false, false, 9);
  for(unsigned int i = 0; i < blockDiagVec.size(); i++) {
    blockDiagVec[i] = 0.0;
  }

  double *blockDiagArr;
  /*Nodal,Non-Ghosted,Write,9 dof*/
  da->vecGetBuffer<double>(blockDiagVec, blockDiagArr, false, false, false, 9);

  if(da->iAmActive()) {
    maxD = (da->getMaxDepth());
    hFac = 1.0/((double)(1u << (maxD-1)));
    /*Loop through All Elements including ghosted*/
    for(da->init<ot::DA_FLAGS::ALL>();
        da->curr() < da->end<ot::DA_FLAGS::ALL>();
        da->next<ot::DA_FLAGS::ALL>()) {
      unsigned int idx = da->curr();
      unsigned int lev = da->getLevel(idx);
      double h = hFac*(1u << (maxD - lev));
      double fac = h/2.0;
      unsigned int indices[8];
      da->getNodeIndices(indices);
      unsigned char childNum = da->getChildNumber();
      unsigned char hnMask = da->getHangingNodeIndex(idx);
      unsigned char elemType = 0;
      GET_ETYPE_BLOCK(elemType,hnMask,childNum)
        for(int k = 0; k < 8; k++) {
          if(bdyArr[indices[k]]) {
            /*Dirichlet Node*/
            for(int dof = 0; dof < 3; dof++) {
              blockDiagArr[(9*indices[k]) + (3*dof) + dof] = 1.0;
            } /*end dof*/
          } else { 
            for(int dof = 0; dof < 3; dof++) {
              blockDiagArr[(9*indices[k])+(3*dof) + dof] += (mu*fac*
                  LaplacianType2Stencil[childNum][elemType][k][k]);
            } /*end dof*/
            for(int dofOut = 0; dofOut < 3; dofOut++) {
              for(int dofIn = 0; dofIn < 3; dofIn++) {
                blockDiagArr[(9*indices[k]) + (3*dofOut) + dofIn] +=
                  ((mu+lambda)*fac*
                   GradDivType2Stencil[childNum][elemType][(3*k) + dofOut][(3*k) + dofIn]);
              }/*end dofIn*/
            } /*end dofOut*/
          }
        } /*end k*/
    } /*end i*/
  } /*end if active*/

  da->vecRestoreBuffer<double>(blockDiagVec,blockDiagArr,false,false,false,9);

  for(unsigned int i = 0; i < nodeSize; i++) {
    double a11 = blockDiagVec[(9*i)];
    double a12 = blockDiagVec[(9*i)+1];
    double a13 = blockDiagVec[(9*i)+2];
    double a21 = blockDiagVec[(9*i)+3];
    double a22 = blockDiagVec[(9*i)+4];
    double a23 = blockDiagVec[(9*i)+5];
    double a31 = blockDiagVec[(9*i)+6];
    double a32 = blockDiagVec[(9*i)+7];
    double a33 = blockDiagVec[(9*i)+8];

    double detA = ((a11*a22*a33)-(a11*a23*a32)-(a21*a12*a33)
        +(a21*a13*a32)+(a31*a12*a23)-(a31*a13*a22));

    invBlockDiagEntries[3*i][0] = (a22*a33-a23*a32)/detA;

    invBlockDiagEntries[3*i][1] = -(a12*a33-a13*a32)/detA;

    invBlockDiagEntries[3*i][2] = (a12*a23-a13*a22)/detA;

    invBlockDiagEntries[(3*i)+1][0] = -(a21*a33-a23*a31)/detA;

    invBlockDiagEntries[(3*i)+1][1] = (a11*a33-a13*a31)/detA;

    invBlockDiagEntries[(3*i)+1][2] = -(a11*a23-a13*a21)/detA;

    invBlockDiagEntries[(3*i)+2][0] = (a21*a32-a22*a31)/detA;

    invBlockDiagEntries[(3*i)+2][1] = -(a11*a32-a12*a31)/detA;

    invBlockDiagEntries[(3*i)+2][2] = (a11*a22-a12*a21)/detA;
  }//end for i

  blockDiagVec.clear();
}


void getDofAndNodeSizeForElasticityMat(Mat J, unsigned int & dof, unsigned int & nodeSize) {
  ot::DAMG damg;
  MatShellGetContext(J, (void**)(&damg));

  dof = 3;
  nodeSize = damg->da->getNodeSize();
}

