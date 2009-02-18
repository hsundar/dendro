
/**
  @file vecMass.C
  @author Rahul Sampath, rahul.sampath@gmail.com
  */

#include "petsc.h"
#include "omg.h"
#include "oda.h"
#include "vecMass.h"
#include "elasticityJac.h"

#ifdef PETSC_USE_LOG
extern int vecMassDiagEvent;
extern int vecMassMultEvent;
extern int vecMassFinestDiagEvent;
extern int vecMassFinestMultEvent;
#endif

extern double**** MassType2Stencil; 

PetscErrorCode CreateConstVecMass(ot::DAMG damg, Mat *jac) {
  PetscFunctionBegin;
  PetscInt buildFullCoarseMat;
  int totalLevels;
  PetscTruth flg;
  PetscOptionsGetInt(PETSC_NULL,"-buildFullCoarseMat",&buildFullCoarseMat,&flg);
  totalLevels = damg->nlevels;
  ot::DA* da = damg->da;
  if((totalLevels == damg->nlevels) && buildFullCoarseMat) {
    //This is the coarsest.
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
    if(!myRank) {
      std::cout<<"Building Full vecMass Mat."<<std::endl;
    }
    if(!(da->computedLocalToGlobal())) {
      da->computeLocalToGlobalMappings();
    }
    char matType[30];
    PetscTruth typeFound;
    PetscOptionsGetString(PETSC_NULL,"-fullJacMatType",matType,30,&typeFound);
    if(!typeFound) {
      std::cout<<"I need a MatType for the full matrix!"<<std::endl;
      assert(false);
    } 
    da->createMatrix(*jac, matType, 3);
    if(!myRank) {
      std::cout<<"Finished Building Full vecMass Mat."<<std::endl;
    }
  }else {  
    //The size this processor owns ( without ghosts).
    unsigned int  m,n;
    m=n=(3*(da->getNodeSize()));
    MatCreateShell(damg->comm, m ,n, PETSC_DETERMINE, PETSC_DETERMINE, damg, jac);
    MatShellSetOperation(*jac ,MATOP_MULT, (void (*)(void)) ConstVecMassMatMult);
    MatShellSetOperation(*jac ,MATOP_GET_DIAGONAL, (void (*)(void)) ConstVecMassMatGetDiagonal);
    MatShellSetOperation(*jac ,MATOP_DESTROY, (void (*)(void)) ConstVecMassMatDestroy);
  }
  PetscFunctionReturn(0);
}//end fn.

PetscErrorCode ConstVecMassMatDestroy(Mat J) {
  PetscFunctionBegin;
  //Nothing to be done here. No new pointers were created during creation. 
  PetscFunctionReturn(0);
}

#define VEC_MASS_ELEM_DIAG_BLOCK {\
  unsigned int idx = da->curr();\
  unsigned int lev = da->getLevel(idx);\
  double h = hFac*(1u << (maxD - lev));\
  double fac = h*h*h/8.0;\
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
        diagArr[(3*indices[k])+dof] += (fac*(MassType2Stencil[childNum][elemType][k][k]));\
      } /*end dof*/\
    }\
  } /*end k*/\
}

#define VEC_MASS_DIAG_BLOCK {\
  ot::DA* da = damg->da;\
  iC(VecZeroEntries(diag));\
  PetscScalar *diagArr;\
  unsigned char* bdyArr = data->bdyArr;\
  /*Nodal,Non-Ghosted,Write,3 dof*/\
  da->vecGetBuffer(diag, diagArr, false, false, false, 3);\
  unsigned int maxD;\
  double hFac;\
  if(da->iAmActive()) {\
    maxD = (da->getMaxDepth());\
    hFac = 1.0/((double)(1u << (maxD-1)));\
    /*Loop through All Elements including ghosted*/\
    for(da->init<ot::DA_FLAGS::ALL>(); da->curr() < da->end<ot::DA_FLAGS::ALL>(); da->next<ot::DA_FLAGS::ALL>()) {\
      VEC_MASS_ELEM_DIAG_BLOCK\
    } /*end i*/\
  } /*end if active*/\
  da->vecRestoreBuffer(diag, diagArr, false, false, false, 3);\
  /*2 IOP = 1 FLOP. Loop counters are included too.*/\
  PetscLogFlops(93*(da->getGhostedElementSize()));\
}

PetscErrorCode ConstVecMassMatGetDiagonal(Mat J, Vec diag) {
  PetscFunctionBegin;

  PetscLogEventBegin(vecMassDiagEvent,diag,0,0,0);

  ot::DAMG damg;
  iC(MatShellGetContext(J, (void**)(&damg)));

  ElasticityData* data = (static_cast<ElasticityData*>(damg->user));

  bool isFinestLevel = (damg->nlevels == 1);

  if(isFinestLevel) {
    PetscLogEventBegin(vecMassFinestDiagEvent,diag,0,0,0);
  }

  VEC_MASS_DIAG_BLOCK 

    if(isFinestLevel) {
      PetscLogEventEnd(vecMassFinestDiagEvent,diag,0,0,0);
    }

  PetscLogEventEnd(vecMassDiagEvent,diag,0,0,0);

  PetscFunctionReturn(0);
}

#define VEC_MASS_ELEM_MULT_BLOCK {\
  unsigned int idx = da->curr();\
  unsigned int lev = da->getLevel(idx);\
  double h = hFac*(1u << (maxD - lev));\
  double fac = h*h*h/8.0;\
  unsigned int indices[8];\
  da->getNodeIndices(indices);\
  unsigned char childNum = da->getChildNumber();\
  unsigned char hnMask = da->getHangingNodeIndex(idx);\
  unsigned char elemType = 0;\
  GET_ETYPE_BLOCK(elemType, hnMask, childNum)\
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
            outArr[(3*indices[k]) + dof] += (fac*(MassType2Stencil[childNum][elemType][k][j])\
                *inArr[(3*indices[j]) + dof]);\
          }/*end for dof*/\
        }\
      }/*end for j*/\
    }\
  }/*end for k*/\
}

#define VEC_MASS_MULT_BLOCK {\
  ot::DA* da = damg->da;\
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
  /*Nodal,Non-Ghosted,Read,3 dof*/\
  da->vecGetBuffer(in, inArr, false, false, true, 3);\
  if(da->iAmActive()) {\
    da->ReadFromGhostsBegin<PetscScalar>(inArr, 3);\
  }\
  /*Nodal,Non-Ghosted,Write,3 dof*/\
  da->vecGetBuffer(out, outArr, false, false, false, 3);\
  /*Loop through All Independent Elements*/\
  if(da->iAmActive()) {\
    for(da->init<ot::DA_FLAGS::INDEPENDENT>();\
        da->curr() < da->end<ot::DA_FLAGS::INDEPENDENT>(); da->next<ot::DA_FLAGS::INDEPENDENT>() ) {\
      VEC_MASS_ELEM_MULT_BLOCK\
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
        da->curr() < da->end<ot::DA_FLAGS::DEPENDENT>();  da->next<ot::DA_FLAGS::DEPENDENT>() ) {\
      VEC_MASS_ELEM_MULT_BLOCK\
    } /*end loop for dependent elems*/\
  } /*end if active*/\
  da->vecRestoreBuffer(in, inArr, false, false, true, 3);\
  da->vecRestoreBuffer(out, outArr, false, false, false, 3);\
  /*2 IOP = 1 FLOP. Loop counters are included too.*/\
  PetscLogFlops(1095*(da->getGhostedElementSize()));\
}

PetscErrorCode ConstVecMassMatMult(Mat J, Vec in, Vec out)
{
  PetscFunctionBegin;

  PetscLogEventBegin(vecMassMultEvent,in,out,0,0);

  ot::DAMG damg;
  iC(MatShellGetContext(J, (void**)(&damg)));

  ElasticityData* data = (static_cast<ElasticityData*>(damg->user));

  bool isFinestLevel = (damg->nlevels == 1);

  if(isFinestLevel) {
    PetscLogEventBegin(vecMassFinestMultEvent,in,out,0,0);
  }

  VEC_MASS_MULT_BLOCK 

    if(isFinestLevel) {
      PetscLogEventEnd(vecMassFinestMultEvent,in,out,0,0);
    }

  PetscLogEventEnd(vecMassMultEvent,in,out,0,0);

  PetscFunctionReturn(0);
}

#define BUILD_FULL_VEC_MASS_BLOCK {\
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
  if(da->iAmActive()) {\
    for(da->init<ot::DA_FLAGS::WRITABLE>();\
        da->curr() < da->end<ot::DA_FLAGS::WRITABLE>();  da->next<ot::DA_FLAGS::WRITABLE>()) {\
      unsigned int idx = da->curr();\
      unsigned int lev = da->getLevel(idx);\
      double h = hFac*(1u << (maxD - lev));\
      double fac = h/2.0;\
      unsigned int indices[8];\
      da->getNodeIndices(indices);\
      unsigned char childNum = da->getChildNumber();\
      unsigned char hnMask = da->getHangingNodeIndex(idx);\
      unsigned char elemType = 0;\
      GET_ETYPE_BLOCK(elemType, hnMask, childNum)\
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
                currRec.val = (fac*(MassType2Stencil[childNum][elemType][k][j]));\
                records.push_back(currRec);\
              } /*end for dof*/\
            }\
          } /*end for j*/\
        }\
      } /*end for k*/\
      if(records.size() > 1000) {\
        /*records will be cleared inside the function*/\
        da->setValuesInMatrix(B, records, 3, ADD_VALUES);\
      }\
    } /*end writable*/\
  } /*end if active*/\
  da->setValuesInMatrix(B, records, 3, ADD_VALUES);\
  MatAssemblyBegin(B, MAT_FLUSH_ASSEMBLY);\
  MatAssemblyEnd(B, MAT_FLUSH_ASSEMBLY);\
  if(da->iAmActive()) {\
    for(da->init<ot::DA_FLAGS::WRITABLE>();\
        da->curr() < da->end<ot::DA_FLAGS::WRITABLE>();  da->next<ot::DA_FLAGS::WRITABLE>()) {\
      unsigned int idx = da->curr();\
      unsigned int lev = da->getLevel(idx);\
      double h = hFac*(1u << (maxD - lev));\
      double fac = h/2.0;\
      unsigned int indices[8];\
      da->getNodeIndices(indices);\
      unsigned char childNum = da->getChildNumber();\
      unsigned char hnMask = da->getHangingNodeIndex(idx);\
      unsigned char elemType = 0;\
      GET_ETYPE_BLOCK(elemType, hnMask, childNum)\
      for(int k = 0; k < 8; k++) {\
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
    } /*end writable*/\
  } /*end if active*/\
  da->setValuesInMatrix(B, records, 3, INSERT_VALUES);\
  MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);\
  MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);\
}

PetscErrorCode ComputeConstVecMass(ot::DAMG damg, Mat J, Mat B) {
  //For matShells nothing to be done here.
  PetscFunctionBegin;

  PetscTruth isshell;
  PetscTypeCompare((PetscObject)B, MATSHELL, &isshell);
  if(isshell) {
    PetscFunctionReturn(0);
  }

  ElasticityData* data = (static_cast<ElasticityData*>(damg->user));

  //Assuming that B and J are the same.

  BUILD_FULL_VEC_MASS_BLOCK 

    PetscFunctionReturn(0);
}

#undef BUILD_FULL_VEC_MASS_BLOCK 
#undef VEC_MASS_MULT_BLOCK 
#undef VEC_MASS_DIAG_BLOCK 
#undef VEC_MASS_ELEM_DIAG_BLOCK 
#undef VEC_MASS_ELEM_MULT_BLOCK 


