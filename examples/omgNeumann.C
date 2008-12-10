/**
  @file omgNeumann.C
  @author Ilya Lashuk, ilya.lashuk@gmail.com
  */
#include <mpi.h>
#include <cstdio>
#include "oda.h"
#include "omg.h"
#include "Point.h"
#include "parUtils.h"
#include "octUtils.h"
#include "TreeNode.h"
#include "handleStencils.h"
#include <cstdlib>
#include "sys.h"
#include <sstream>
#include "omgNeumann.h"
#include "dendro.h"

#ifdef PETSC_USE_LOG
extern int Jac2MultEvent;
extern int Jac2DiagEvent;
extern int Jac2FinestMultEvent;
extern int Jac2FinestDiagEvent;
#endif

static double**** LaplacianType2Stencil; 
static double**** MassType2Stencil; 
static double*** RHSType2Stencil; 
static std::vector<double> force_values; // value of RHS at centers of "ALL" octants on this process;  global since ComputeRHS needs it

struct Jac2MFreeData {
  std::vector<double>* matProp;
  bool isFinestLevel;
  Mat Jmat_private;
  Vec inTmp;
  Vec outTmp;
};

#define CHK_AND_SCATTER_FINE_TO_COARSE_BLOCK {\
  if( damg_i->da_aux == NULL ) {\
    daf = da;\
    fMatPropVec = matPropVecPtr;\
    changedPartition  = false;\
  }else {\
    daf = damg_i->da_aux;\
    /*Need to Scatter Values*/\
    fMatPropVec = new std::vector<double>;\
    changedPartition  = true;\
    /*elemental - non-ghosted*/\
    std::vector<double> tmpVec1;\
    da->createVector<double>(tmpVec1,true,false,2);\
    double *vec1Arr = NULL;\
    /*Elemental,non-Ghosted,Write,2 Dof.*/\
    da->vecGetBuffer<double>(tmpVec1,vec1Arr,true,false,false,2);\
    matPropArr = NULL;\
    /*Elemental,Ghosted,Read-only,2 Dof.*/\
    da->vecGetBuffer<double>((*matPropVecPtr), matPropArr, true, true, true, 2);\
    if(da->iAmActive()) {\
      for(da->init<ot::DA_FLAGS::WRITABLE>(); da->curr() < da->end<ot::DA_FLAGS::WRITABLE>(); da->next<ot::DA_FLAGS::WRITABLE>()) {\
        unsigned int idx = da->curr();\
        vec1Arr[2*idx] = matPropArr[2*idx];\
        vec1Arr[2*idx+1] = matPropArr[2*idx+1];\
      }\
    }\
    da->vecRestoreBuffer<double>((*matPropVecPtr), matPropArr, true, true, true, 2);\
    da->vecRestoreBuffer<double>(tmpVec1,vec1Arr,true,false,false,2);\
    par::scatterValues<double>(tmpVec1, (*fMatPropVec), (2*(daf->getElementSize())), da->getComm());\
    tmpVec1.clear();\
  }\
}

#define ASSIGN_MAT_PROP_FINE_TO_COARSE_ELEM_BLOCK {\
  /*The fine loop is always Writable, but the coarse loop*/\
  /*could be Independent or W_Dependent. Hence the fine counter must*/\
  /*be incremented properly to align with the coarse.*/\
  unsigned int idxC= da->curr();\
  Point Cpt = da->getCurrentOffset();\
  assert(daf->curr() < daf->end<ot::DA_FLAGS::WRITABLE>());\
  while(daf->getCurrentOffset() != Cpt) {\
    daf->next<ot::DA_FLAGS::WRITABLE>();\
    assert(daf->curr() < daf->end<ot::DA_FLAGS::WRITABLE>());\
  }\
  if(daf->getLevel(daf->curr()) == da->getLevel(idxC)) {\
    /*The coarse and fine elements are the same,*/\
    assert(daf->curr() < daf->end<ot::DA_FLAGS::WRITABLE>());\
    unsigned int idxF = daf->curr();\
    matPropArr[2*idxC] = fMatArr[2*idxF];\
    matPropArr[2*idxC+1] = fMatArr[2*idxF+1];\
    daf->next<ot::DA_FLAGS::WRITABLE>();\
  }else {\
    double cLapVal = 0.0;\
    double cMassVal = 0.0;\
    for(unsigned char cNumFine = 0; cNumFine < 8; cNumFine++) {\
      /*The coarse and fine elements are NOT the same. */\
      /*Loop over each of the 8 children of the coarse element.*/\
      /*These are the underlying fine elements.*/\
      assert(daf->curr() < daf->end<ot::DA_FLAGS::WRITABLE>());\
      unsigned int idxF = daf->curr();\
      cLapVal += fMatArr[2*idxF];\
      cMassVal += fMatArr[2*idxF+1];\
      daf->next<ot::DA_FLAGS::WRITABLE>();\
    }\
    matPropArr[2*idxC] = (cLapVal/8.0);\
    matPropArr[2*idxC+1] = (cMassVal/8.0);\
  }\
  if(matPropArr[2*idxC] > maxCoeff) {\
    maxCoeff = matPropArr[2*idxC];\
  }\
  if(matPropArr[2*idxC] < minCoeff) {\
    minCoeff = matPropArr[2*idxC];\
  }\
}

#define FINE_TO_COARSE_BLOCK {\
  damg_i = damg[i];\
  damg_i->user = ctx;\
  da = damg_i->da;\
  assert(da->iAmActive() == daf->iAmActive());\
  da->createVector<double>((*matPropVecPtr), true, true, 2);\
  for(unsigned int j = 0; j < matPropVecPtr->size(); j++) {\
    (*matPropVecPtr)[j] = 0.0;\
  }\
  comm = da->getCommActive();\
  MPI_Comm_rank(comm,&rank);\
  /*Elemental,Ghosted,Write,2 Dof.*/\
  matPropArr = NULL;\
  da->vecGetBuffer<double>((*matPropVecPtr), matPropArr,\
      true, true, false, 2);\
  double *fMatArr = NULL;\
  if(changedPartition) {\
    /*Elemental, non-Ghosted, Read-only, 2 Dof.*/\
    daf->vecGetBuffer<double>((*fMatPropVec), fMatArr,\
        true, false, true, 2);\
  }else {\
    /*Elemental, Ghosted, Read-only, 2 Dof.*/\
    daf->vecGetBuffer<double>((*fMatPropVec), fMatArr,\
        true, true, true, 2);\
  }\
  if(da->iAmActive()) {\
    double maxCoeff = 0.0;\
    double minCoeff = 1.0e+8;\
    double globalMaxCoeff;\
    double globalMinCoeff;\
    /*Loop through the coarse and fine simultaneously*/\
    /*Note: If Coarse is Independent, then the*/\
    /*corresponding Fine is also independent.*/\
    /*Hence, overlapping comm with comp is possible.*/\
    /*First, we loop though the dependent elements.*/\
    /*Then we begin the communication and simulatenously*/\
    /*loop over the independent elements.*/\
    for(da->init<ot::DA_FLAGS::W_DEPENDENT>(),\
        daf->init<ot::DA_FLAGS::WRITABLE>();\
        da->curr() < da->end<ot::DA_FLAGS::W_DEPENDENT>();\
        da->next<ot::DA_FLAGS::W_DEPENDENT>()) {\
      ASSIGN_MAT_PROP_FINE_TO_COARSE_ELEM_BLOCK \
    } /*end dependent loop*/\
    da->ReadFromGhostElemsBegin<double>(matPropArr,2);\
    for(da->init<ot::DA_FLAGS::INDEPENDENT>(),\
        daf->init<ot::DA_FLAGS::WRITABLE>();\
        da->curr() < da->end<ot::DA_FLAGS::INDEPENDENT>();\
        da->next<ot::DA_FLAGS::INDEPENDENT>()) {\
      ASSIGN_MAT_PROP_FINE_TO_COARSE_ELEM_BLOCK \
    } /*end Independent loop */\
    da->ReadFromGhostElemsEnd<double>(matPropArr);\
    par::Mpi_Reduce<double>(&maxCoeff, &globalMaxCoeff, 1, MPI_MAX, 0, comm);\
    par::Mpi_Reduce<double>(&minCoeff, &globalMinCoeff, 1, MPI_MIN, 0, comm);\
    if(!rank) {\
      std::cout<<"Level: "<<i<<" Max Lap. Coeff: "\
      <<globalMaxCoeff<<" Min Lap. Coeff: "\
      <<globalMinCoeff<<std::endl;\
    }\
    MPI_Barrier(comm);\
  } /*end check if active*/\
  da->vecRestoreBuffer<double>((*matPropVecPtr),\
      matPropArr, true, true, false, 2);\
  if(changedPartition) {\
    /*Elemental, non-Ghosted, Read-only, 2 Dof.*/\
    daf->vecRestoreBuffer<double>((*fMatPropVec),\
        fMatArr, true, false, true, 2);\
  }else {\
    /*Elemental, Ghosted, Read-only, 2 Dof.*/\
    daf->vecRestoreBuffer<double>((*fMatPropVec),\
        fMatArr, true, true, true, 2);\
  }\
}

void SetPDECoefFromPts(
    ot::DAMG* damg,
    const std::vector<double>& centers,
    void (*CalcCoef)(const std::vector<double> & pts, std::vector<double> & values)
    )
{
  int nlevels = damg[0]->nlevels; //number of mg levels
  ot::DAMG damg_i = damg[nlevels-1];

  //Set Mat Props Finest to Coarsest...
  void * ctx;

  // Set for the finest level first
  ctx = new Jac2MFreeData;

  std::vector<double> * matPropVecPtr = new std::vector<double>;

  (static_cast<Jac2MFreeData*>(ctx))->matProp = matPropVecPtr;
  (static_cast<Jac2MFreeData*>(ctx))->isFinestLevel = true;
  (static_cast<Jac2MFreeData*>(ctx))->Jmat_private = NULL;
  (static_cast<Jac2MFreeData*>(ctx))->inTmp = NULL;
  (static_cast<Jac2MFreeData*>(ctx))->outTmp = NULL;

  double *matPropArr = NULL;
  ot::DA* da;
  int rank;
  MPI_Comm comm;

  damg_i->user = ctx;
  da = damg_i->da;
  comm = da->getCommActive();
  MPI_Comm_rank(comm,&rank);
  /*Elem,Ghosted, 2-dof vector.*/
  /*Note: I am creating a ghosted vector only*/
  /*because the mat-vec will need it.*/
  /*So this way, I can avoid mallocs inside the mat-vec.*/
  da->createVector<double>((*matPropVecPtr), true, true, 2);
  for(unsigned int i = 0; i < matPropVecPtr->size(); i++) {
    (*matPropVecPtr)[i] = 0.0;
  }
  /*Elemental,Ghosted,Write,2 Dof.*/
  da->vecGetBuffer<double>((*matPropVecPtr), matPropArr,
      true, true, false, 2);
  if(da->iAmActive()) {
    double maxCoeff = 0.0;
    double minCoeff = 1.0e+80;
    double globalMaxCoeff;
    double globalMinCoeff;

    // call user-provided routine to calculate the coefficients
    // we can use here *matPropVectPtr since this vector is both elemental and ghosted
    CalcCoef(centers, *matPropVecPtr);

    // calculate min and max "alpha" for debugging
    for(da->init<ot::DA_FLAGS::ALL>();
	da->curr() < da->end<ot::DA_FLAGS::ALL>();
	da->next<ot::DA_FLAGS::ALL>()) {
      unsigned int idx = da->curr();

      if(matPropArr[2*idx] > maxCoeff) 
	maxCoeff = matPropArr[2*idx];

      if(matPropArr[2*idx] < minCoeff) 
	minCoeff = matPropArr[2*idx];
    }

    par::Mpi_Reduce<double>(&maxCoeff, &globalMaxCoeff, 1, MPI_MAX, 0, comm);
    par::Mpi_Reduce<double>(&minCoeff, &globalMinCoeff, 1, MPI_MIN, 0, comm);
    if(!rank) {
      std::cout<<"Level: "<<(nlevels-(damg_i->nlevels))<<
	" Max Lap. Coeff: "<<globalMaxCoeff<<
	" Min Lap. Coeff: "<<globalMinCoeff<<std::endl;
    }
  } /*end if active*/
  da->vecRestoreBuffer<double>((*matPropVecPtr), matPropArr,
      true, true, false, 2);

  //The coarser levels...
  ot::DA* daf;
  std::vector<double>* fMatPropVec = NULL;
  bool changedPartition;

  if(nlevels > 1) {
    CHK_AND_SCATTER_FINE_TO_COARSE_BLOCK 
  }

  //Coarser levels
  for(int i = (nlevels-2); i >= 0; i--) {
      ctx = new Jac2MFreeData;

    matPropVecPtr = new std::vector<double>;

    (static_cast<Jac2MFreeData*>(ctx))->matProp = matPropVecPtr;
    (static_cast<Jac2MFreeData*>(ctx))->isFinestLevel = false;
    (static_cast<Jac2MFreeData*>(ctx))->Jmat_private = NULL;
    (static_cast<Jac2MFreeData*>(ctx))->inTmp = NULL;
    (static_cast<Jac2MFreeData*>(ctx))->outTmp = NULL;

    FINE_TO_COARSE_BLOCK

      if(changedPartition) {
	fMatPropVec->clear();
	delete fMatPropVec;
      }

    if(i) {	
      CHK_AND_SCATTER_FINE_TO_COARSE_BLOCK 
    }
  }//end for i
}//end fn.


void DestroyUserContexts(ot::DAMG* damg) {
  
  int       nlevels = damg[0]->nlevels; //number of multigrid levels

  for(int i = 0; i < nlevels; i++) {
    Jac2MFreeData* ctx = (static_cast<Jac2MFreeData*>(damg[i]->user));
    ctx->matProp->clear();
    delete ctx->matProp;
    if(ctx->Jmat_private) {
      MatDestroy(ctx->Jmat_private);
      ctx->Jmat_private = NULL;
    }
    if(ctx->inTmp) {
      VecDestroy(ctx->inTmp);
      ctx->inTmp = NULL;
    }
    if(ctx->outTmp) {
      VecDestroy(ctx->outTmp);
      ctx->outTmp = NULL;
    }
    delete ctx;
  }
}//end fn.

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
    for(da->init<ot::DA_FLAGS::ALL>(); da->curr() < da->end<ot::DA_FLAGS::ALL>(); da->next<ot::DA_FLAGS::ALL>()) {\
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
  da->ReadFromGhostsBegin<PetscScalar>(inArr,1);\
  /*Nodal,Non-Ghosted,Write,1 dof*/\
  da->vecGetBuffer(out,outArr,false,false,false,1);\
  /*Loop through All Independent Elements*/\
  if(da->iAmActive()) {\
    for(da->init<ot::DA_FLAGS::INDEPENDENT>(); da->curr() < da->end<ot::DA_FLAGS::INDEPENDENT>(); da->next<ot::DA_FLAGS::INDEPENDENT>() ) {\
      JAC_TYPE2_ELEM_MULT_BLOCK\
    } /*end independent*/\
  } /*end if active*/\
  da->ReadFromGhostsEnd<PetscScalar>(inArr);\
  /*Loop through All Dependent Elements,*/\
  /*i.e. Elements which have atleast one*/\
  /*vertex owned by this processor and at least one*/\
  /*vertex not owned by this processor.*/\
  if(da->iAmActive()) {\
    for(da->init<ot::DA_FLAGS::DEPENDENT>(); da->curr() < da->end<ot::DA_FLAGS::DEPENDENT>();  da->next<ot::DA_FLAGS::DEPENDENT>() ) {\
      JAC_TYPE2_ELEM_MULT_BLOCK\
    } /*end loop for dependent elems*/\
  } /*end if active*/\
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
    for(da->init<ot::DA_FLAGS::WRITABLE>(); da->curr() < da->end<ot::DA_FLAGS::WRITABLE>();  da->next<ot::DA_FLAGS::WRITABLE>()) {\
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

  //Assuming that B and J are the same.

  BUILD_FULL_JAC_TYPE2_BLOCK(B) 

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

/*
 * This function computes the right-hand side for FEM system. The load (force) is sampled at the center of each octant. These sample values are read from the global variable force_values */
PetscErrorCode ComputeRHS_omgNeumann(ot::DAMG damg,Vec in) {
  PetscFunctionBegin;            
  ot::DA* da = damg->da;
  PetscScalar *inarray;
  VecZeroEntries(in);
  da->vecGetBuffer(in,inarray,false,false,false,1);

  if(da->iAmActive()) {
    for(da->init<ot::DA_FLAGS::ALL>();
	da->curr() < da->end<ot::DA_FLAGS::ALL>();
	da->next<ot::DA_FLAGS::ALL>())  
    {
      unsigned int idx = da->curr();
      unsigned levelhere = (da->getLevel(idx) - 1);
      
      double hxOct = ldexp(1.0,-levelhere);
      assert(hxOct==1.0/(1u<<levelhere));
      double fac = ((hxOct*hxOct*hxOct)/8.0);
      
      unsigned int indices[8];
      da->getNodeIndices(indices); 
      
      unsigned char childNum = da->getChildNumber();
      unsigned char hnMask = da->getHangingNodeIndex(idx);
      unsigned char elemType = 0;
      GET_ETYPE_BLOCK(elemType,hnMask,childNum)
      
      for(unsigned int j = 0; j < 8; j++) 
	inarray[indices[j]] += force_values[idx]*RHSType2Stencil[childNum][elemType][j]*fac;
    
    }//end looping over octants
    
  }//end if active

  da->vecRestoreBuffer(in,inarray,false,false,false,1);
  PetscFunctionReturn(0);
}


static void CalculateCenters(ot::DA* da, std::vector<double> & centers)
{
  // this only resizes "centers"
  // also, since vector is both elemental and ghosted, we don't need the "getbuffer" function for access
  da->createVector<double>(centers, true/*elemental*/, true/*ghosted*/, 3/*values per octant*/);

  // initialize all entries to zero (loop over elements might not touch some entries)
  // maybe this can be skipped. but then some function will have to calculate garbage from garbage
  centers.assign(centers.size(),0.0);
  
  if(da->iAmActive()) {
    unsigned maxD = da->getMaxDepth();
    
    for(da->init<ot::DA_FLAGS::ALL>();
	da->curr() < da->end<ot::DA_FLAGS::ALL>();
	da->next<ot::DA_FLAGS::ALL>())  
    {
      Point pt = da->getCurrentOffset();
      unsigned int idx = da->curr();
      
      double x = ldexp( static_cast<double>(pt.xint()) , 1-maxD );
      double y = ldexp( static_cast<double>(pt.yint()) , 1-maxD );
      double z = ldexp( static_cast<double>(pt.zint()) , 1-maxD );
      unsigned levelhere = (da->getLevel(idx) - 1);
      double halfSide = ldexp(0.5, -levelhere);
      
      centers[3*idx]=x+halfSide;
      centers[3*idx+1]=y+halfSide;
      centers[3*idx+2]=z+halfSide;
    }//end looping over octants
    
  }//end if active
}

/**
 * @author Ilya Lashuk
 * @brief This function calculates the approximate solution to the scalar equation  -div(alpha*grad u) + beta * u = f in the unit cube subject to homogeneous Neumann boundary conditions. 
 * @param pts coordinates of points. Each point must be in [0,1)^3.  Octree will be built based on these points. No octant of the octree will contain more than 1 point. Coordinates are stored in a sequence X1, Y1, Z1, X2, Y2, Z2, ...
 * @param CalcCoef function pointer to function which will calculate the values of "alpha" and "beta" (see the equation above) for the center of each octant. Coordinates of all centers are provided in the first parameter (layout is X1,Y1,Z1,X2,Y2,Z2,...). Values should be returned in second parameter. The following layout should be used: alpha1,beta1,alpha2,beta2,.... Caller will set the array size for the second parameter (to be equal to 2/3 size of the first parameter). 
 * @param  CalcRHS function pointer to function which will calculate the value of "f" (see the equation above) for the center of each octant. Coordinates of all centers are provided in the first parameter (layout is the same as above). Values of "f" should be returned in second parameter. Caller will set the array size for the second parameter (to be equal to 1/3 size of the first parameter).
 * @param numMultigridLevels desired number of multigrid levels.
 * @param sol in this parameter the solution is returned. This object is automatically destroyed when the "damg" object (see next parameter) is destroyed.
 * @param damg in this parameter the multigrid context (used to calculate solution) is returned. User must destroy this object when no longer needed by calling DAMGDestroy
 */
void solve_neumann(
    /* input parameters: */
     std::vector<double>& pts,
     void (*CalcCoef)(const std::vector<double> & pts, std::vector<double> & values),
     void (*CalcRHS)(const std::vector<double> & pts, std::vector<double> & values),
     int numMultigridLevels,
     /* output parameters */
     Vec& sol,
     ot::DAMG * & damg
    )
{
  using namespace std;
  
  const int dim = 3;
  double gSize[3]={1.0, 1.0, 1.0};
  const double mgLoadFac = 1.5;
  const bool compressLut = false;
  const bool incCorner = true; 
  const int maxNumPts=1;  // we want at most 1 point per octant
  const int maxDepth = 30;
  const int dof=1;

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  
  // enforce that all points are inside [0,1)^3
  for(size_t i=0;i<pts.size();i++) 
    if (pts[i]<0)
      pts[i]=0;
    else if (pts[i]>=1)
      pts[i]=1-ldexp(0.5,-maxDepth);

  // now convert points to octree 
  // "pts" is cleared inside this function
  vector<ot::TreeNode> linOct;
  ot::points2Octree(pts, gSize, linOct, dim, maxDepth, maxNumPts, MPI_COMM_WORLD);

  //now balance the octree; "linOct" is cleared inside this function
  vector<ot::TreeNode> balOct;
  ot::balanceOctree (linOct, balOct, dim, maxDepth, true, MPI_COMM_WORLD, NULL, NULL);
  
  /*
  for(int i = 0; i < 1; i++) {
    std::vector<ot::TreeNode> tmpOct = balOct;
    balOct.clear();
    ot::refineOctree(tmpOct, balOct); 
  }
  */

  // createRegularOctree(balOct,5,3,maxDepth,MPI_COMM_WORLD);

  // print how many finest-level octants we have
  DendroIntL locBalSize = balOct.size();
  DendroIntL totBalSize;
  par::Mpi_Reduce<DendroIntL>(&locBalSize, &totBalSize, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  if(!rank) {
    cout << "# of octants on finest level: "<< totBalSize << endl; 
  }
  
  // create multigrid solver object, coarser octrees, finite element meshes and the interpolation
  // balOct is cleared inside this function
  // numMultigridLevels might be modified inside this function
  ot::DAMGCreateAndSetDA(MPI_COMM_WORLD,numMultigridLevels, NULL, &damg, balOct,
      dof, mgLoadFac, compressLut, incCorner);
 
  // set some options to make "CreateJacobian2" and "ComputeJacobian2" work. These routines are related to matrix-vector multiplication on all levels.
  PetscOptionsSetValue("-jacType","2");
  
  // load stencils, i.e. integrals of: products of gradients of shape functions, product of shape functions and shape functions themselves
  createLmatType2(LaplacianType2Stencil);
  createMmatType2(MassType2Stencil);
  createRHSType2(RHSType2Stencil);

  // here will be stored centers of "ALL" octants on this process;
  vector<double> centers;   
  
  // calculate centers for "ALL" octants (both pre-ghost and writable)
  CalculateCenters(damg[numMultigridLevels-1]->da, centers);
  
  // call user-provided function to sample the right-hand-side
  force_values.resize(centers.size()/3);
  CalcRHS(centers, force_values);
  
  // calculate PDE coef. for all levels. This calls user-provided function CalcCoef
  SetPDECoefFromPts(damg, centers,CalcCoef); 
  
  // set up functions which do matvecs and compute right hand side
  ot::DAMGSetKSP(damg, CreateJacobian2,ComputeJacobian2,ComputeRHS_omgNeumann);
  
  // solve (using zero initial guess -- this is the default)
  ot::DAMGSolve(damg);

  // get the solution
  sol = DAMGGetx(damg);
  
  // destroy objects 
  destroyLmatType2(LaplacianType2Stencil);
  destroyMmatType2(MassType2Stencil);
  destroyRHSType2(RHSType2Stencil);
  DestroyUserContexts(damg);

  return;
}

/**
 * @author Ilya Lashuk
 * @brief This function calculates the approximate solution to the scalar equation  -div(alpha*grad u) + beta * u = f in the unit cube subject to homogeneous Neumann boundary conditions. 
 * @param octs sequence of octants. Octree will be built by completing this sequence of octants.
 * @param CalcCoef function pointer to function which will calculate the values of "alpha" and "beta" (see the equation above) for the center of each octant. Coordinates of all centers are provided in the first parameter (layout is X1,Y1,Z1,X2,Y2,Z2,...). Values should be returned in second parameter. The following layout should be used: alpha1,beta1,alpha2,beta2,.... Caller will set the array size for the second parameter (to be equal to 2/3 size of the first parameter). 
 * @param  CalcRHS function pointer to function which will calculate the value of "f" (see the equation above) for the center of each octant. Coordinates of all centers are provided in the first parameter (layout is the same as above). Values of "f" should be returned in second parameter. Caller will set the array size for the second parameter (to be equal to 1/3 size of the first parameter).
 * @param numMultigridLevels desired number of multigrid levels.
 * @param sol in this parameter the solution is returned. This object is automatically destroyed when the "damg" object (see next parameter) is destroyed.
 * @param damg in this parameter the multigrid context (used to calculate solution) is returned. User must destroy this object when no longer needed by calling DAMGDestroy
 */

void solve_neumann_oct(
    /* input parameters: */
     std::vector<ot::TreeNode>& octs,
     void (*CalcCoef)(const std::vector<double> & pts, std::vector<double> & values),
     void (*CalcRHS)(const std::vector<double> & pts, std::vector<double> & values),
     int numMultigridLevels,
     /* output parameters */
     Vec& sol,
     ot::DAMG * & damg
    )
{
  using namespace std;
  
  const int dim = 3;
  double gSize[3]={1.0, 1.0, 1.0};
  const double mgLoadFac = 1.5;
  const bool compressLut = false;
  const bool incCorner = true; 
  const int maxDepth = octs[0].getMaxDepth();
  const int dof=1;


  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  // complete the octree 
  vector<ot::TreeNode> linOct;
  ot::completeOctree(octs, linOct, 3 /*dim*/,
      maxDepth,
      true, /* isUnique */
      false, /* isSorted */
      false, /* assertNoEmptyProcs */
      MPI_COMM_WORLD);
  octs.clear();
    
  //now balance the octree; "linOct" is cleared inside this function
  vector<ot::TreeNode> balOct;
  ot::balanceOctree (linOct, balOct, dim, maxDepth, true, MPI_COMM_WORLD, NULL, NULL);
  
  // debug 
  // balOct.clear();
  // createRegularOctree(balOct,7,3,maxDepth,MPI_COMM_WORLD);
  
  // debug -- print octree to file
  // writeNodesToFile ("the_tree.ot", balOct);
  
  // print how many finest-level octants we have
  DendroIntL locBalSize = balOct.size();
  DendroIntL totBalSize;
  par::Mpi_Reduce<DendroIntL>(&locBalSize, &totBalSize, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  if(!rank) {
    cout << "# of octants on finest level: "<< totBalSize << endl; 
  }
  
  // create multigrid solver object, coarser octrees, finite element meshes and the interpolation
  // balOct is cleared inside this function
  // numMultigridLevels might be modified inside this function
  ot::DAMGCreateAndSetDA(MPI_COMM_WORLD,numMultigridLevels, NULL, &damg, balOct,
      dof, mgLoadFac, compressLut, incCorner);

  // set some options to make "CreateJacobian2" and "ComputeJacobian2" work. These routines are related to matrix-vector multiplication on all levels.
  PetscOptionsSetValue("-jacType","2");
  
  // load stencils, i.e. integrals of: products of gradients of shape functions, product of shape functions and shape functions themselves
  createLmatType2(LaplacianType2Stencil);
  createMmatType2(MassType2Stencil);
  createRHSType2(RHSType2Stencil);

  // here will be stored centers of "ALL" octants on this process;
  vector<double> centers;   
  
  // calculate centers for "ALL" octants (both pre-ghost and writable)
  CalculateCenters(damg[numMultigridLevels-1]->da, centers);
  
  // call user-provided function to sample the right-hand-side
  force_values.resize(centers.size()/3);
  CalcRHS(centers, force_values);
  
  // calculate PDE coef. for all levels. This calls user-provided function CalcCoef
  SetPDECoefFromPts(damg, centers,CalcCoef); 
  
  // set up functions which do matvecs and compute right hand side
  ot::DAMGSetKSP(damg, CreateJacobian2,ComputeJacobian2,ComputeRHS_omgNeumann);
  
  // solve (using zero initial guess -- this is the default)
  ot::DAMGSolve(damg);

  // get the solution
  sol = DAMGGetx(damg);
  
  // destroy objects 
  destroyLmatType2(LaplacianType2Stencil);
  destroyMmatType2(MassType2Stencil);
  destroyRHSType2(RHSType2Stencil);
  DestroyUserContexts(damg);

  return;
}

