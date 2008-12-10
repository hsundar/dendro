
/**
  @file odaJac.C
  @author Rahul Sampath, rahul.sampath@gmail.com
  */

#include "petsc.h"
#include "petscvec.h"
#include "petscmat.h"
#include "odaJac.h"

#ifdef PETSC_USE_LOG
extern int Jac1DiagEvent;
extern int Jac1MultEvent;
extern int Jac1FinestDiagEvent;
extern int Jac1FinestMultEvent;
#endif

extern double**** LaplacianType2Stencil; 
extern double**** MassType2Stencil; 

PetscErrorCode Jacobian1ShellMatMult(Mat J, Vec in, Vec out) {
  PetscFunctionBegin;

  Jac1MFreeData* ctx;

  MatShellGetContext(J, (void**)(&ctx));

  assert(ctx != NULL);
  assert(ctx->da != NULL);  

  if(ctx->da->iAmActive()) {      
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

PetscErrorCode CreateJacobian1(ot::DA *da, Mat *jac)
{
  PetscInt  m,n;
  PetscFunctionBegin;
  //The size this processor owns ( without ghosts).
  m=n=da->getNodeSize();
  Jac1MFreeData* data = new Jac1MFreeData;
  data->da = da;
  data->Jmat_private = NULL;
  data->inTmp = NULL;
  data->outTmp = NULL;
  data->isFinestLevel = true;//single grid is always the finest.
  iC(MatCreateShell(da->getComm(), m ,n,PETSC_DETERMINE,PETSC_DETERMINE,
        (void*)(data), jac));
  iC(MatShellSetOperation(*jac ,MATOP_MULT, (void(*)(void)) Jacobian1MatMult));
  iC(MatShellSetOperation(*jac ,MATOP_GET_DIAGONAL, (void(*)(void)) Jacobian1MatGetDiagonal));
  iC(MatShellSetOperation(*jac ,MATOP_DESTROY, (void(*)(void)) Jacobian1MatDestroy));
  PetscFunctionReturn(0);
}

PetscErrorCode CreateActiveJacobian1(ot::DA *da, Mat *jac)
{
  PetscInt  m,n;
  PetscFunctionBegin;
  if(da->iAmActive()) {
    //The size this processor owns ( without ghosts).
    m=n=da->getNodeSize();
    Jac1MFreeData* data = new Jac1MFreeData;
    data->da = da;
    data->Jmat_private = NULL;
    data->inTmp = NULL;
    data->outTmp = NULL;
    data->isFinestLevel = true;//single grid is always the finest.
    iC(MatCreateShell(da->getCommActive(), m ,n,PETSC_DETERMINE,PETSC_DETERMINE,
          (void*)(data), jac));
    iC(MatShellSetOperation(*jac ,MATOP_MULT, (void(*)(void)) Jacobian1MatMult));
    iC(MatShellSetOperation(*jac ,MATOP_GET_DIAGONAL, (void(*)(void)) Jacobian1MatGetDiagonal));
    iC(MatShellSetOperation(*jac ,MATOP_DESTROY, (void(*)(void)) Jacobian1MatDestroy));
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ComputeJacobian1(ot::DA* da, Mat J)
{
  PetscFunctionBegin;
  //Do nothing
  PetscFunctionReturn(0);
}

PetscErrorCode CreateAndComputeMassMatrix(ot::DA* da, Mat* jac) {
  PetscInt  m,n;
  PetscFunctionBegin;
  //The size this processor owns ( without ghosts).
  m=n=da->getNodeSize();
  Jac1MFreeData* data = new Jac1MFreeData;
  data->da = da;
  data->Jmat_private = NULL;
  data->inTmp = NULL;
  data->outTmp = NULL;
  iC(MatCreateShell(da->getComm(), m ,n,PETSC_DETERMINE,PETSC_DETERMINE, (void*)(data), jac));
  iC(MatShellSetOperation(*jac ,MATOP_MULT, (void(*)(void)) MassMatMult));
  iC(MatShellSetOperation(*jac ,MATOP_GET_DIAGONAL, (void(*)(void)) MassMatGetDiagonal));
  iC(MatShellSetOperation(*jac ,MATOP_DESTROY, (void(*)(void)) MassMatDestroy));
  PetscFunctionReturn(0);
}

PetscErrorCode Jacobian1MatDestroy(Mat J) {
  PetscFunctionBegin;
  Jac1MFreeData *data;
  iC(MatShellGetContext( J, (void **)&data));
  if(data != NULL) {
    if(J != data->Jmat_private) {
      if(data->Jmat_private) {
        MatDestroy(data->Jmat_private);
        data->Jmat_private = NULL;
      }
    }
    if(data->inTmp) {
      VecDestroy(data->inTmp);
      data->inTmp = NULL;
    }
    if(data->outTmp) {
      VecDestroy(data->outTmp);
      data->outTmp = NULL;
    }
    delete data;
    data = NULL;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode MassMatDestroy(Mat J) {
  PetscFunctionBegin;    
  Jac1MFreeData *data;
  iC(MatShellGetContext( J, (void **)&data));
  if(data) {
    if(data->Jmat_private) {
      MatDestroy(data->Jmat_private);
      data->Jmat_private = NULL;
    }
    delete data;
    data = NULL;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode Jacobian1MatGetDiagonal(Mat J, Vec diag) {
  PetscFunctionBegin;

  PetscLogEventBegin(Jac1DiagEvent,diag,0,0,0);
  Jac1MFreeData *data;
  iC(MatShellGetContext(J, (void **)&data));
  if(data->isFinestLevel) {
    PetscLogEventBegin(Jac1FinestDiagEvent,diag,0,0,0);
  }
  iC(VecZeroEntries(diag));

  PetscScalar *diagArr;
  //Nodal,Non-Ghosted,Write,1 dof
  data->da->vecGetBuffer(diag,diagArr,false,false,false,1);

  //If Required Use Some Problem Specific information and set xFac, yFac, zFac.
  //Here, I'm simply assuming that the domain is of unit size.
  //To get the physical dimension of any element in a particular direction say
  //'x', xFac will be multiplied by (1u << (maxD-level)) of that element.
  unsigned int maxD;
  double hFac;

  PetscReal lapFac = 1.0;
  PetscReal massFac = 1.0;
  PetscTruth optFound;
  PetscOptionsGetReal("lap","-MatPropFac",&lapFac,&optFound);
  PetscOptionsGetReal("mass","-MatPropFac",&massFac,&optFound);

  if(data->da->iAmActive()) {
    maxD = (data->da->getMaxDepth());
    hFac = 1.0/((double)(1u << (maxD-1)));

    //Loop through All Elements including ghosted
    for(data->da->init<ot::DA_FLAGS::ALL>();
        data->da->curr() < data->da->end<ot::DA_FLAGS::ALL>(); 
        data->da->next<ot::DA_FLAGS::ALL>()) {
      //This returns the 8 vertices in the Morton order.
      unsigned int lev = data->da->getLevel(data->da->curr());
      double h = hFac*(1u << (maxD - lev));
      double fac1 = lapFac*h/2.0;
      double fac2 = massFac*h*h*h/8.0;
      unsigned int indices[8];
      data->da->getNodeIndices(indices); 
      unsigned char childNum = data->da->getChildNumber();
      unsigned char hnMask = data->da->getHangingNodeIndex(data->da->curr()); 
      unsigned char elemType = 0;
      GET_ETYPE_BLOCK(elemType,hnMask,childNum)
        for(int k = 0;k < 8;k++) {
          diagArr[indices[k]] +=  ((fac1*(LaplacianType2Stencil[childNum][elemType][k][k])) +  
              (fac2*(MassType2Stencil[childNum][elemType][k][k])));
        }//end k
    }//end i
  }

  data->da->vecRestoreBuffer(diag,diagArr,false,false,false,1);

  PetscLogFlops(44*(data->da->getGhostedElementSize()));
  PetscLogEventEnd(Jac1DiagEvent,diag,0,0,0);
  if(data->isFinestLevel) {
    PetscLogEventEnd(Jac1FinestDiagEvent,diag,0,0,0);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode MassMatGetDiagonal(Mat J, Vec diag) {
  PetscFunctionBegin;  
  Jac1MFreeData *data;
  iC(MatShellGetContext(J, (void **)&data));
  iC(VecZeroEntries(diag));

  //If Required Use Some Problem Specific information and set xFac, yFac, zFac.
  //Here, I'm simply assuming that the domain is of unit size.
  //To get the physical dimension of any element in a particular direction say
  //'x', xFac will be multiplied by (1u << (maxD-level)) of that element.
  unsigned int maxD;
  double hFac;
  if(data->da->iAmActive()) {
    maxD = (data->da->getMaxDepth());
    hFac = 1.0/((double)(1u << (maxD-1)));
  }

  PetscScalar *diagArr;
  //Nodal,Non-Ghosted,Write,1 dof
  data->da->vecGetBuffer(diag,diagArr,false,false,false,1);

  if(data->da->iAmActive()) {
    //Loop through All Elements including ghosted
    for(data->da->init<ot::DA_FLAGS::ALL>();
        data->da->curr() < data->da->end<ot::DA_FLAGS::ALL>();
        data->da->next<ot::DA_FLAGS::ALL>()) {
      //This returns the 8 vertices in the Morton order.
      unsigned int lev = data->da->getLevel(data->da->curr());
      double h = hFac*(1u << (maxD - lev));
      double fac2 = h*h*h/8.0;
      unsigned int indices[8];
      data->da->getNodeIndices(indices); 
      unsigned char childNum = data->da->getChildNumber();
      unsigned char hnMask = data->da->getHangingNodeIndex(data->da->curr()); 
      unsigned char elemType = 0;
      GET_ETYPE_BLOCK(elemType,hnMask,childNum)
        for(int k = 0;k < 8;k++) {
          diagArr[indices[k]] +=  (fac2*(MassType2Stencil[childNum][elemType][k][k]));
        }//end k        
    }//end i
  }

  data->da->vecRestoreBuffer(diag,diagArr,false,false,false,1);  
  PetscFunctionReturn(0);
}

#define JACOBIAN_TYPE1_MULT_BLOCK {\
  unsigned int lev = data->da->getLevel(data->da->curr());\
  double h = hFac*(1u << (maxD - lev));\
  double fac1 = lapFac*h/2.0;\
  double fac2 = massFac*h*h*h/8.0;\
  unsigned int indices[8];\
  data->da->getNodeIndices(indices);\
  unsigned char childNum = data->da->getChildNumber();\
  unsigned char hnMask = data->da->getHangingNodeIndex(data->da->curr());\
  unsigned char elemType = 0;\
  GET_ETYPE_BLOCK(elemType,hnMask,childNum)\
  for(int k = 0;k < 8;k++) {\
    for(int j=0;j<8;j++) {\
      outArr[indices[k]] +=  (((fac1*(LaplacianType2Stencil[childNum][elemType][k][j])) +\
            (fac2*(MassType2Stencil[childNum][elemType][k][j])))*inArr[indices[j]]);\
    } /*end for j*/\
  } /*end for k*/\
}

PetscErrorCode Jacobian1MatMult(Mat J, Vec in, Vec out)
{
  PetscFunctionBegin;

  Jac1MFreeData *data;
  PetscLogEventBegin(Jac1MultEvent,in,out,0,0);
  iC(MatShellGetContext( J, (void **)&data));  
  if(data->isFinestLevel) {
    PetscLogEventBegin(Jac1FinestMultEvent,in,out,0,0);
  }
  iC(VecZeroEntries(out));

  unsigned int maxD;
  double hFac;
  if(data->da->iAmActive()) {
    maxD = data->da->getMaxDepth(); 
    hFac = 1.0/((double)(1u << (maxD-1)));
  }

  PetscReal lapFac = 1.0;
  PetscReal massFac = 1.0;
  PetscTruth optFound;
  PetscOptionsGetReal("lap","-MatPropFac",&lapFac,&optFound);
  PetscOptionsGetReal("mass","-MatPropFac",&massFac,&optFound);

  PetscScalar *outArr=NULL;
  PetscScalar *inArr=NULL; 
  //Nodal,Non-Ghosted,Read,1 dof
  data->da->vecGetBuffer(in,inArr,false,false,true,1);

  //std::cout << __func__ << ": Updating Ghost..." << std::endl;
  if(data->da->iAmActive()) {
    data->da->ReadFromGhostsBegin<PetscScalar>(inArr,1);
  }

  //Nodal,Non-Ghosted,Write,1 dof
  data->da->vecGetBuffer(out,outArr,false,false,false,1);

  //Loop through All Independent Elements     

  if(data->da->iAmActive()) {
    for(data->da->init<ot::DA_FLAGS::INDEPENDENT>();
        data->da->curr() < data->da->end<ot::DA_FLAGS::INDEPENDENT>();
        data->da->next<ot::DA_FLAGS::INDEPENDENT>() ) {
      JACOBIAN_TYPE1_MULT_BLOCK
    }//end independent
  }

  if(data->da->iAmActive()) {
    data->da->ReadFromGhostsEnd<PetscScalar>(inArr);
  }
  //Loop through All Dependent Elements, i.e. Elements which have atleast one
  //vertex owned by this processor and at least one vertex not owned by this
  //processor. For simplicity, I may write to nodes that I do not own. These
  //can be discarded later. Every processor is only supposed to write to the
  //nodes it owns. An alternate strategy would require me to check if i own the
  //node before i write, I think this might be inefficient. So, if possible
  //allow me to write to nodes I do not own. You can discard them while
  //restoring the buffers.
  //Note: Whenever, I refer to nodes I refer to regular nodes only. If any of the
  //vertices are hanging then the regular node that will be used instead is used
  //to make the classifications such as ghost or non-ghost elem/node.

  if(data->da->iAmActive()) {
    for(data->da->init<ot::DA_FLAGS::DEPENDENT>();
        data->da->curr() < data->da->end<ot::DA_FLAGS::DEPENDENT>();
        data->da->next<ot::DA_FLAGS::DEPENDENT>() ) {
      JACOBIAN_TYPE1_MULT_BLOCK
    }//end loop for dependent elems
  }

  data->da->vecRestoreBuffer(in,inArr,false,false,true,1);
  data->da->vecRestoreBuffer(out,outArr,false,false,false,1);

  //Only useful flops, i.e. flops for the FEM calculation itself, are included.
  //The overhead due to looping through the octants and uncompressing the
  //nodelist and the transformations for the maps are not included. 
  PetscLogFlops(332*(data->da->getGhostedElementSize()));
  PetscLogEventEnd(Jac1MultEvent,in,out,0,0);
  if(data->isFinestLevel) {
    PetscLogEventEnd(Jac1FinestMultEvent,in,out,0,0);
  }

  PetscFunctionReturn(0);
}

#undef JACOBIAN_TYPE1_MULT_BLOCK

#define MASS_MULT_BLOCK {\
  unsigned int lev = data->da->getLevel(data->da->curr());\
  double h = hFac*(1u << (maxD - lev));\
  double fac2 = h*h*h/8.0;\
  unsigned int indices[8];\
  data->da->getNodeIndices(indices);\
  unsigned char childNum = data->da->getChildNumber();\
  unsigned char hnMask = data->da->getHangingNodeIndex(data->da->curr());\
  unsigned char elemType = 0;\
  GET_ETYPE_BLOCK(elemType,hnMask,childNum)\
  for(int k = 0;k < 8;k++) {\
    for(int j=0;j<8;j++) {\
      outArr[indices[k]] +=  ((fac2*(MassType2Stencil[childNum][elemType][k][j]))*inArr[indices[j]]);\
    } /*end for j*/\
  } /*end for k*/\
}

PetscErrorCode MassMatMult(Mat J, Vec in, Vec out)
{
  PetscFunctionBegin;
  Jac1MFreeData *data;

  iC(MatShellGetContext( J, (void **)&data));  

  iC (VecSet(out, 0));

  unsigned int maxD;
  double hFac;

  if(data->da->iAmActive()) {
    maxD = data->da->getMaxDepth(); 
    hFac = 1.0/((double)(1u << (maxD-1)));
  }

  PetscScalar *outArr=NULL;
  PetscScalar *inArr=NULL; 
  //Nodal,Non-Ghosted,Read,1 dof
  data->da->vecGetBuffer(in,inArr,false,false,true,1);

  if(data->da->iAmActive()) {
    data->da->ReadFromGhostsBegin<PetscScalar>(inArr,1);
  }

  //Nodal,Non-Ghosted,Write,1 dof
  data->da->vecGetBuffer(out,outArr,false,false,false,1);

  //Loop through All Independent Elements     
  if(data->da->iAmActive()) {
    for(data->da->init<ot::DA_FLAGS::INDEPENDENT>();
        data->da->curr() < data->da->end<ot::DA_FLAGS::INDEPENDENT>();
        data->da->next<ot::DA_FLAGS::INDEPENDENT>() ) {
      MASS_MULT_BLOCK
    }//end independent
  }

  // std::cout << "Finished Independent" << std::endl;
  if(data->da->iAmActive()) {
    data->da->ReadFromGhostsEnd<PetscScalar>(inArr);
  }
  //Loop through All Dependent Elements, i.e. Elements which have atleast one
  //vertex owned by this processor and at least one vertex not owned by this
  //processor. For simplicity, I may write to nodes that I do not own. These
  //can be discarded later. Every processor is only supposed to write to the
  //nodes it owns. An alternate strategy would require me to check if i own the
  //node before i write, I think this might be inefficient. So, if possible
  //allow me to write to nodes I do not own. You can discard them while
  //restoring the buffers.
  //Note: Whenever, I refer to nodes I refer to regular nodes only. If any of the
  //vertices are hanging then the regular node that will be used instead is used
  //to make the classifications such as ghost or non-ghost elem/node.

  if(data->da->iAmActive()) {
    for(data->da->init<ot::DA_FLAGS::DEPENDENT>();
        data->da->curr() < data->da->end<ot::DA_FLAGS::DEPENDENT>();
        data->da->next<ot::DA_FLAGS::DEPENDENT>() ) {
      MASS_MULT_BLOCK
    }//end loop for dependent elems
  }

  data->da->vecRestoreBuffer(in,inArr,false,false,true,1);

  data->da->vecRestoreBuffer(out,outArr,false,false,false,1);

  PetscFunctionReturn(0);
}

#undef MASS_MULT_BLOCK 

PetscErrorCode CreateAndComputeFullJacobian1(ot::DA* da,Mat * J)
{
  PetscFunctionBegin;

  assert(da->computedLocalToGlobal());

  char matType[30];
  PetscTruth typeFound;
  PetscOptionsGetString(PETSC_NULL,"-fullJacMatType",matType,30,&typeFound);
  if(!typeFound) {
    std::cout<<"I need a MatType for the full Jacobian matrix!"<<std::endl;
    MPI_Finalize();
    exit(0);		
  } 
  da->createMatrix(*J, matType, 1);
  MatZeroEntries(*J);
  std::vector<ot::MatRecord> records;

  unsigned int maxD;
  double hFac;
  if(da->iAmActive()) {
    maxD = da->getMaxDepth(); 
    hFac = 1.0/((double)(1u << (maxD-1)));
  }


  PetscReal lapFac = 1.0;
  PetscReal massFac = 1.0;
  PetscTruth optFound;
  PetscOptionsGetReal("lap","-MatPropFac",&lapFac,&optFound);
  PetscOptionsGetReal("mass","-MatPropFac",&massFac,&optFound);

  if(da->iAmActive()) {
    for(da->init<ot::DA_FLAGS::WRITABLE>();
        da->curr() < da->end<ot::DA_FLAGS::WRITABLE>();
        da->next<ot::DA_FLAGS::WRITABLE>()) {
      unsigned int lev = da->getLevel(da->curr());
      double h = hFac*(1u << (maxD - lev));
      double fac1 = lapFac*h/2.0;
      double fac2 = massFac*h*h*h/8.0;
      unsigned int indices[8];
      da->getNodeIndices(indices);
      unsigned char childNum = da->getChildNumber();
      unsigned char hnMask = da->getHangingNodeIndex(da->curr());
      unsigned char elemType = 0;
      GET_ETYPE_BLOCK(elemType,hnMask,childNum)
        for(int k = 0;k < 8;k++) {
          for(int j=0;j<8;j++) {
            ot::MatRecord currRec;
            currRec.rowIdx = indices[k];
            currRec.colIdx = indices[j];
            currRec.rowDim = 0;
            currRec.colDim = 0;
            currRec.val = ((fac1*(LaplacianType2Stencil[childNum][elemType][k][j])) +
                (fac2*(MassType2Stencil[childNum][elemType][k][j])));
            records.push_back(currRec);
          }//end for j
        }//end for k
      if(records.size() > 1000) {
        da->setValuesInMatrix(*J, records, 1, ADD_VALUES);
      }
    }//end writable
  }//end if active

  da->setValuesInMatrix(*J, records, 1, ADD_VALUES);

  MatAssemblyBegin(*J, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*J, MAT_FINAL_ASSEMBLY);

  PetscFunctionReturn(0);
}

PetscErrorCode CreateAndComputeFullActiveJacobian1(ot::DA* da,Mat * J)
{
  PetscFunctionBegin;

  assert(da->computedLocalToGlobal());

  if(da->iAmActive()) {
    char matType[30];
    PetscTruth typeFound;
    PetscOptionsGetString(PETSC_NULL,"-fullJacMatType",matType,30,&typeFound);
    if(!typeFound) {
      std::cout<<"I need a MatType for the full Jacobian matrix!"<<std::endl;
      MPI_Finalize();
      exit(0);		
    } 
    da->createActiveMatrix(*J, matType, 1);
    MatZeroEntries(*J);
    std::vector<ot::MatRecord> records;

    unsigned int maxD;
    double hFac;
    maxD = da->getMaxDepth(); 
    hFac = 1.0/((double)(1u << (maxD-1)));

    PetscReal lapFac = 1.0;
    PetscReal massFac = 1.0;
    PetscTruth optFound;
    PetscOptionsGetReal("lap","-MatPropFac",&lapFac,&optFound);
    PetscOptionsGetReal("mass","-MatPropFac",&massFac,&optFound);

    for(da->init<ot::DA_FLAGS::WRITABLE>();
        da->curr() < da->end<ot::DA_FLAGS::WRITABLE>();
        da->next<ot::DA_FLAGS::WRITABLE>()) {
      unsigned int lev = da->getLevel(da->curr());
      double h = hFac*(1u << (maxD - lev));
      double fac1 = lapFac*h/2.0;
      double fac2 = massFac*h*h*h/8.0;
      unsigned int indices[8];
      da->getNodeIndices(indices);
      unsigned char childNum = da->getChildNumber();
      unsigned char hnMask = da->getHangingNodeIndex(da->curr());
      unsigned char elemType = 0;
      GET_ETYPE_BLOCK(elemType,hnMask,childNum)
        for(int k = 0;k < 8;k++) {
          for(int j=0;j<8;j++) {
            ot::MatRecord currRec;
            currRec.rowIdx = indices[k];
            currRec.colIdx = indices[j];
            currRec.rowDim = 0;
            currRec.colDim = 0;
            currRec.val =
              ((fac1*(LaplacianType2Stencil[childNum][elemType][k][j])) +
               (fac2*(MassType2Stencil[childNum][elemType][k][j])));
            records.push_back(currRec);
          }//end for j
        }//end for k
      if(records.size() > 1000) {
        da->setValuesInMatrix(*J, records, 1, ADD_VALUES);
      }
    }//end writable

    da->setValuesInMatrix(*J, records, 1, ADD_VALUES);

    MatAssemblyBegin(*J, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*J, MAT_FINAL_ASSEMBLY);

  }//end if active

  PetscFunctionReturn(0);
}



