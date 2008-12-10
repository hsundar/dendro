
/**
  @file blockDiag.C
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  */

#include "petscmat.h"
#include "petscpc.h"
#include "blockDiag.h"
#include <cassert>

namespace ot {

  extern void (*getDofAndNodeSizeForPC_BlockDiag)(Mat pcMat,
      unsigned int & dof, unsigned int & nodeSize);

  extern void (*computeInvBlockDiagEntriesForPC_BlockDiag)(Mat pcMat,
      double **invBlockDiagEntries);


  PetscErrorCode PCSetUp_BlockDiag(PC pc) {

    PetscFunctionBegin;
    PC_BlockDiag* data = static_cast<PC_BlockDiag*>(pc->data);

    Mat pcMat = pc->pmat;
    PetscTruth isshell;
    PetscTypeCompare((PetscObject)pcMat, MATSHELL, &isshell);

    if(!isshell) {
      SETERRQ(PETSC_ERR_SUP," Expected a MATSHELL.");
      assert(false);
    }

    if(pc->setupcalled == 0) {
      if(getDofAndNodeSizeForPC_BlockDiag) {
        (*getDofAndNodeSizeForPC_BlockDiag)(pcMat, data->dof, data->nodeSize);
      } else {
        SETERRQ(PETSC_ERR_USER," Expected function to be set: getDofAndNodeSizeForPC_BlockDiag");
        assert(false);
      }

      //Allocate memory
      assert(data->invBlockDiagEntries == NULL);
      if((data->dof) && (data->nodeSize)) {
        typedef double* doublePtr;
        data->invBlockDiagEntries = new doublePtr[(data->dof)*(data->nodeSize)];
        for(int i = 0; i < ((data->dof)*(data->nodeSize)); i++) {
          data->invBlockDiagEntries[i] = new double[data->dof];
        }
        PetscLogObjectMemory(pc, (((data->dof)*(data->nodeSize))*sizeof(double)));
      }

      if(computeInvBlockDiagEntriesForPC_BlockDiag) {
        (*computeInvBlockDiagEntriesForPC_BlockDiag)(pcMat, data->invBlockDiagEntries);
      } else {
        SETERRQ(PETSC_ERR_USER,
            " Expected function to be set: computeInvBlockDiagEntriesForPC_BlockDiag");
        assert(false);
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode PCApply_BlockDiag(PC pc, Vec x, Vec y) {

    PetscFunctionBegin;
    PC_BlockDiag* data = static_cast<PC_BlockDiag*>(pc->data);

    if(!(data->invBlockDiagEntries)) {
      PCSetUp_BlockDiag(pc);
    }

    //y = invBlockDiagEntries*x
    VecZeroEntries(y);
    PetscScalar *yArr = NULL;
    PetscScalar *xArr = NULL;
    VecGetArray(y, &yArr);
    VecGetArray(x, &xArr);
    for(int i = 0; i < data->nodeSize; i++) {
      for(int j = 0; j < data->dof; j++) {
        for(int k = 0; k < data->dof; k++) {
          yArr[((data->dof)*i) + j] += 
            ((data->invBlockDiagEntries[((data->dof)*i) + j][k])*
             xArr[((data->dof)*i) + k]);
        }
      }
    }
    VecRestoreArray(y, &yArr);
    VecRestoreArray(x, &xArr);

    PetscFunctionReturn(0);
  }

  PetscErrorCode PCCreate_BlockDiag(PC pc) {

    PetscFunctionBegin;
    PC_BlockDiag* data = new PC_BlockDiag;

    pc->data = (void*)(data);

    PetscLogObjectMemory(pc, sizeof(PC_BlockDiag));

    //Initialize Data
    data->dof = 0;
    data->nodeSize = 0;
    data->invBlockDiagEntries = NULL;

    pc->ops->apply = PCApply_BlockDiag;
    pc->ops->setup = PCSetUp_BlockDiag;
    pc->ops->destroy = PCDestroy_BlockDiag;
    pc->ops->setfromoptions = PCSetFromOptions_BlockDiag;
    pc->ops->applytranspose = NULL;
    pc->ops->view = NULL;
    pc->ops->applyrichardson = NULL;
    pc->ops->applysymmetricleft = NULL;
    pc->ops->applysymmetricright = NULL;

    PetscFunctionReturn(0);
  }

  PetscErrorCode PCDestroy_BlockDiag(PC pc) {

    PetscFunctionBegin;
    PC_BlockDiag* data = static_cast<PC_BlockDiag*>(pc->data);

    if(data) {

      if(data->invBlockDiagEntries) {
        for(int i = 0; i < ((data->dof)*(data->nodeSize)); i++) {
          delete [] (data->invBlockDiagEntries[i]);
          data->invBlockDiagEntries[i] = NULL;
        }
        delete [] (data->invBlockDiagEntries);
        data->invBlockDiagEntries = NULL;
      }

      delete data;
      data = NULL;

    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode PCSetFromOptions_BlockDiag(PC pc) {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

} //end namespace 


