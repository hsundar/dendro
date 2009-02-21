
/**
  @file blockDiag.h
  @brief A Block Diagonal Precondioner. This can be used as an
  example for adding user-defined preconditioners. 
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  */

#ifndef _BLOCKDIAG_H_
#define _BLOCKDIAG_H_

#include "petscpc.h"
#include "petscvec.h"

namespace ot {

  /**
    @struct PC_BlockDiag	
    @brief Context for using the Block Diagonal preconditioner.
    @author Rahul Sampath

    This can be used as an example of adding your own user-defined preconditioners.
    */
  typedef struct {
    unsigned int dof;
    unsigned int nodeSize;
    double **invBlockDiagEntries;
  } PC_BlockDiag;

  PetscErrorCode PCSetUp_BlockDiag(PC pc);

  PetscErrorCode PCApply_BlockDiag(PC pc, Vec x, Vec y);

  PetscErrorCode PCDestroy_BlockDiag(PC pc);

  PetscErrorCode PCSetFromOptions_BlockDiag(PC pc);

  PetscErrorCode PCCreate_BlockDiag(PC pc);

} //end namespace 

#endif

