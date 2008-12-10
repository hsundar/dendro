
/**
@file sys.h
@brief Initialization and Cleaning up of some private DENDRO variables.
@author Rahul Sampath, rahul.sampath@gmail.com
*/

#ifndef _SYS_H_
#define _SYS_H_

#include "petscsys.h"

namespace ot {

  /**
    @brief Call this function at the beginning of the main program, just after calling PetscInitialize()
    @return an error flag	
  */
  PetscErrorCode RegisterEvents();

}

#endif

