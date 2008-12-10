
#ifndef __UPDATE_CTX_H__
#define __UPDATE_CTX_H__

#include "mpi.h"
#include <vector>

namespace ot {

  class updateContext {
    public:
      void *                          buffer;
      void *                          keys;
      std::vector<MPI_Request*>       requests;

      bool operator== (updateContext other) {
        return( buffer == other.buffer ); 
      }

      ~updateContext() {
        requests.clear();
      }
  }; 

} //end namespace

#endif

