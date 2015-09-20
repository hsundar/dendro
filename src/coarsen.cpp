
/**
  @file Coarsen.C
  @brief A set of functions for coarsening octrees.
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  */

#include "TreeNode.h"
#include "parUtils.h"
#include "dendro.h"

#ifdef __DEBUG__
#ifndef __DEBUG_OCT__
#define __DEBUG_OCT__
#endif
#endif

namespace ot {

  int simpleCoarsen(const std::vector<TreeNode > &in, std::vector<TreeNode> &out, MPI_Comm comm) {
    PROF_SIMPLE_COARSE_BEGIN

      out.clear();

    std::vector<TreeNode> tmpIn = in;

    MPI_Comm   new_comm;
    par::splitComm2way(tmpIn.empty(), &new_comm, comm);

    if(!tmpIn.empty()) {
      for(int i = 0; i < tmpIn.size(); i++) {
        out.push_back(tmpIn[i].getParent());
      }//end for i

      std::vector<TreeNode> tmpOut;
      unsigned int dim = tmpIn[0].getDim();
      unsigned int maxDepth = tmpIn[0].getMaxDepth();

      //Although out is not empty on any processor in new_comm, we can't
      //guarantee that there will be sufficient number of elements in out after
      //removing duplicates across processors
      completeOctree(out, tmpOut,  dim, maxDepth, false, false, false, new_comm);

      out = tmpOut;
      tmpOut.clear();
    }

    PROF_SIMPLE_COARSE_END
  }//end function

  //New implementation. Written on April 21, 2008.
  //Assumption: in is sorted, linear, complete (not necessarily balanced)
  int coarsenOctree(const std::vector<TreeNode > &in, std::vector<TreeNode> &out) {
    PROF_COARSE_SEQ_BEGIN
      unsigned int inSz = static_cast<unsigned int>(in.size());
    out.clear();

    unsigned int dim; 
    if(inSz) {
      dim = in[0].getDim();
    }else {
      PROF_COARSE_SEQ_END
    }

    if (inSz < (1 << dim)) {
      out = in;
      std::cout<<"WARNING: Too Few Octants to Coarse."<<std::endl;	
      PROF_COARSE_SEQ_END
    }

    assert(in[0].getChildNumber() == 0);

    int prevFCidx = 0;

    while (prevFCidx < inSz) {
      int currFCidx = (prevFCidx + 1);
      //The Order of the operands is important here. Must Check for array out of
      //bounds first.
      while( (currFCidx < inSz) && (in[currFCidx].getChildNumber() != 0) ) {
        currFCidx++;
      }      
      if(currFCidx >= (prevFCidx + (1 << dim))) {
        //Coarsen the first 8 octants including prevFCidx and insert the rest
        //as they are. Only the first 8 octants will be coarsened
        out.push_back(in[prevFCidx].getParent());
        if( (prevFCidx + (1 << dim)) < currFCidx ) {
          out.insert(out.end(), 
              (in.begin() + (prevFCidx + (1 << dim))), (in.begin() + currFCidx));
        }
      } else {
        //Can't coarsen any octant in between prevFCidx and currFCidx
        if(prevFCidx < currFCidx) {
          out.insert(out.end(),
              (in.begin() + prevFCidx), (in.begin() + currFCidx));
        }
      }
      prevFCidx = currFCidx;
    }

    PROF_COARSE_SEQ_END
  }//end function

  //New Implementation. Written on April 21, 2008
  //Assumption: in is globally sorted, linear and complete (not necessarily
  //balanced).
  int coarsenOctree(const std::vector<TreeNode > &in, std::vector<TreeNode> &out,
      unsigned int dim, unsigned int maxDepth, MPI_Comm comm,
      bool skipPartition, MPI_Comm* newCommPtr, bool* iAmActive) {
    PROF_COARSE_BEGIN

      int size;
    MPI_Comm_size(comm, &size);
    out.clear();

    if(newCommPtr) {
      *newCommPtr = comm;
    }

    if(iAmActive) {
      *iAmActive = true;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    //Check if it is a single processor...
    if (size == 1) {

#ifndef __SILENT_MODE__
      std::cout<<"Coarsen Octree. inpSize: "<<in.size()<<" activeNpes: 1"<<std::endl; 
#endif

      coarsenOctree(in,out);
      PROF_COARSE_END
    }

    std::vector<TreeNode> tmpIn = in;
    DendroIntL inSz = tmpIn.size();

    //////////////////////////////////////////////////////////////////////////////////////////////////
    //Check if the input is too small...
    DendroIntL globSize;
    par::Mpi_Allreduce<DendroIntL>(&inSz, &globSize, 1, MPI_SUM, comm);

    int rank;
    MPI_Comm_rank(comm, &rank);

    //min grain size = 1000
    const DendroIntL THOUSAND = 1000;
    if (globSize < (THOUSAND*size)) {
      int splittingSize = (globSize/THOUSAND); 
      if(splittingSize == 0) {
        splittingSize = 1; 
      }

      unsigned int avgLoad = (globSize/splittingSize);
      int leftOvers = (globSize - (splittingSize*avgLoad));

      std::vector<ot::TreeNode> tmpIn2;
      if(rank >= splittingSize) {
        par::scatterValues<ot::TreeNode>(tmpIn, tmpIn2, 0, comm);
      }else if(rank < leftOvers) {
        par::scatterValues<ot::TreeNode>(tmpIn, tmpIn2, (avgLoad+1), comm);
      }else {
        par::scatterValues<ot::TreeNode>(tmpIn, tmpIn2, avgLoad, comm);
      }
      tmpIn.clear();

      MPI_Comm newComm;
      par::splitCommUsingSplittingRank(splittingSize, &newComm, comm);

#ifndef __SILENT_MODE__
      if(!rank) {
        std::cout<<"Input to Coarsen is small ("<<globSize
          <<"). npes = "<<size<<" Splitting Comm. "<<std::endl;
      }
#endif

      if(rank < splittingSize) {
        coarsenOctree(tmpIn2, out, dim, maxDepth, newComm, true, NULL, NULL);
      } else {
        if(iAmActive) {
          *iAmActive = false;
        }
      }
      tmpIn2.clear();

      if(newCommPtr) {
        *newCommPtr = newComm;
      }

      PROF_COARSE_END
    }//end the special case of too few elements.

#ifndef __SILENT_MODE__
    if(!rank) {
      std::cout<<"Coarsen Octree. inpSize: "<<globSize <<" activeNpes: "<<size<<std::endl; 
    }
#endif

    //////////////////////////////////////////////////////////////////////////////////////////////////
    //Normal case all processors have enough number of elements for overlap.
    //This partition will guarantee that each proc has the min grain size 
    if(!skipPartition) {
      par::partitionW<ot::TreeNode>(tmpIn, NULL, comm);
    }

    inSz = tmpIn.size();

    int firstFCidx = 0;

    while( (firstFCidx < inSz) && (tmpIn[firstFCidx].getChildNumber() != 0) ) {
      firstFCidx++;
    }

    if(firstFCidx == inSz) {
      firstFCidx = -1;
    }

    int lastIdxToSend = -1;

    int lastFCidx = -1;

    if(firstFCidx >= 0) {
      if(rank < (size-1)) {
        lastFCidx = (inSz - 1);
        while( (lastFCidx >= 0) && (tmpIn[lastFCidx].getChildNumber() != 0) ) {
          lastFCidx--;
        }

        if(lastFCidx >= 0) {
          lastIdxToSend = (inSz - lastFCidx);
        }
      } else {
        lastFCidx = inSz;
      }//end if PN
    }

    MPI_Request sendFirstIdxRequest;
    MPI_Request sendLastIdxRequest;
    MPI_Request recvFirstIdxRequest;
    MPI_Request recvLastIdxRequest;

    int lastFCidxOfPrevP;
    int firstFCidxOfNextP;

    if(rank) {
      //Recv lastFCidx from the prev processor
      par::Mpi_Irecv<int>(&lastFCidxOfPrevP, 1, (rank-1),
          1, comm, &recvLastIdxRequest);

      //Send firstFCidx to the prev processor
      par::Mpi_Issend<int>(&firstFCidx, 1, (rank-1),
          2, comm, &sendFirstIdxRequest);
    }

    if(rank < (size-1)) {
      //Recv firstFCidx from the next processor
      par::Mpi_Irecv<int>(&firstFCidxOfNextP, 1, (rank+1),
          2, comm, &recvFirstIdxRequest);

      //Send lastFCidx to the next processor
      par::Mpi_Issend<int>(&lastIdxToSend, 1, (rank+1),
          1, comm, &sendLastIdxRequest);
    }

    //Process local elements, overlapping with communication
    int prevFCidx = firstFCidx;

    //If firstFCidx == lastFCidx, it will be processed after recieving the
    //messages from the previous and next processors
    //On the last processor, we want to process all elements after the firstFC.
    //So, on the last processor lastFCidx is set to the end of the list (inSz). 
    //On all other processors, we must wait for the firstFCofNextP before
    //processing elements after (including) our lastFC. Hence lastFC is not
    //processed here for other processors.
    while (prevFCidx < lastFCidx) {
      int currFCidx = (prevFCidx + 1);
      while( (currFCidx < lastFCidx) &&
          (tmpIn[currFCidx].getChildNumber() != 0) ) {
        currFCidx++;
      }      
      if(currFCidx >= (prevFCidx + (1 << dim))) {
        //Coarsen the first 8 octants including prevFCidx and insert the rest
        //as they are. Only the first 8 octants will be coarsened
        out.push_back(tmpIn[prevFCidx].getParent());
        if((prevFCidx + (1 << dim)) < currFCidx) {
          out.insert(out.end(), 
              (tmpIn.begin() + (prevFCidx + (1 << dim))),
              (tmpIn.begin() + currFCidx));
        }
      } else {
        //Can't coarsen any octant in between prevFCidx and currFCidx
        if(prevFCidx < currFCidx) {
          out.insert(out.end(),
              (tmpIn.begin() + prevFCidx), (tmpIn.begin() + currFCidx));
        }
      }
      prevFCidx = currFCidx;
    }//end main while loop

    if(rank) {
      MPI_Status statusWait;
      MPI_Wait(&recvLastIdxRequest, &statusWait);

      if ( (lastFCidxOfPrevP >= 0) && (firstFCidx >= 0) ) {
        if ( (lastFCidxOfPrevP + firstFCidx) >= (1 << dim) ) {
          int stIdx = ((1 << dim) - lastFCidxOfPrevP);
          if(stIdx < 0) {
            stIdx = 0;
          }
          //The elements 0 to stIdx are siblings of an octant on the prev
          //processor. The prev. processor will add the parent
          if(stIdx < firstFCidx) {
            out.insert(out.begin(),
                (tmpIn.begin() + stIdx), (tmpIn.begin() + firstFCidx) );
          }
        } else {
          //No octant in the interval is coarsened 
          if(firstFCidx) {
            out.insert(out.begin(), tmpIn.begin(),
                (tmpIn.begin() + firstFCidx) );
          }
        }
      } else {
        //Special case
        if (lastFCidxOfPrevP < 0) {
          //There is no FC on the previous processor
          //Since, every processor has atleast 8 elements, immaterial of
          //whether or not our processor has any FCs we can add all the
          //elements until our firstFCidx (not inclusive). If we don't have any
          //FCs then we can add all our elements, since they will neither be
          //combined with the previous processor nor the next processor
          if (firstFCidx >= 0) {
            if (firstFCidx) {
              out.insert(out.begin(), tmpIn.begin(),
                  (tmpIn.begin() + firstFCidx));
            }
          } else {
            //Since there are no FCs on this processor, the main loop would
            //have done nothing
            assert(out.empty());
            out = tmpIn;
          }
        } else {
          //Since the previous processor has a FC and we have none and since
          //every processor has atleast 8 elements (in particular this
          //processor owns atleast 8 elements), there are more than 8 elements
          //between the FC on the prev. processor and the next FC (there may be
          //none. It is ok. That is to be treated just as if the next FC idx is
          //past the end of octree itself)
          //the last FC on the prev. processor and 7 successive elements will
          //be coarsened and the remaining elements until the next FC will be
          //added as it is. Since, the next FC (if it exists) is not on our
          //processor, we must simply add all our elements, except those that
          //are the siblings of the last FC on the prev. proc.
          assert(out.empty());
          if (lastFCidxOfPrevP < (1 << dim) ) {
            int stIdx = ((1 << dim) - lastFCidxOfPrevP); 
            out.insert(out.begin(), (tmpIn.begin() + stIdx), tmpIn.end());
          } else {
            out = tmpIn;
          }
        }
      }//end if special case
    }//end if not P0

    if (rank < (size-1)) {
      MPI_Status statusWait;
      MPI_Wait(&recvFirstIdxRequest, &statusWait);

      if ( (lastIdxToSend >= 0) && (firstFCidxOfNextP >= 0) ) {
        if ( (lastIdxToSend + firstFCidxOfNextP) >= (1 << dim) ) {
          out.push_back(tmpIn[lastFCidx].getParent());
          if ( (lastFCidx + (1 << dim)) < inSz ) {
            out.insert(out.end(),
                (tmpIn.begin() + (lastFCidx + (1 << dim))),
                (tmpIn.begin() + inSz));
          }
        } else {
          out.insert(out.end(), (tmpIn.begin() + lastFCidx),
              (tmpIn.begin() + inSz));
        }
      } else {
        //Special case
        if (lastIdxToSend < 0) {
          //This processor does not have any FCs. This case will be handled
          //while checking our firstFC with the lastFC of the prev P.
          //Only P0 will not handle this case, since it will not have any prev P.
          //However, P0 is guaranteed to have atleast one FC (the 0,0,0 node)
          //Nothing to be done here.
        } else {
          //This processor has a FC, but the next processor does not have any
          //FCs. Since the next processor has atleast 8 elements. There are
          //more than 8 elements between the lastFC on this processor and the
          //next FC in the linear octree.
          out.push_back(tmpIn[lastFCidx].getParent());          
          if ( (lastFCidx + (1 << dim)) < inSz ) {
            out.insert(out.end(),
                (tmpIn.begin() + (lastFCidx + (1 << dim))),
                (tmpIn.begin() + inSz));
          }
        }
      }//end if special case
    }//end if not PN

    if (rank) {
      MPI_Status statusWait;
      MPI_Wait(&sendFirstIdxRequest, &statusWait);
    }

    if (rank < (size-1)) {
      MPI_Status statusWait;
      MPI_Wait(&sendLastIdxRequest, &statusWait);
    }

    PROF_COARSE_END
  }//end function

}//end namespace


