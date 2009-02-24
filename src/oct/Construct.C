
/**
  @file Construct.C
  @brief A set of functions for Octree Construction
  @author Hari Sundar, hsundar@gmail.com
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  */

#include "parUtils.h"
#include "TreeNode.h"
#include <cassert>
#include "nodeAndValues.h"
#include "binUtils.h"
#include "dendro.h"

#ifdef __DEBUG__
#ifndef __DEBUG_OCT__
#define __DEBUG_OCT__
#endif
#endif

namespace ot {

  int regularGrid2Octree(const std::vector<double>& elementValues,
      unsigned int N, unsigned int nx, unsigned int ny, unsigned int nz,
      unsigned int xs, unsigned int ys, unsigned int zs, std::vector<TreeNode>& linOct,
      unsigned int dim, unsigned int maxDepth, double thresholdFac, MPI_Comm comm) {
    PROF_RG2O_BEGIN

      int rank;
    int npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    bool isNvalid = binOp::isPowerOfTwo(N);

    assert(isNvalid);

    unsigned int rgLevel = binOp::fastLog2(N);
    unsigned int elemLen = (1u << (maxDepth - rgLevel));

    //1. Convert input to regular octree.
    std::vector<ot::NodeAndValues<double, 1> > tmpList(nx*ny*nz);

    for(int k = 0; k < nz; k++) {
      for(int j = 0; j < ny; j++) {
        for(int i = 0; i < nx; i++) {
          unsigned int currX = ((i + xs)*elemLen); 
          unsigned int currY = ((j + ys)*elemLen); 
          unsigned int currZ = ((k + zs)*elemLen);
          unsigned int idx = ((nx*(j + (ny*k))) + i);
          tmpList[idx].node = ot::TreeNode(currX, currY, currZ, rgLevel, dim, maxDepth);
          tmpList[idx].values[0] = elementValues[idx];
        }//end for i
      }//end for j
    }//end for k

    //2. Convert to Morton Ordering
    std::vector<ot::NodeAndValues<double, 1> > tnAndValsList;
    par::sampleSort<ot::NodeAndValues<double, 1> >(tmpList, tnAndValsList, comm);
    tmpList.clear();

    DendroIntL inSz = tnAndValsList.size();
    DendroIntL globInSize;
    par::Mpi_Allreduce<DendroIntL>(&inSz, &globInSize, 1, MPI_SUM, comm);

    assert( globInSize > (10*npes) );

    int npesCurr = npes;
    MPI_Comm commCurr = comm;

    bool repeatLoop = true;
    while(repeatLoop) {

      inSz = tnAndValsList.size();

      MPI_Request requests[4];

      //Send first 7 to previous processor
      //Send last 7 to next processor
      ot::NodeAndValues<double, 1> prevOcts[7];
      if( rank ) {
        par::Mpi_Irecv<ot::NodeAndValues<double, 1> >(prevOcts, 7, (rank - 1), 
            1, commCurr, requests);
      }

      ot::NodeAndValues<double, 1> nextOcts[7];
      if( rank < (npesCurr - 1) ) {
        par::Mpi_Irecv<ot::NodeAndValues<double, 1> >(nextOcts, 7, (rank + 1),
            1, commCurr, (requests + 1));
      }

      if( rank ) {
        par::Mpi_Issend<ot::NodeAndValues<double, 1> >((&(*(tnAndValsList.begin()))), 7,
            (rank - 1), 1, commCurr, (requests + 2));
      }

      if( rank < (npesCurr - 1) ) {
        par::Mpi_Issend<ot::NodeAndValues<double, 1> >((&(*(tnAndValsList.end() - 7))), 7,
            (rank + 1), 1, commCurr, (requests + 3));
      }

      MPI_Status statuses[4];

      if(npesCurr > 1) {
        MPI_Waitall(4, requests, statuses);
      }

      //Check and coarsen
      int idxOfPrevFC = -1;
      if( rank ) {
        for(int i = 0; i < 7; i++) {
          //There may be more than 1 FC, we only need the last one here
          if(prevOcts[i].node.getChildNumber() == 0) {
            idxOfPrevFC = i;
          }
        }//end for i
      }

      int idxOfNextFC = -1;
      if( rank < (npesCurr - 1) ) {
        for(int i = 0; i < 7; i++) {
          //There may be more than 1 FC, we only need the first one here
          if(nextOcts[i].node.getChildNumber() == 0) {
            idxOfNextFC = i;
            break;
          }
        }//end for i
      }

      int myFirstFC = -1;
      for(int i = 0; i < inSz; i++) {
        if(tnAndValsList[i].node.getChildNumber() == 0) {
          myFirstFC = i;
          break;
        }
      }//end for i

      int myLastFC = -1;
      if( myFirstFC >= 0 ) {
        for(int i = (inSz - 1); i >= 0; i--) {
          if(tnAndValsList[i].node.getChildNumber() == 0) {
            myLastFC = i;
            break;
          }
        }//end for i
      }

      //Process upto myFirstFC (exclusive)
      //Note, we check the threshold on both the processors so that the
      //partition is not changed unnecessarily. 
      if( (myFirstFC >= 0) && (idxOfPrevFC >= 0) ) {
        int fcGap = (myFirstFC + 7 - idxOfPrevFC);
        if( fcGap < 8 ) {
          //Can not coarsen
          if(myFirstFC) {
            tmpList.insert(tmpList.end(), tnAndValsList.begin(),
                tnAndValsList.begin() + myFirstFC);
          }
        } else if( fcGap > 8 ) {
          //Can coarsen upto (excluding) (1 + idxOfPrevFC) 
          double minVal = prevOcts[idxOfPrevFC].values[0];
          double maxVal = prevOcts[idxOfPrevFC].values[0];
          for(int i = idxOfPrevFC; i < 7; i++){
            double currVal = prevOcts[i].values[0];
            if(currVal < minVal) {
              minVal = currVal;
            }
            if(currVal > maxVal) {
              maxVal = currVal;
            }
          }//end for i
          for(int i = 0; i < (1 + idxOfPrevFC); i++){
            double currVal = tnAndValsList[i].values[0];
            if(currVal < minVal) {
              minVal = currVal;
            }
            if(currVal > maxVal) {
              maxVal = currVal;
            }
          }//end for i
          if((maxVal - minVal) >= thresholdFac) {
            //Can't coarsen
            if(myFirstFC) {
              tmpList.insert(tmpList.end(), tnAndValsList.begin(),
                  tnAndValsList.begin() + myFirstFC);
            }
          } else {
            //Can coarsen. The previous processor will add the parent
            if((1 + idxOfPrevFC) < myFirstFC) {
              tmpList.insert(tmpList.end(), tnAndValsList.begin() + (1 + idxOfPrevFC), 
                  tnAndValsList.begin() + myFirstFC);
            }
          }
        } else {
          //Can coarsen upto (excluding) myFirstFC
          double minVal = prevOcts[idxOfPrevFC].values[0];
          double maxVal = prevOcts[idxOfPrevFC].values[0];
          for(int i = idxOfPrevFC; i < 7; i++){
            double currVal = prevOcts[i].values[0];
            if(currVal < minVal) {
              minVal = currVal;
            }
            if(currVal > maxVal) {
              maxVal = currVal;
            }
          }//end for i
          for(int i = 0; i < myFirstFC; i++){
            double currVal = tnAndValsList[i].values[0];
            if(currVal < minVal) {
              minVal = currVal;
            }
            if(currVal > maxVal) {
              maxVal = currVal;
            }
          }//end for i 
          if((maxVal - minVal) >= thresholdFac) {
            //Can't coarsen
            if(myFirstFC) {
              tmpList.insert(tmpList.end(), tnAndValsList.begin(),
                  tnAndValsList.begin() + myFirstFC);
            }
          } else {
            //Can coarsen. The previous processor will add the parent 
          }
        }
      } else {
        if(myFirstFC >= 0) {
          //Can not coarsen
          if(myFirstFC) {
            tmpList.insert(tmpList.end(), tnAndValsList.begin(),
                tnAndValsList.begin() + myFirstFC);
          }
        } else if (idxOfPrevFC >= 0) {
          //Can coarsen upto (excluding) (1 + idxOfPrevFC) 
          double minVal = prevOcts[idxOfPrevFC].values[0];
          double maxVal = prevOcts[idxOfPrevFC].values[0];
          for(int i = idxOfPrevFC; i < 7; i++){
            double currVal = prevOcts[i].values[0];
            if(currVal < minVal) {
              minVal = currVal;
            }
            if(currVal > maxVal) {
              maxVal = currVal;
            }
          }//end for i
          for(int i = 0; i < myFirstFC; i++){
            double currVal = tnAndValsList[i].values[0];
            if(currVal < minVal) {
              minVal = currVal;
            }
            if(currVal > maxVal) {
              maxVal = currVal;
            }
          }//end for i
          if((maxVal - minVal) >= thresholdFac) {
            //Can't coarsen.
            tmpList = tnAndValsList;
          } else {
            //Can coarsen. The previous processor will add the parent 
            tmpList.insert(tmpList.end(), tnAndValsList.begin() + (1 + idxOfPrevFC), 
                tnAndValsList.end());
          }
        } else {
          //Can not coarsen
          tmpList = tnAndValsList;
        }
      }

      //Process myFirstFC (inclusive) to myLastFC (exclusive)
      int prevFCidx = myFirstFC;
      for(int idx = myFirstFC + 1; idx <= myLastFC; idx++) {
        if(tnAndValsList[idx].node.getChildNumber() == 0) {
          int fcGap = (idx - prevFCidx);
          if(fcGap < 8) {
            //Can't coarsen
            tmpList.insert(tmpList.end(), tnAndValsList.begin() + prevFCidx,
                tnAndValsList.begin() + idx);
          } else if(fcGap > 8) {
            //Can coarsen upto (excluding) (8 + prevFCidx)
            double minVal = tnAndValsList[prevFCidx].values[0];
            double maxVal = tnAndValsList[prevFCidx].values[0];
            double sumVal = tnAndValsList[prevFCidx].values[0];
            for(int j = prevFCidx; j < (8 + prevFCidx); j++) {
              double currVal = tnAndValsList[j].values[0];
              if(currVal < minVal) {
                minVal = currVal;
              }
              if(currVal > maxVal) {
                maxVal = currVal;
              }
              sumVal += currVal;
            }//end for j 
            if((maxVal - minVal) >= thresholdFac) {
              //Can't coarsen.
              tmpList.insert(tmpList.end(), tnAndValsList.begin() + prevFCidx,
                  tnAndValsList.begin() + idx);
            } else {
              //Can coarsen. 
              ot::NodeAndValues<double, 1> tmpObj;
              tmpObj.node = tnAndValsList[prevFCidx].node.getParent();
              tmpObj.values[0] = (sumVal/8.0);
              tmpList.push_back(tmpObj);
              tmpList.insert(tmpList.end(), tnAndValsList.begin() + prevFCidx + 8,
                  tnAndValsList.begin() + idx);
            }
          } else {
            //Can coarsen upto (excluding) idx
            double minVal = tnAndValsList[prevFCidx].values[0];
            double maxVal = tnAndValsList[prevFCidx].values[0];
            double sumVal = tnAndValsList[prevFCidx].values[0];
            for(int j = prevFCidx; j < idx; j++) {
              double currVal = tnAndValsList[j].values[0];
              if(currVal < minVal) {
                minVal = currVal;
              }
              if(currVal > maxVal) {
                maxVal = currVal;
              }
              sumVal += currVal;
            }//end for j 
            if((maxVal - minVal) >= thresholdFac) {
              //Can't coarsen.
              tmpList.insert(tmpList.end(), tnAndValsList.begin() + prevFCidx,
                  tnAndValsList.begin() + idx);
            } else {
              //Can coarsen. 
              ot::NodeAndValues<double, 1> tmpObj;
              tmpObj.node = tnAndValsList[prevFCidx].node.getParent();
              tmpObj.values[0] = (sumVal/8.0);
              tmpList.push_back(tmpObj);
            }
          }
          prevFCidx = idx;
        }//end if found FC
      }//end for idx

      //Process myLastFC (inclusive) to end
      if( (myLastFC >= 0) && (idxOfNextFC >= 0) ) {
        int fcGap = inSz + idxOfNextFC - myLastFC;
        if(fcGap < 8) {
          //Can't coarsen
          tmpList.insert(tmpList.end(), tnAndValsList.begin() + myLastFC,
              tnAndValsList.end());
        } else if (fcGap > 8) {
          //Can coarsen upto (excluding) (myLastFC + 8)
          double minVal = tnAndValsList[myLastFC].values[0];
          double maxVal = tnAndValsList[myLastFC].values[0];
          double sumVal = tnAndValsList[myLastFC].values[0];
          if((myLastFC + 8) > inSz) {
            for(int i = myLastFC; i < inSz; i++){
              double currVal = tnAndValsList[i].values[0];
              if(currVal < minVal) {
                minVal = currVal;
              }
              if(currVal > maxVal) {
                maxVal = currVal;
              }
              sumVal += currVal;
            }//end for i
            for(int i = 0; i < (myLastFC + 8 - inSz); i++){
              double currVal = nextOcts[i].values[0];
              if(currVal < minVal) {
                minVal = currVal;
              }
              if(currVal > maxVal) {
                maxVal = currVal;
              }
              sumVal += currVal;
            }//end for i         
          } else {
            for(int i = myLastFC; i < (myLastFC + 8); i++){
              double currVal = tnAndValsList[i].values[0];
              if(currVal < minVal) {
                minVal = currVal;
              }
              if(currVal > maxVal) {
                maxVal = currVal;
              }
              sumVal += currVal;
            }//end for i
          }
          if((maxVal - minVal) >= thresholdFac) {
            //Can't coarsen
            tmpList.insert(tmpList.end(), tnAndValsList.begin() + myLastFC,
                tnAndValsList.end());
          } else {
            //Can coarsen
            ot::NodeAndValues<double, 1> tmpObj;
            tmpObj.node = tnAndValsList[myLastFC].node.getParent();
            tmpObj.values[0] = (sumVal/8.0);
            tmpList.push_back(tmpObj);
            if((myLastFC + 8) < inSz) {
              tmpList.insert(tmpList.end(), tnAndValsList.begin() + myLastFC + 8,
                  tnAndValsList.end());
            }
          }
        } else {
          //Can coarsen  
          double minVal = tnAndValsList[myLastFC].values[0];
          double maxVal = tnAndValsList[myLastFC].values[0];
          double sumVal = tnAndValsList[myLastFC].values[0];
          if((myLastFC + 8) > inSz) {
            for(int i = myLastFC; i < inSz; i++){
              double currVal = tnAndValsList[i].values[0];
              if(currVal < minVal) {
                minVal = currVal;
              }
              if(currVal > maxVal) {
                maxVal = currVal;
              }
              sumVal += currVal;
            }//end for i
            for(int i = 0; i < (myLastFC + 8 - inSz); i++){
              double currVal = nextOcts[i].values[0];
              if(currVal < minVal) {
                minVal = currVal;
              }
              if(currVal > maxVal) {
                maxVal = currVal;
              }
              sumVal += currVal;
            }//end for i         
          } else {
            for(int i = myLastFC; i < (myLastFC + 8); i++){
              double currVal = tnAndValsList[i].values[0];
              if(currVal < minVal) {
                minVal = currVal;
              }
              if(currVal > maxVal) {
                maxVal = currVal;
              }
              sumVal += currVal;
            }//end for i
          }
          if((maxVal - minVal) >= thresholdFac) {
            //Can't coarsen
            tmpList.insert(tmpList.end(), tnAndValsList.begin() + myLastFC,
                tnAndValsList.end());
          } else {
            //Can coarsen
            ot::NodeAndValues<double, 1> tmpObj;
            tmpObj.node = tnAndValsList[myLastFC].node.getParent();
            tmpObj.values[0] = (sumVal/8.0);
            tmpList.push_back(tmpObj);
          }
        }
      } else {
        if(myLastFC >= 0) {
          //Can coarsen upto (excluding) (myLastFC + 8)
          double minVal = tnAndValsList[myLastFC].values[0];
          double maxVal = tnAndValsList[myLastFC].values[0];
          double sumVal = tnAndValsList[myLastFC].values[0];
          if((myLastFC + 8) > inSz) {
            for(int i = myLastFC; i < inSz; i++){
              double currVal = tnAndValsList[i].values[0];
              if(currVal < minVal) {
                minVal = currVal;
              }
              if(currVal > maxVal) {
                maxVal = currVal;
              }
              sumVal += currVal;
            }//end for i
            if(rank < (npesCurr - 1)) {
              for(int i = 0; i < (myLastFC + 8 - inSz); i++){
                double currVal = nextOcts[i].values[0];
                if(currVal < minVal) {
                  minVal = currVal;
                }
                if(currVal > maxVal) {
                  maxVal = currVal;
                }
                sumVal += currVal;
              }//end for i         
            }
          } else {
            for(int i = myLastFC; i < (myLastFC + 8); i++){
              double currVal = tnAndValsList[i].values[0];
              if(currVal < minVal) {
                minVal = currVal;
              }
              if(currVal > maxVal) {
                maxVal = currVal;
              }
              sumVal += currVal;
            }//end for i
          } 
          if((maxVal - minVal) >= thresholdFac) {
            //Can't coarsen
            tmpList.insert(tmpList.end(), tnAndValsList.begin() + myLastFC,
                tnAndValsList.end());
          } else {
            //Can coarsen
            ot::NodeAndValues<double, 1> tmpObj;
            tmpObj.node = tnAndValsList[myLastFC].node.getParent();
            tmpObj.values[0] = (sumVal/8.0);
            tmpList.push_back(tmpObj);
            if((myLastFC + 8) < inSz) {
              tmpList.insert(tmpList.end(), tnAndValsList.begin() + myLastFC + 8,
                  tnAndValsList.end());
            }
          }
        } else {
          //Can't coarsen.
          //This case is the same as myFirstFC < 0 and it has been taken care
          //of earlier. 
        }
      }

      DendroIntL outSz = tmpList.size();
      DendroIntL globOutSize;
      par::Mpi_Allreduce<DendroIntL>(&outSz, &globOutSize, 1, MPI_SUM, commCurr);

      if(globOutSize == globInSize) {
        repeatLoop = false;
      }

      //prepare for next iteration...
      globInSize = globOutSize;
      tnAndValsList = tmpList;
      tmpList.clear();

      if(globInSize < (10*npesCurr)) {
        int splittingSize = (globInSize/10); 
        if(splittingSize == 0) {
          splittingSize = 1; 
        }

        unsigned int avgLoad = (globInSize/splittingSize);
        int leftOvers = (globInSize - (splittingSize*avgLoad));

        if(rank >= splittingSize) {
          par::scatterValues<ot::NodeAndValues<double, 1> >(tnAndValsList, tmpList, 0, commCurr);
        }else if(rank < leftOvers) {
          par::scatterValues<ot::NodeAndValues<double, 1> >(tnAndValsList, tmpList, (avgLoad + 1), commCurr);
        }else {
          par::scatterValues<ot::NodeAndValues<double, 1> >(tnAndValsList, tmpList, avgLoad, commCurr);
        }
        tnAndValsList = tmpList;
        tmpList.clear();

        MPI_Comm newComm;
        par::splitCommUsingSplittingRank(splittingSize, &newComm, commCurr);
        commCurr = newComm;
        npesCurr = splittingSize;

        if(rank >= splittingSize) {
          repeatLoop = false;
        }
      }//end if splitting comm

    }//end while

    linOct.clear();
    for(int i = 0; i < tnAndValsList.size(); i++) {
      linOct.push_back(tnAndValsList[i].node);
    }//end for i

    PROF_RG2O_END
  }//end function

  //New Implementation. Written on April 20, 2008
  int completeOctree(const std::vector<TreeNode > & inp,
      std::vector<TreeNode > & out, 
      unsigned int dim, unsigned int maxDepth, bool isUnique,
      bool isSorted, bool assertNoEmptyProcs, MPI_Comm comm) {

    PROF_N2O_BEGIN

      TreeNode root (dim, maxDepth);

    int size;
    MPI_Comm_size(comm, &size);

    if(size == 1) {
      completeSubtree(root, inp, out, dim, maxDepth, isUnique, isSorted);
      PROF_N2O_END
    }//end single proc case	

    out = inp;

    //Sort and remove duplicate leaves.
    if (!isUnique) {
      par::removeDuplicates<ot::TreeNode>(out, isSorted, comm);
    } else if (!isSorted) {
      std::vector<TreeNode> tmpOut;
      par::sampleSort<ot::TreeNode>(out, tmpOut,comm); 
      out = tmpOut;
      tmpOut.clear();
    }

    //Remove empty processors...    
    MPI_Comm   new_comm;
    if(assertNoEmptyProcs) {
      new_comm = comm;
      assert(!out.empty());
    } else {
      par::splitComm2way(out.empty(), &new_comm, comm);
    }

    if(!(out.empty())) {	
      int new_rank, new_size; 

      MPI_Comm_rank (new_comm, &new_rank);
      MPI_Comm_size (new_comm, &new_size);

      MPI_Request requestSend;    
      MPI_Request requestRecv;

      ot::TreeNode begBuf;
      ot::TreeNode lastElem;

      if(new_rank) {
        //Recv
        par::Mpi_Irecv<ot::TreeNode>(&begBuf, 1, (new_rank-1) , 1, new_comm, &requestRecv);
      }//end if not P0

      if (new_rank < (new_size-1)) {
        lastElem = out[out.size() - 1];
        //Send
        par::Mpi_Issend<ot::TreeNode>(&lastElem, 1, (new_rank+1), 1, new_comm, &requestSend);
      }//end if not PN 


      //Add missing corners to complete the region.
      //Add the first corner leaf on the first processor.
      if (new_rank == 0) {
        ot::TreeNode minCorner(0,0,0,maxDepth,dim,maxDepth);
#ifdef __DEBUG_OCT__
        assert(areComparable(out[0], minCorner));
#endif
        if ( (out[0] != minCorner) && (!out[0].isAncestor(minCorner)) ) {
          ot::TreeNode ncaTmp = getNCA(out[0],minCorner);
          std::vector<ot::TreeNode>  kids;
          ncaTmp.addChildren(kids);
          out.insert(out.begin(),kids[0]);
          kids.clear();
        }//end if
      }//end if

      //Add the last corner leaf on the last processor.
      if (new_rank == (new_size-1) ) {
        ot::TreeNode maxCorner(((1u << maxDepth)-1),
            ((1u << maxDepth)-1),((1u << maxDepth)-1),maxDepth,dim,maxDepth);
#ifdef __DEBUG_OCT__
        assert(areComparable(out[out.size()-1], maxCorner));
#endif
        if ( (out[out.size()-1] != maxCorner) && 
            (!out[out.size()-1].isAncestor(maxCorner)) ) {
          ot::TreeNode ncaTmp =  getNCA(out[out.size()-1],maxCorner);
          std::vector<ot::TreeNode> kids;
          ncaTmp.addChildren(kids);
          out.insert(out.end(),kids[(1 << dim)-1]);
          kids.clear();
        }//end if
      }//end if

      std::vector<ot::TreeNode> tmpList;
      for(unsigned int i = 0; i < (out.size()-1); i++) {
#ifdef __DEBUG_OCT__
        assert(areComparable(out[i], out[i+1]));
#endif
        if(out[i].isAncestor(out[i+1])) {
          appendCompleteRegion(out[i], out[i+1], tmpList, false, false);
        } else {
          appendCompleteRegion(out[i], out[i+1], tmpList, true, false);
        }
      }//end for i

      //Only the last processor adds the last element. All the other processors would have 
      //sent it to the next processor, which will add it if it is not an ancestor of 
      //the first element on that processor
      if(new_rank == (new_size - 1) ) {
        tmpList.push_back(out[out.size() - 1]);
      }

      if(new_rank) {
        MPI_Status statusWait;
        MPI_Wait(&requestRecv, &statusWait);

        std::vector<ot::TreeNode> begList;
#ifdef __DEBUG_OCT__
        assert(areComparable(begBuf, out[0]));
#endif
        if(begBuf.isAncestor(out[0])) {
          appendCompleteRegion(begBuf, out[0], begList, false, false);
        } else {
          appendCompleteRegion(begBuf, out[0], begList, true, false);
        }
        out = begList;
        begList.clear();
        out.insert(out.end(), tmpList.begin(), tmpList.end());
      } else {
        out = tmpList;
      }

      tmpList.clear();

      if(new_rank < (new_size-1)) {
        MPI_Status statusWait;
        MPI_Wait(&requestSend, &statusWait);
      }

    } //out not empty


    PROF_N2O_END

  }//end function

  //New Implementation. Written on April 19, 2008
  int completeSubtree(TreeNode block, const std::vector<TreeNode > & inp, std::vector<TreeNode > & out,
      unsigned int dim, unsigned int maxDepth, bool isUnique, bool isSorted) {

    PROF_N2O_SEQ_BEGIN

      out = inp;
    //Sort and remove duplicate leaves.
    if (!isUnique) {
      seq::makeVectorUnique<TreeNode>(out, isSorted) ;
    } else if (!isSorted) {
      sort(out.begin(),out.end());
    }

    if (!out.empty()) {
      //Add missing corners to complete the region.
      //Add the first corner .
      ot::TreeNode minCorner = block.getDFD();
#ifdef __DEBUG_OCT__
      assert(areComparable(out[0], minCorner));
#endif
      if ( (out[0] != minCorner) && (!out[0].isAncestor(minCorner)) ) {
        ot::TreeNode ncaTmp = getNCA(out[0],minCorner);
        std::vector<ot::TreeNode> kids;
        ncaTmp.addChildren(kids);
        out.insert(out.begin(),kids[0]);
        kids.clear();
      }//end if

      //Add the last corner.
      ot::TreeNode maxCorner = block.getDLD();
#ifdef __DEBUG_OCT__
      assert(areComparable(out[out.size()-1], maxCorner));
#endif
      if ( (out[out.size()-1] != maxCorner) && (!out[out.size()-1].isAncestor(maxCorner)) ) {
        ot::TreeNode ncaTmp = getNCA(out[out.size()-1],maxCorner);
        std::vector<ot::TreeNode> kids;
        ncaTmp.addChildren(kids);
        out.insert(out.end(),kids[(1 << dim)-1]);
        kids.clear();
      }//end if

      std::vector<ot::TreeNode> tmpList;
      for(unsigned int i = 0; i < (out.size()-1); i++) {
#ifdef __DEBUG_OCT__
        assert(areComparable(out[i], out[i+1]));
#endif
        if(out[i].isAncestor(out[i+1])) {
          appendCompleteRegion(out[i], out[i+1], tmpList, false, false);
        } else {
          appendCompleteRegion(out[i], out[i+1], tmpList, true, false);
        }
      }//end for i

      tmpList.push_back(out[out.size() - 1]);

      out = tmpList;
      tmpList.clear();
    }//end if out empty

    PROF_N2O_SEQ_END

  }//end function

  int points2Octree(std::vector<double>& pts, double * gLens, std::vector<TreeNode> & nodes,
      unsigned int dim, unsigned int maxDepth, unsigned int maxNumPts, MPI_Comm comm ) {

    PROF_P2O_BEGIN

      int size;
    MPI_Comm_size(comm, &size);

    //Sequential...
    if (size == 1) {
      points2OctreeSeq(pts, gLens, nodes, dim, maxDepth, maxNumPts); 
      PROF_P2O_END
    }//end if sequential

    int rank;
    MPI_Comm_rank(comm, &rank);

    TreeNode root (dim, maxDepth);

    if (maxDepth == 0) {
      if (!rank) {
        nodes.resize(1);
        nodes[0] = root;
      } else {
        nodes.resize(0);
      }
      PROF_P2O_END
    }

    unsigned int ptsLen = pts.size();
    assert( (dim == 1) || (dim == 2) || (dim == 3 ));
    assert( (ptsLen%dim) == 0);

    DendroIntL numNodes = (ptsLen/dim);

    DendroIntL totSize;
    par::Mpi_Allreduce<DendroIntL>(&numNodes, &totSize, 1, MPI_SUM, comm);

    if (totSize <= static_cast<DendroIntL>(maxNumPts)) {
      if (!rank) {
        nodes.resize(1);
        nodes[0] = root;
      } else {
        nodes.resize(0);
      }//end if-else
      PROF_P2O_END
    }//end if

    //Tackle small problems separately....
    //min Grain size = 1000 
    const DendroIntL THOUSAND = 1000;
    if (totSize < (THOUSAND*size)) {
      int splittingSize = (totSize/THOUSAND); 

      if(splittingSize == 0) {
        splittingSize = 1; 
      }

      unsigned int avgLoad = (totSize/splittingSize);
      int leftOvers = (totSize - (splittingSize*avgLoad));

      std::vector<double> tmpPts;
      if(rank >= splittingSize) {
        par::scatterValues<double>(pts, tmpPts, 0, comm);
      }else if(rank < leftOvers) {
        par::scatterValues<double>(pts, tmpPts, (dim*(avgLoad+1)), comm);
      }else {
        par::scatterValues<double>(pts, tmpPts, (dim*avgLoad), comm);
      }
      pts.clear();

      MPI_Comm newComm;
      par::splitCommUsingSplittingRank(splittingSize, &newComm, comm);

#ifndef __SILENT_MODE__
      if(!rank) {
        std::cout<<"Input to p2o is small ("<<totSize
          <<"). npes = "<<size<<" Splitting Comm. "<<std::endl;
      }
#endif

      if(rank < splittingSize) {
        points2Octree(tmpPts, gLens, nodes, dim, maxDepth, maxNumPts, newComm); 
      }
      tmpPts.clear();

      PROF_P2O_END
    }//reduce procs for small problems

    //Tackle large problems.... 

    nodes.resize(numNodes);

    for (DendroIntL i = 0; i < numNodes; i++) {
      //The constructor will ignore unnecessary arguments (for lower
      //dimensions).
      unsigned int px = (unsigned int)(pts[i*dim]*((double)(1u << maxDepth))/(gLens[0]));
      unsigned int py, pz = 0;
      if(dim > 1) {
        py = (unsigned int)(pts[(i*dim)+1]*((double)(1u << maxDepth))/gLens[1]);
        if(dim > 2) {
          pz = (unsigned int)(pts[(i*dim)+2]*((double)(1u << maxDepth))/gLens[2]);
        }
      }
      nodes[i] = TreeNode(px, py, pz, maxDepth, dim, maxDepth);
    }//end for

    pts.clear(); 

    std::vector<ot::TreeNode> tmpNodes ;

    //Sort nodes (pts.) and partition them.
    par::sampleSort<ot::TreeNode>(nodes, tmpNodes, comm); 
    nodes = tmpNodes;
    tmpNodes.clear();

    std::vector<ot::TreeNode> leaves;
    std::vector<ot::TreeNode> minsAllBlocks;

    blockPartStage1_p2o(nodes, leaves, dim, maxDepth, comm);
    blockPartStage2_p2o(nodes, leaves, minsAllBlocks, dim, maxDepth, comm);
    //leaves will be sorted.

    p2oLocal(nodes, leaves, maxNumPts, dim, maxDepth);

    PROF_P2O_END

  }//end function

  //Added on April 19, 2008
  int points2OctreeSeq(std::vector<double>& pts, double * gLens, std::vector<TreeNode> & nodes,
      unsigned int dim, unsigned int maxDepth, unsigned int maxNumPts) {

    PROF_P2O_SEQ_BEGIN

      if (maxDepth == 0) {    
        nodes.resize(1);
        nodes[0] = TreeNode(dim,maxDepth);
        PROF_P2O_SEQ_END
      }

    unsigned int ptsLen = pts.size();
    assert( (dim == 1) || (dim == 2) || (dim == 3 ));
    assert( (ptsLen%dim) == 0);

    TreeNode root (dim,maxDepth);
    unsigned int numNodes = ptsLen/dim;
    nodes.resize(numNodes);

    for (unsigned int i = 0; i < numNodes; i++) {
      //The constructor will ignore unnecessary arguments (for lower
      //dimensions).
      unsigned int px = (unsigned int)(pts[i*dim]*((double)(1u << maxDepth))/(gLens[0]));
      unsigned int py, pz = 0;
      if(dim > 1) {
        py = (unsigned int)(pts[(i*dim)+1]*((double)(1u << maxDepth))/gLens[1]);
        if(dim > 2) {
          pz = (unsigned int)(pts[(i*dim)+2]*((double)(1u << maxDepth))/gLens[2]);
        }
      }
      nodes[i] = TreeNode(px, py, pz, maxDepth, dim, maxDepth);
    }//end for

    pts.clear(); 

    unsigned int totSize = nodes.size();

    if (totSize <= maxNumPts) {
      nodes.resize(1);
      nodes[0] = root;
      PROF_P2O_SEQ_END
    }//end if

    //Sort nodes (pts.) and partition them.  
    sort(nodes.begin(), nodes.end());

    std::vector<TreeNode> leaves;
    leaves.push_back(root);

    p2oLocal(nodes, leaves, maxNumPts, dim, maxDepth);

    PROF_P2O_SEQ_END

  }//end function

  int p2oLocal(std::vector<TreeNode> & nodes, std::vector<TreeNode>& leaves,
      unsigned int maxNumPts, unsigned int dim, unsigned int maxDepth) {
    PROF_P2O_LOCAL_BEGIN

      std::vector<TreeNode> nodesDup = nodes;
    nodes.clear();

    std::vector<TreeNode> addList;
    std::vector<TreeNode> nodesTmp;

    std::vector<TreeNode> *nodeSrc = &nodes;
    std::vector<TreeNode> *nodeDest = &nodesTmp;
    std::vector<TreeNode> *leafSrc = &leaves;
    std::vector<TreeNode> *leafDest = &addList;

    //The source and destination alternate for each iteration.
    do {
      //Leaves remain sorted inside this loop!
      //Reset all wts to 0.
      for (unsigned int i = 0; i < leafSrc->size(); i++) {
        (*leafSrc)[i].setWeight(0);
      }

      //nodesDup and leaves are both sorted at this point.
      unsigned int nextNode = 0;
      unsigned int nextPt = 0;
      //Some elements of nodesDup are not inside any element in leaves.
      //Note, this is different from other similar loops.
      //Here, nodesDup is merely  a counter and remains unchanged for all
      //outer do-while iterations, if a box had less than maxNumpts, it will
      //be removed from leaves and placed in nodeDest. However, all its
      //decendants in nodesDup still remain.
      //Also, note the first pt. need not be in any block. Hence, there is an
      //extra if-else wrapper to skip pts, that lie outside all blocks and
      //are lesser than them in the Morton order.
      while (nextPt < nodesDup.size()) {
        if (nodesDup[nextPt] >= (*leafSrc)[nextNode]) {
#ifdef __DEBUG_OCT__
          assert(areComparable((*leafSrc)[nextNode], nodesDup[nextPt]));
#endif
          if (((*leafSrc)[nextNode].isAncestor(nodesDup[nextPt])) ||
              ((*leafSrc)[nextNode] == nodesDup[nextPt])) {
            (*leafSrc)[nextNode].addWeight(1);
            nextPt++;
          } else {
            nextNode++;
            if (nextNode == leafSrc->size()) {
              //Note: No assert(false) here. 
              break;
            }
          }
        } else {
          nextPt++;
        }
      }//end while

      leafDest->resize((1 << dim)*(leafSrc->size()));
      unsigned int addListSz=0;     
      unsigned int nodesPrevSz = nodeSrc->size();
      nodeSrc->resize(nodesPrevSz + (leafSrc->size()));
      nodeDest->resize(nodeSrc->size());
      unsigned int tmpCtr=0;
      unsigned int nodesCtr=0;
      // int cntX=0;
      for (unsigned int i = 0; i < leafSrc->size(); i++) {
        if ( ((*leafSrc)[i].getWeight() > maxNumPts) &&
            ((*leafSrc)[i].getLevel() < maxDepth) ) {
          TreeNode thisNode = (*leafSrc)[i];
          std::vector<TreeNode> children;
          thisNode.addChildren(children);
          for (int ci = 0; ci < (1 << dim); ci++) {
            (*leafDest)[addListSz] = children[ci];
            addListSz++;
          }
          children.clear();
        } else {
          while ((nodesCtr < nodesPrevSz) &&
              ((*nodeSrc)[nodesCtr] < (*leafSrc)[i]) ) {
            (*nodeDest)[tmpCtr++] = (*nodeSrc)[nodesCtr++];
          }
          (*nodeDest)[tmpCtr++] = (*leafSrc)[i];
        }//end if-else
      }//end for i

      for (unsigned int i = nodesCtr; i < nodesPrevSz; i++) {
        (*nodeDest)[tmpCtr++] = (*nodeSrc)[i];
      }

      leafDest->resize(addListSz);
      nodeDest->resize(tmpCtr);

      // swap pointers ...
      std::vector<TreeNode> *tmpPtr = nodeSrc;
      nodeSrc = nodeDest;
      nodeDest = tmpPtr;
      tmpPtr = leafSrc;
      leafSrc = leafDest;
      leafDest = tmpPtr;
      leafDest->clear();
      nodeDest->clear();  
    } while (!leafSrc->empty());

    if ( nodeSrc != &nodes ) {
      nodes = nodesTmp;
    }
    nodesTmp.clear();

    //Reset All weights to 1.
    for (unsigned int i = 0; i < nodes.size(); i++) {
      nodes[i].setWeight(1);
    }

    PROF_P2O_LOCAL_END
  }

  //New Implementation. Written on April 19th, 2008
  //Both ends are inclusive. The output is sorted.
  int appendCompleteRegion(TreeNode first, TreeNode second, 
      std::vector<ot::TreeNode>& newNodes, bool includeMin, bool includeMax) {

    PROF_COMPLETE_REGION_BEGIN

      unsigned int dim = first.getDim();
    unsigned int maxDepth = first.getMaxDepth();

    TreeNode min = ((first < second) ? first : second);

    if(includeMin) {
      newNodes.push_back(min);
    }

    if (first == second) {
      PROF_COMPLETE_REGION_END
    }

    TreeNode max = ((first > second) ? first : second);

    //Add nodes > min and < max
    TreeNode nca = getNCA(min,max);

    if(min == nca) {
      //special case. Top down approach
      ot::TreeNode tmpAncestor = min;
      bool repeatLoop;
      do {
        repeatLoop = false;
        std::vector<ot::TreeNode> tmpChildList;
        tmpAncestor.addChildren(tmpChildList);
        for(unsigned int j = 0; j < tmpChildList.size(); j++) {
#ifdef __DEBUG_OCT__
          assert(areComparable(tmpChildList[j], max));
#endif
          if( (tmpChildList[j] < max) &&
              (!(tmpChildList[j].isAncestor(max))) ) {
            newNodes.push_back(tmpChildList[j]);
          } else if(tmpChildList[j].isAncestor(max)) {
            tmpAncestor = tmpChildList[j];
            repeatLoop = true;
            break;
          } else {
            assert(tmpChildList[j] == max);
            break;
          }
        }//end for j
      } while (repeatLoop);
    } else {
      TreeNode currentNode = min;    
      while(currentNode > nca) {
        TreeNode parentOfCurrent = currentNode.getParent();
        std::vector<ot::TreeNode> myBros;
        parentOfCurrent.addChildren(myBros);
        for(unsigned int i = 0; i < myBros.size(); i++) {
#ifdef __DEBUG_OCT__
          assert(areComparable(myBros[i], max));
#endif
          if( (myBros[i] > min) &&
              (myBros[i] < max) && (!(myBros[i].isAncestor(max))) ) {
            //Bottom-up here
            newNodes.push_back(myBros[i]);
          } else if(myBros[i].isAncestor(max)) {
            //Top-down now
            //nca will be the parentOfCurrent 
            //If myBros[i] is an acestor of max then
            //it is automatically > min and < max
            //The octants are automatically inserted in
            //the sorted order, due to properties of space-filling curves.
            ot::TreeNode tmpAncestor = myBros[i];
            bool repeatLoop;
            do {
              repeatLoop = false;
              std::vector<ot::TreeNode> tmpChildList;
              tmpAncestor.addChildren(tmpChildList);
              for(unsigned int j = 0; j < tmpChildList.size(); j++) {
#ifdef __DEBUG_OCT__
                assert(areComparable(tmpChildList[j], max));
#endif
                if( (tmpChildList[j] < max) &&
                    (!(tmpChildList[j].isAncestor(max))) ) {
                  newNodes.push_back(tmpChildList[j]);
                } else if(tmpChildList[j].isAncestor(max)) {
                  tmpAncestor = tmpChildList[j];
                  repeatLoop = true;
                  break;
                } else {
                  assert(tmpChildList[j] == max);
                  break;
                }
              }//end for j
            } while (repeatLoop);
            break;	
          }//end select between bottom-up/ top-down
        }//end for i		
        currentNode = parentOfCurrent;
      }//end while	
    }//end if special case

    if(includeMax) {
      newNodes.push_back(max);
    }

    PROF_COMPLETE_REGION_END
  }//end function

  /*
  //Original Code
  //Both ends are inclusive. The output is sorted.
  int completeRegion(TreeNode first, TreeNode second, 
  std::vector<ot::TreeNode>& newNodes, bool includeMin, bool includeMax) {

  PROF_COMPLETE_REGION_BEGIN

  unsigned int dim = first.getDim();
  unsigned int maxDepth = first.getMaxDepth();

  TreeNode min = ( (first < second) ? first : second );
  TreeNode max = ( (first > second) ? first : second );

  newNodes.clear();

  if(includeMin) {
  newNodes.push_back(min);
  }

  if (first == second) {
  return 1;
  }

  if(includeMax) {
  newNodes.push_back(max);
  }

//Add nodes > min and < max
TreeNode nca = getNCA(min,max);

std::vector<TreeNode> workingNodes;
std::vector<TreeNode> addList;
std::vector<TreeNode> * wSrc = &(workingNodes);
std::vector<TreeNode> * wDest = &(addList);
workingNodes.push_back(nca);
while (!(wSrc->empty())) {
wDest->resize((wSrc->size())*(1 << dim));
unsigned int addListCtr=0;
for (unsigned int k=0;k < wSrc->size();k++) {
if ( ((*wSrc)[k] > min) && ((*wSrc)[k] < max) && 
(!(*wSrc)[k].isAncestor(max)) ) {
newNodes.push_back((*wSrc)[k]);
} else if ((*wSrc)[k].getLevel() < maxDepth) {
if ((*wSrc)[k].isAncestor(min) || (*wSrc)[k].isAncestor(max) ) {
std::vector<TreeNode>  tmpChildren;
(*wSrc)[k].addChildren(tmpChildren);
for (int l=0;l<(1 << dim);l++) {
if (tmpChildren[l]<max) {
(*wDest)[addListCtr++] = tmpChildren[l];
}
}//end for l
tmpChildren.clear();
}
}//end if-else
}//end for k
wDest->resize(addListCtr);
std::vector<TreeNode> * tmpPtr = wSrc;
wSrc = wDest;
wDest = tmpPtr;
}//end while

sort(newNodes.begin(),newNodes.end());

PROF_COMPLETE_REGION_END
}//end function
*/

}//end namespace

