/**
  @file Construct.C
  @brief A set of functions for Octree Construction
  @author Hari Sundar, hsundar@gmail.com
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  */

#include "parUtils.h"
#include "TreeNode.h"
#include <cassert>
#include <list>
#include <bits/algorithmfwd.h>
#include "nodeAndValues.h"
#include "binUtils.h"
#include "dendro.h"
#include "testUtils.h"

#include "treenode2vtk.h"

#ifdef __DEBUG__
#ifndef __DEBUG_OCT__
#define __DEBUG_OCT__
#endif
#endif

namespace ot {

  int regularGrid2Octree(const std::vector<double> &elementValues,
                         unsigned int N, unsigned int nx, unsigned int ny, unsigned int nz,
                         unsigned int xs, unsigned int ys, unsigned int zs, std::vector<TreeNode> &linOct,
                         unsigned int dim, unsigned int maxDepth, double thresholdFac, MPI_Comm comm) {
    PROF_RG2O_BEGIN

    int rank;
    int npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    const int MIN_GRAIN_SIZE = 10;

    bool isNvalid = binOp::isPowerOfTwo(N);

    assert(isNvalid);

    unsigned int rgLevel = binOp::fastLog2(N);
    unsigned int elemLen = (1u << (maxDepth - rgLevel));

    //1. Convert input to regular octree.
    std::vector<ot::NodeAndValues<double, 1> > tmpList(nx * ny * nz);

    for (int k = 0; k < nz; k++) {
      for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
          unsigned int currX = ((i + xs) * elemLen);
          unsigned int currY = ((j + ys) * elemLen);
          unsigned int currZ = ((k + zs) * elemLen);
          unsigned int idx = ((nx * (j + (ny * k))) + i);
          tmpList[idx].node = ot::TreeNode(currX, currY, currZ, rgLevel, dim, maxDepth);
          tmpList[idx].values[0] = elementValues[idx];
        } //end for i
      } //end for j
    } //end for k

    //2. Convert to Morton Ordering
    std::vector<ot::NodeAndValues<double, 1> > tnAndValsList;
    par::sampleSort<ot::NodeAndValues<double, 1> >(tmpList, tnAndValsList, comm);
    tmpList.clear();

    DendroIntL inSz = tnAndValsList.size();
    DendroIntL globInSize;
    par::Mpi_Allreduce<DendroIntL>(&inSz, &globInSize, 1, MPI_SUM, comm);

    bool repeatLoop = true;
    int npesCurr = npes;
    MPI_Comm commCurr = comm;
    if (globInSize < (MIN_GRAIN_SIZE * npes)) {
      int splittingSize = (globInSize / MIN_GRAIN_SIZE);
      if (splittingSize == 0) {
        splittingSize = 1;
      }

      unsigned int avgLoad = (globInSize / splittingSize);
      int leftOvers = (globInSize - (splittingSize * avgLoad));

      if (rank >= splittingSize) {
        par::scatterValues<ot::NodeAndValues<double, 1> >(tnAndValsList, tmpList, 0, commCurr);
      } else if (rank < leftOvers) {
        par::scatterValues<ot::NodeAndValues<double, 1> >(tnAndValsList, tmpList, (avgLoad + 1), commCurr);
      } else {
        par::scatterValues<ot::NodeAndValues<double, 1> >(tnAndValsList, tmpList, avgLoad, commCurr);
      }
      tnAndValsList = tmpList;
      tmpList.clear();

      MPI_Comm newComm;
      par::splitCommUsingSplittingRank(splittingSize, &newComm, commCurr);
      commCurr = newComm;
      npesCurr = splittingSize;

      if (rank >= splittingSize) {
        repeatLoop = false;
      }
    } else {
      par::partitionW<ot::NodeAndValues<double, 1> >(tnAndValsList, NULL, comm);
    }

    if (globInSize < MIN_GRAIN_SIZE) {
      repeatLoop = false;
    }

    while (repeatLoop) {

      inSz = tnAndValsList.size();

      assert(inSz >= MIN_GRAIN_SIZE);

      MPI_Request requests[4];

      //Send first 7 to previous processor
      //Send last 7 to next processor
      ot::NodeAndValues<double, 1> prevOcts[7];
      if (rank) {
        par::Mpi_Irecv<ot::NodeAndValues<double, 1> >(prevOcts, 7, (rank - 1),
                                                      1, commCurr, &(requests[0]));
      }

      ot::NodeAndValues<double, 1> nextOcts[7];
      if (rank < (npesCurr - 1)) {
        par::Mpi_Irecv<ot::NodeAndValues<double, 1> >(nextOcts, 7, (rank + 1),
                                                      1, commCurr, &(requests[1]));
      }

      if (rank) {
        par::Mpi_Issend<ot::NodeAndValues<double, 1> >((&(*(tnAndValsList.begin()))), 7,
                                                       (rank - 1), 1, commCurr, &(requests[2]));
      }

      if (rank < (npesCurr - 1)) {
        par::Mpi_Issend<ot::NodeAndValues<double, 1> >((&(*(tnAndValsList.end() - 7))), 7,
                                                       (rank + 1), 1, commCurr, &(requests[3]));
      }

      MPI_Status statuses[4];

      if (rank) {
        MPI_Wait(&(requests[0]), &(statuses[0]));
        MPI_Wait(&(requests[2]), &(statuses[2]));
      }

      if (rank < (npesCurr - 1)) {
        MPI_Wait(&(requests[1]), &(statuses[1]));
        MPI_Wait(&(requests[3]), &(statuses[3]));
      }

      //Check and coarsen
      int idxOfPrevFC = -1;
      if (rank) {
        for (int i = 0; i < 7; i++) {
          //There may be more than 1 FC, we only need the last one here
          if (prevOcts[i].node.getChildNumber() == 0) {
            idxOfPrevFC = i;
          }
        } //end for i
      }

      int idxOfNextFC = -1;
      if (rank < (npesCurr - 1)) {
        for (int i = 0; i < 7; i++) {
          //There may be more than 1 FC, we only need the first one here
          if (nextOcts[i].node.getChildNumber() == 0) {
            idxOfNextFC = i;
            break;
          }
        } //end for i
      }

      int myFirstFC = -1;
      for (int i = 0; i < inSz; i++) {
        if (tnAndValsList[i].node.getChildNumber() == 0) {
          myFirstFC = i;
          break;
        }
      } //end for i

      int myLastFC = -1;
      if (myFirstFC >= 0) {
        for (int i = (inSz - 1); i >= 0; i--) {
          if (tnAndValsList[i].node.getChildNumber() == 0) {
            myLastFC = i;
            break;
          }
        } //end for i
      }

      //Process upto myFirstFC (exclusive)
      //Note, we check the threshold on both the processors so that the
      //partition is not changed unnecessarily.
      if ((myFirstFC >= 0) && (idxOfPrevFC >= 0)) {
        int fcGap = (myFirstFC + 7 - idxOfPrevFC);
        if (fcGap < 8) {
          //Can not coarsen
          if (myFirstFC) {
            tmpList.insert(tmpList.end(), tnAndValsList.begin(),
                           tnAndValsList.begin() + myFirstFC);
          }
        } else {
          //fcGap >= 8
          //Can coarsen upto (excluding) (1 + idxOfPrevFC)
          double minVal = prevOcts[idxOfPrevFC].values[0];
          double maxVal = prevOcts[idxOfPrevFC].values[0];
          for (int i = idxOfPrevFC; i < 7; i++) {
            double currVal = prevOcts[i].values[0];
            if (currVal < minVal) {
              minVal = currVal;
            }
            if (currVal > maxVal) {
              maxVal = currVal;
            }
          } //end for i
          for (int i = 0; i < (1 + idxOfPrevFC); i++) {
            double currVal = tnAndValsList[i].values[0];
            if (currVal < minVal) {
              minVal = currVal;
            }
            if (currVal > maxVal) {
              maxVal = currVal;
            }
          } //end for i
          if ((maxVal - minVal) >= thresholdFac) {
            //Can't coarsen
            if (myFirstFC) {
              tmpList.insert(tmpList.end(), tnAndValsList.begin(),
                             tnAndValsList.begin() + myFirstFC);
            }
          } else {
            //Can coarsen. The previous processor will add the parent
            if ((1 + idxOfPrevFC) < myFirstFC) {
              tmpList.insert(tmpList.end(), tnAndValsList.begin() + (1 + idxOfPrevFC),
                             tnAndValsList.begin() + myFirstFC);
            } //if to skip equality
          } //end check threshold
        } //end check fcGap
      } else {
        if (myFirstFC >= 0) {
          //Can not coarsen
          //The prev processor has no FC. Since the prev processor has more
          //than 8 elements, our elements will remain.
          if (myFirstFC) {
            tmpList.insert(tmpList.end(), tnAndValsList.begin(),
                           tnAndValsList.begin() + myFirstFC);
          }
        } else if (idxOfPrevFC >= 0) {
          //Can coarsen upto (excluding) (1 + idxOfPrevFC)
          double minVal = prevOcts[idxOfPrevFC].values[0];
          double maxVal = prevOcts[idxOfPrevFC].values[0];
          for (int i = idxOfPrevFC; i < 7; i++) {
            double currVal = prevOcts[i].values[0];
            if (currVal < minVal) {
              minVal = currVal;
            }
            if (currVal > maxVal) {
              maxVal = currVal;
            }
          } //end for i
          for (int i = 0; i < (1 + idxOfPrevFC); i++) {
            double currVal = tnAndValsList[i].values[0];
            if (currVal < minVal) {
              minVal = currVal;
            }
            if (currVal > maxVal) {
              maxVal = currVal;
            }
          } //end for i
          if ((maxVal - minVal) >= thresholdFac) {
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
      for (int idx = myFirstFC + 1; idx <= myLastFC; idx++) {
        if (tnAndValsList[idx].node.getChildNumber() == 0) {
          int fcGap = (idx - prevFCidx);
          if (fcGap < 8) {
            //Can't coarsen
            tmpList.insert(tmpList.end(), tnAndValsList.begin() + prevFCidx,
                           tnAndValsList.begin() + idx);
          } else {
            //fcGap >= 8
            //Can coarsen upto (excluding) (8 + prevFCidx)
            double minVal = tnAndValsList[prevFCidx].values[0];
            double maxVal = tnAndValsList[prevFCidx].values[0];
            double sumVal = 0;
            for (int j = prevFCidx; j < (8 + prevFCidx); j++) {
              double currVal = tnAndValsList[j].values[0];
              if (currVal < minVal) {
                minVal = currVal;
              }
              if (currVal > maxVal) {
                maxVal = currVal;
              }
              sumVal += currVal;
            } //end for j
            if ((maxVal - minVal) >= thresholdFac) {
              //Can't coarsen.
              tmpList.insert(tmpList.end(), tnAndValsList.begin() + prevFCidx,
                             tnAndValsList.begin() + idx);
            } else {
              //Can coarsen.
              ot::NodeAndValues<double, 1> tmpObj;
              tmpObj.node = tnAndValsList[prevFCidx].node.getParent();
              tmpObj.values[0] = (sumVal / 8.0);
              tmpList.push_back(tmpObj);
              if ((prevFCidx + 8) < idx) {
                tmpList.insert(tmpList.end(), tnAndValsList.begin() + prevFCidx + 8,
                               tnAndValsList.begin() + idx);
              } //if to skip equality
            }
          } //end check fcGap
          prevFCidx = idx;
        } //end if found FC
      } //end for idx

      //Process myLastFC (inclusive) to end
      if ((myLastFC >= 0) && (idxOfNextFC >= 0)) {
        int fcGap = inSz + idxOfNextFC - myLastFC;
        if (fcGap < 8) {
          //Can't coarsen
          tmpList.insert(tmpList.end(), tnAndValsList.begin() + myLastFC,
                         tnAndValsList.end());
        } else {
          //fcGap >= 8
          //Can coarsen upto (excluding) (myLastFC + 8)
          double minVal = tnAndValsList[myLastFC].values[0];
          double maxVal = tnAndValsList[myLastFC].values[0];
          double sumVal = 0;
          if ((myLastFC + 8) > inSz) {
            for (int i = myLastFC; i < inSz; i++) {
              double currVal = tnAndValsList[i].values[0];
              if (currVal < minVal) {
                minVal = currVal;
              }
              if (currVal > maxVal) {
                maxVal = currVal;
              }
              sumVal += currVal;
            } //end for i
            for (int i = 0; i < (myLastFC + 8 - inSz); i++) {
              double currVal = nextOcts[i].values[0];
              if (currVal < minVal) {
                minVal = currVal;
              }
              if (currVal > maxVal) {
                maxVal = currVal;
              }
              sumVal += currVal;
            } //end for i
          } else {
            for (int i = myLastFC; i < (myLastFC + 8); i++) {
              double currVal = tnAndValsList[i].values[0];
              if (currVal < minVal) {
                minVal = currVal;
              }
              if (currVal > maxVal) {
                maxVal = currVal;
              }
              sumVal += currVal;
            } //end for i
          }
          if ((maxVal - minVal) >= thresholdFac) {
            //Can't coarsen
            tmpList.insert(tmpList.end(), tnAndValsList.begin() + myLastFC,
                           tnAndValsList.end());
          } else {
            //Can coarsen
            ot::NodeAndValues<double, 1> tmpObj;
            tmpObj.node = tnAndValsList[myLastFC].node.getParent();
            tmpObj.values[0] = (sumVal / 8.0);
            tmpList.push_back(tmpObj);
            if ((myLastFC + 8) < inSz) {
              tmpList.insert(tmpList.end(), tnAndValsList.begin() + myLastFC + 8,
                             tnAndValsList.end());
            }
          } //end check threshold
        } //end check fcGap
      } else {
        if (myLastFC >= 0) {
          //Can coarsen upto (excluding) (myLastFC + 8)
          double minVal = tnAndValsList[myLastFC].values[0];
          double maxVal = tnAndValsList[myLastFC].values[0];
          double sumVal = 0;
          if ((myLastFC + 8) > inSz) {
            for (int i = myLastFC; i < inSz; i++) {
              double currVal = tnAndValsList[i].values[0];
              if (currVal < minVal) {
                minVal = currVal;
              }
              if (currVal > maxVal) {
                maxVal = currVal;
              }
              sumVal += currVal;
            } //end for i
            if (rank < (npesCurr - 1)) {
              for (int i = 0; i < (myLastFC + 8 - inSz); i++) {
                double currVal = nextOcts[i].values[0];
                if (currVal < minVal) {
                  minVal = currVal;
                }
                if (currVal > maxVal) {
                  maxVal = currVal;
                }
                sumVal += currVal;
              } //end for i
            } //end if last proc
          } else {
            for (int i = myLastFC; i < (myLastFC + 8); i++) {
              double currVal = tnAndValsList[i].values[0];
              if (currVal < minVal) {
                minVal = currVal;
              }
              if (currVal > maxVal) {
                maxVal = currVal;
              }
              sumVal += currVal;
            } //end for i
          }
          if ((maxVal - minVal) >= thresholdFac) {
            //Can't coarsen
            tmpList.insert(tmpList.end(), tnAndValsList.begin() + myLastFC,
                           tnAndValsList.end());
          } else {
            //Can coarsen
            ot::NodeAndValues<double, 1> tmpObj;
            tmpObj.node = tnAndValsList[myLastFC].node.getParent();
            tmpObj.values[0] = (sumVal / 8.0);
            tmpList.push_back(tmpObj);
            if ((myLastFC + 8) < inSz) {
              tmpList.insert(tmpList.end(), tnAndValsList.begin() + myLastFC + 8,
                             tnAndValsList.end());
            }
          }
        } else {
          //This case is the same as myFirstFC < 0 and it has been taken care
          //of earlier.
        }
      }

      DendroIntL outSz = tmpList.size();
      DendroIntL globOutSize;
      par::Mpi_Allreduce<DendroIntL>(&outSz, &globOutSize, 1, MPI_SUM, commCurr);

      if (globOutSize == globInSize) {
        repeatLoop = false;
      }

      if (globOutSize < MIN_GRAIN_SIZE) {
        repeatLoop = false;
      }

      //prepare for next iteration...
      globInSize = globOutSize;
      tnAndValsList = tmpList;
      tmpList.clear();

      if (globInSize < (MIN_GRAIN_SIZE * npesCurr)) {
        int splittingSize = (globInSize / MIN_GRAIN_SIZE);
        if (splittingSize == 0) {
          splittingSize = 1;
        }

        unsigned int avgLoad = (globInSize / splittingSize);
        int leftOvers = (globInSize - (splittingSize * avgLoad));

        if (rank >= splittingSize) {
          par::scatterValues<ot::NodeAndValues<double, 1> >(tnAndValsList, tmpList, 0, commCurr);
        } else if (rank < leftOvers) {
          par::scatterValues<ot::NodeAndValues<double, 1> >(tnAndValsList, tmpList, (avgLoad + 1), commCurr);
        } else {
          par::scatterValues<ot::NodeAndValues<double, 1> >(tnAndValsList, tmpList, avgLoad, commCurr);
        }
        tnAndValsList = tmpList;
        tmpList.clear();

        MPI_Comm newComm;
        par::splitCommUsingSplittingRank(splittingSize, &newComm, commCurr);
        commCurr = newComm;
        npesCurr = splittingSize;

        if (rank >= splittingSize) {
          repeatLoop = false;
        }
      } else {
        par::partitionW<ot::NodeAndValues<double, 1> >(tnAndValsList, NULL, commCurr);
      } //end if splitting comm

    } //end while

    linOct.clear();
    for (int i = 0; i < tnAndValsList.size(); i++) {
      linOct.push_back(tnAndValsList[i].node);
    } //end for i

    PROF_RG2O_END
  } //end function

//New Implementation. Written on April 20, 2008
  int completeOctree(const std::vector<TreeNode> &inp,
                     std::vector<TreeNode> &out,
                     unsigned int dim, unsigned int maxDepth, bool isUnique,
                     bool isSorted, bool assertNoEmptyProcs, MPI_Comm comm) {

    PROF_N2O_BEGIN

    TreeNode root(dim, maxDepth);

    int size;
    MPI_Comm_size(comm, &size);

    if (size == 1) {
      completeSubtree(root, inp, out, dim, maxDepth, isUnique, isSorted);
      PROF_N2O_END
    } //end single proc case

    out = inp;

    //Sort and remove duplicate leaves.
    if (!isUnique) {
      par::removeDuplicates<ot::TreeNode>(out, isSorted, comm);
    } else if (!isSorted) {
      std::vector<TreeNode> tmpOut;
      par::sampleSort<ot::TreeNode>(out, tmpOut, comm);
      out = tmpOut;
      tmpOut.clear();
    }

    //Remove empty processors...
    MPI_Comm new_comm;
    // quick and dirty fix to avoid repetetive creation of communicators (which would exhaust MPI resources)
    if (true /* assertNoEmptyProcs */) {
      new_comm = comm;
      assert(!out.empty());
    } else {
      par::splitComm2way(out.empty(), &new_comm, comm);
    }

    if (!(out.empty())) {
      int new_rank, new_size;

      MPI_Comm_rank(new_comm, &new_rank);
      MPI_Comm_size(new_comm, &new_size);

      MPI_Request requestSend;
      MPI_Request requestRecv;

      ot::TreeNode begBuf;
      ot::TreeNode lastElem;

      if (new_rank) {
        //Recv
        par::Mpi_Irecv<ot::TreeNode>(&begBuf, 1, (new_rank - 1), 1, new_comm, &requestRecv);
      } //end if not P0

      if (new_rank < (new_size - 1)) {
        lastElem = out[out.size() - 1];
        //Send
        par::Mpi_Issend<ot::TreeNode>(&lastElem, 1, (new_rank + 1), 1, new_comm, &requestSend);
      } //end if not PN


      //Add missing corners to complete the region.
      //Add the first corner leaf on the first processor.
      if (new_rank == 0) {
        // @milinda is this correct?
        ot::TreeNode minCorner(0, 0, 0, maxDepth, dim, maxDepth);
#ifdef __DEBUG_OCT__
        assert(areComparable(out[0], minCorner));
#endif
        if ((out[0] != minCorner) && (!out[0].isAncestor(minCorner))) {
          ot::TreeNode ncaTmp = getNCA(out[0], minCorner);
          std::vector<ot::TreeNode> kids;
          ncaTmp.addChildren(kids);
          out.insert(out.begin(), kids[0]);
          kids.clear();
        } //end if
      } //end if

      //Add the last corner leaf on the last processor.
      if (new_rank == (new_size - 1)) {
        // @author hari sundar: only works for Morton, so changed to root.getDLD()
        // ot::TreeNode maxCorner(((1u << maxDepth) - 1), ((1u << maxDepth) - 1), ((1u << maxDepth) - 1), maxDepth, dim, maxDepth);
        ot::TreeNode maxCorner = root.getDLD();
#ifdef __DEBUG_OCT__
        assert(areComparable(out[out.size() - 1], maxCorner));
#endif
        if ((out[out.size() - 1] != maxCorner) &&
            (!out[out.size() - 1].isAncestor(maxCorner))) {
          ot::TreeNode ncaTmp = getNCA(out[out.size() - 1], maxCorner);
          std::vector<ot::TreeNode> kids;
          ncaTmp.addChildren(kids);
          out.insert(out.end(), kids[(1 << dim) - 1]);
          kids.clear();
        } //end if
      } //end if

      std::vector<ot::TreeNode> tmpList;
      for (unsigned int i = 0; i < (out.size() - 1); i++) {
#ifdef __DEBUG_OCT__
        assert(areComparable(out[i], out[i + 1]));
#endif
        if (out[i].isAncestor(out[i + 1])) {
          appendCompleteRegion(out[i], out[i + 1], tmpList, false, false);
        } else {
          appendCompleteRegion(out[i], out[i + 1], tmpList, true, false);
        }
      } //end for i

      //Only the last processor adds the last element. All the other processors would have
      //sent it to the next processor, which will add it if it is not an ancestor of
      //the first element on that processor
      if (new_rank == (new_size - 1)) {
        tmpList.push_back(out[out.size() - 1]);
      }

      if (new_rank) {
        MPI_Status statusWait;
        MPI_Wait(&requestRecv, &statusWait);

        std::vector<ot::TreeNode> begList;
#ifdef __DEBUG_OCT__
        assert(areComparable(begBuf, out[0]));
#endif
        if (begBuf.isAncestor(out[0])) {
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

      if (new_rank < (new_size - 1)) {
        MPI_Status statusWait;
        MPI_Wait(&requestSend, &statusWait);
      }

    } //out not empty


    PROF_N2O_END

  } //end function

//New Implementation. Written on April 19, 2008
  int completeSubtree(TreeNode block, const std::vector<TreeNode> &inp, std::vector<TreeNode> &out,
                      unsigned int dim, unsigned int maxDepth, bool isUnique, bool isSorted) {

    PROF_N2O_SEQ_BEGIN

    out = inp;
    //Sort and remove duplicate leaves.
    if (!isUnique) {
      seq::makeVectorUnique<TreeNode>(out, isSorted);
    } else if (!isSorted) {
      omp_par::merge_sort(out.begin(), out.end());
    }

    if (!out.empty()) {
      //Add missing corners to complete the region.
      //Add the first corner .
      ot::TreeNode minCorner = block.getDFD();
#ifdef __DEBUG_OCT__
      assert(areComparable(out[0], minCorner));
#endif
      if ((out[0] != minCorner) && (!out[0].isAncestor(minCorner))) {
        ot::TreeNode ncaTmp = getNCA(out[0], minCorner);
        std::vector<ot::TreeNode> kids;
        ncaTmp.addChildren(kids);
        out.insert(out.begin(), kids[0]);
        kids.clear();
      } //end if

      //Add the last corner.
      ot::TreeNode maxCorner = block.getDLD();
#ifdef __DEBUG_OCT__
      assert(areComparable(out[out.size() - 1], maxCorner));
#endif
      if ((out[out.size() - 1] != maxCorner) && (!out[out.size() - 1].isAncestor(maxCorner))) {
        ot::TreeNode ncaTmp = getNCA(out[out.size() - 1], maxCorner);
        std::vector<ot::TreeNode> kids;
        ncaTmp.addChildren(kids);
        out.insert(out.end(), kids[(1 << dim) - 1]);
        kids.clear();
      } //end if

      std::vector<ot::TreeNode> tmpList;
      for (unsigned int i = 0; i < (out.size() - 1); i++) {
#ifdef __DEBUG_OCT__
        assert(areComparable(out[i], out[i + 1]));
#endif
        if (out[i].isAncestor(out[i + 1])) {
          appendCompleteRegion(out[i], out[i + 1], tmpList, false, false);
        } else {
          appendCompleteRegion(out[i], out[i + 1], tmpList, true, false);
        }
      } //end for i

      tmpList.push_back(out[out.size() - 1]);

      out = tmpList;
      tmpList.clear();
    } //end if out empty

    PROF_N2O_SEQ_END

  } //end function

  int points2Octree(std::vector<double> &pts, double *gLens, std::vector<ot::TreeNode> &nodes,
                    unsigned int dim, unsigned int maxDepth, unsigned int maxNumPts, MPI_Comm comm) {

    PROF_P2O_BEGIN

    int size;
    MPI_Comm_size(comm, &size);

    //Sequential...
    if (size == 1) {
      std::cout << YLW " - Calling p2o Seq " NRM << std::endl;
      points2OctreeSeq(pts, gLens, nodes, dim, maxDepth, maxNumPts);
      PROF_P2O_END
    } //end if sequential

    int rank;
    MPI_Comm_rank(comm, &rank);

    TreeNode root(dim, maxDepth);

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
    assert((dim == 1) || (dim == 2) || (dim == 3));
    assert((ptsLen % dim) == 0);

    DendroIntL numNodes = (ptsLen / dim);

    DendroIntL totSize;
    par::Mpi_Allreduce<DendroIntL>(&numNodes, &totSize, 1, MPI_SUM, comm);

    if (totSize <= static_cast<DendroIntL>(maxNumPts)) {
      if (!rank) {
        nodes.resize(1);
        nodes[0] = root;
      } else {
        nodes.resize(0);
      } //end if-else
      PROF_P2O_END
    } //end if

    //Tackle small problems separately....
    //min Grain size = 1000
    // for now, we shall turn this feature off (since it is a hassle for FMM code)
    const DendroIntL THOUSAND = 1;
    if (totSize < (THOUSAND * size)) {
      int splittingSize = (totSize / THOUSAND);

      if (splittingSize == 0) {
        splittingSize = 1;
      }

      unsigned int avgLoad = (totSize / splittingSize);
      int leftOvers = (totSize - (splittingSize * avgLoad));

      std::vector<double> tmpPts;
      if (rank >= splittingSize) {
        par::scatterValues<double>(pts, tmpPts, 0, comm);
      } else if (rank < leftOvers) {
        par::scatterValues<double>(pts, tmpPts, (dim * (avgLoad + 1)), comm);
      } else {
        par::scatterValues<double>(pts, tmpPts, (dim * avgLoad), comm);
      }
      pts.clear();

      MPI_Comm newComm;
      par::splitCommUsingSplittingRank(splittingSize, &newComm, comm);

#ifndef __SILENT_MODE__
      if (!rank) {
        std::cout << YLW " Input to p2o is small (" << totSize
        << "). npes = " << size << " Splitting Comm. " NRM << std::endl;
      }
#endif

      if (rank < splittingSize) {
        points2Octree(tmpPts, gLens, nodes, dim, maxDepth, maxNumPts, newComm);
      }
      tmpPts.clear();

      PROF_P2O_END
    } //reduce procs for small problems

    //Tackle large problems....

    nodes.resize(numNodes);

    double scale[3] = {0, 0, 0};
    for (int i = 0; i < dim; i++) scale[i] = ((double) (1u << maxDepth)) / (gLens[i]);
#pragma omp parallel for
    for (DendroIntL i = 0; i < numNodes; i++) {
      //The constructor will ignore unnecessary arguments (for lower
      //dimensions).
      unsigned int p_int[3] = {0, 0, 0};
      for (int j = 0; j < dim; j++) p_int[j] = (unsigned int) (pts[i * dim + j] * scale[j]);
      nodes[i] = TreeNode(p_int[0], p_int[1], p_int[2], maxDepth, dim, maxDepth);
    } //end for
    pts.clear();

    //Sort nodes (pts.) and partition them.


    // std::cout << rank << "before sample sort: " << nodes.size() << std::endl;
    std::vector<ot::TreeNode> tmpNodes;
    //treeNodesTovtk(nodes,rank,"bf_SS");
    par::sampleSort<ot::TreeNode>(nodes, tmpNodes, comm);
    // std::cout << rank << "after sample sort: " << tmpNodes.size() << std::endl;
    nodes.clear();

//   for(int i=0;i<tmpNodes.size();i++)
//     nodes[i] = tmpNodes[i];
    nodes = tmpNodes;
    //tmpNodes.clear();
    //treeNodesTovtk(nodes,rank,"af_SS");
    std::vector<ot::TreeNode> leaves;
    std::vector<ot::TreeNode> minsAllBlocks;

    // if (!rank) std::cout << RED " Before BlkPart: " NRM << std::endl;

    //assert(par::test::isUniqueAndSorted(nodes, comm));
    // if (nodes.size() > (1 << dim) ) {
    blockPartStage1_p2o(nodes, leaves, dim, maxDepth, comm);
    blockPartStage2_p2o(nodes, leaves, minsAllBlocks, dim, maxDepth, comm);


    assert(par::test::isUniqueAndSorted(leaves, comm));

    p2oLocal(nodes, leaves, maxNumPts, dim, maxDepth);

    PROF_P2O_END

  } //end function

//Added on April 19, 2008
  int points2OctreeSeq(std::vector<double> &pts, double *gLens, std::vector<TreeNode> &nodes,
                       unsigned int dim, unsigned int maxDepth, unsigned int maxNumPts) {

    PROF_P2O_SEQ_BEGIN

    if (maxDepth == 0) {
      nodes.resize(1);
      nodes[0] = TreeNode(dim, maxDepth);
      PROF_P2O_SEQ_END
    }

    unsigned int ptsLen = pts.size();
    assert((dim == 1) || (dim == 2) || (dim == 3));
    assert((ptsLen % dim) == 0);

    TreeNode root(dim, maxDepth);
    unsigned int numNodes = ptsLen / dim;
    nodes.resize(numNodes);

    for (unsigned int i = 0; i < numNodes; i++) {
      //The constructor will ignore unnecessary arguments (for lower
      //dimensions).
      unsigned int px = (unsigned int) (pts[i * dim] * ((double) (1u << maxDepth)) / (gLens[0]));
      unsigned int py, pz = 0;
      if (dim > 1) {
        py = (unsigned int) (pts[(i * dim) + 1] * ((double) (1u << maxDepth)) / gLens[1]);
        if (dim > 2) {
          pz = (unsigned int) (pts[(i * dim) + 2] * ((double) (1u << maxDepth)) / gLens[2]);
        }
      }
      nodes[i] = TreeNode(px, py, pz, maxDepth, dim, maxDepth);
    } //end for

    pts.clear();

    unsigned int totSize = nodes.size();

    if (totSize <= maxNumPts) {
      nodes.resize(1);
      nodes[0] = root;
      PROF_P2O_SEQ_END
    } //end if

    //Sort nodes (pts.) and partition them.
    omp_par::merge_sort(nodes.begin(), nodes.end());

    // std::cout << CYN " - After mergesort, is sorted and unique: " NRM << seq::test::isUniqueAndSorted(nodes) <<
    // std::endl;

    std::vector<TreeNode> leaves;
    leaves.push_back(root);

    treeNodesTovtk(nodes, 0, "input_p2o");
    p2oLocal(nodes, leaves, maxNumPts, dim, maxDepth);

    PROF_P2O_SEQ_END

  } //end function




/**
 * @author Dhairya Malhotra, dhairya.malhotra88@gmail.com
 * @date 08 Feb 2010
 */
  int p2oLocal(std::vector<TreeNode> &nodes, std::vector<TreeNode> &leaves,
               unsigned int maxNumPts, unsigned int dim, unsigned int maxDepth) {
      PROF_P2O_LOCAL_BEGIN;

      std::cout << "entering p2o_local=====" << std::endl;


      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);

      std::list<TreeNode> leaves_lst;
      std::vector<TreeNode>leaves_lst_v;

      unsigned int init_size = leaves.size();
      unsigned int num_pts = nodes.size();

      TreeNode curr_node = leaves[0];
      TreeNode last_node = leaves[init_size - 1].getDLD();
      TreeNode next_node = curr_node.getNext();

      unsigned int curr_pt = 0;
      unsigned int next_pt = curr_pt + maxNumPts;


      //Note: for debug purposes
//      std::cout<<"While Loop Begining."<<std::endl;
//      std::cout<<YLW<<"current node:\t"<<curr_node<<NRM<<std::endl;
//      std::cout<<YLW<<"next_node:\t"<<next_node<<NRM<<std::endl;
//      std::cout<<YLW<<"current point:\t"<<nodes[curr_pt]<<NRM<<std::endl;
//      std::cout<<YLW<<"next_point:\t"<<nodes[next_pt]<<NRM<<std::endl;



      while (next_pt < num_pts) {
        int lev_count=0;
        while (curr_node.isAncestor(nodes[next_pt-1]) & curr_node.isAncestor(nodes[next_pt]) & curr_node.getLevel()<maxDepth) {
          //@hari: We need sure to make that our logic is correct
          //@hari: Whether it is curr_node.isAncestor(nodes[next_pt]) or combination of both. or these are equal since the nodes are sorted.
            curr_node = curr_node.getFirstChild();
            next_node = curr_node.getNext();
            lev_count++;

        }


        unsigned int inc = maxNumPts;
        while (lev_count==(maxDepth-1)) {
          // We have more than maxNumPts points per octant because the node can
          // not be refined any further.
          inc = inc << 1;
          next_pt += inc;
          if (next_pt > num_pts) {
            next_pt = num_pts-1;
            break;
          }
        }

//        bool state_in_between=((nodes[curr_pt]<next_node) & (next_node<nodes[next_pt]));
//        std::cout<<RED<<"In between current point and next point:"<<state_in_between<<std::endl;

        int found_pt=0;
        found_pt=(std::lower_bound(&nodes[curr_pt], &nodes[next_pt], next_node, std::less<TreeNode>()) - &nodes[curr_pt]);
        next_pt = curr_pt + found_pt;
        leaves_lst.push_back(curr_node);

        curr_node = next_node;
        next_node = curr_node.getNext();

        if (next_pt > curr_pt) curr_pt = next_pt;
        next_pt = curr_pt + maxNumPts;


      }


      while (curr_node < last_node) {
        while (curr_node.getDLD() > last_node && curr_node.getLevel() < maxDepth) curr_node = curr_node.getFirstChild();
        leaves_lst.push_back(curr_node);
        if (curr_node.getDLD() == last_node) break;
        curr_node = curr_node.getNext();
      }


      nodes.resize(leaves_lst.size());
      unsigned int i = 0;
      for (std::list<TreeNode>::iterator it = leaves_lst.begin(); it != leaves_lst.end(); it++) {
        nodes[i] = (*it);
        i++;
      }
      leaves_lst.clear();

      std::cout << rank << ": leaving p2o_local" << std::endl;
      PROF_P2O_LOCAL_END
  } //*/




//New Implementation. Written on April 19th, 2008
//Both ends are inclusive. The output is sorted.
  int appendCompleteRegion(TreeNode first, TreeNode second,
                           std::vector<ot::TreeNode> &newNodes, bool includeMin, bool includeMax) {

    PROF_COMPLETE_REGION_BEGIN

    // std::cout << "entering " << __func__ << std::endl;

    unsigned int dim = first.getDim();
    unsigned int maxDepth = first.getMaxDepth();

    TreeNode min = ((first < second) ? first : second);

    if (includeMin) {
      newNodes.push_back(min);
    }

    if (first == second) {
      PROF_COMPLETE_REGION_END
    }

    TreeNode max = ((first > second) ? first : second);
    //Add nodes > min and < max
    TreeNode nca = getNCA(min, max);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    if (min == nca) {
      //special case. Top down approach
      std::cout << "Special Case: min==nca" << std::endl;
      ot::TreeNode tmpAncestor = min;
      bool repeatLoop;
      do {
        repeatLoop = false;
        std::vector<ot::TreeNode> tmpChildList;
        tmpAncestor.addChildren(tmpChildList);
        for (unsigned int j = 0; j < tmpChildList.size(); j++) {
#ifdef __DEBUG_OCT__
          assert(areComparable(tmpChildList[j], max));
#endif
          if ((tmpChildList[j] < max) &&
              (!(tmpChildList[j].isAncestor(max)))) {
            newNodes.push_back(tmpChildList[j]);
          } else if (tmpChildList[j].isAncestor(max)) {
            tmpAncestor = tmpChildList[j];
            repeatLoop = true;
            break;
          } else {
            assert(tmpChildList[j] == max);
            break;
          }
        } //end for j
      }
      while (repeatLoop);
    } else {

      // std::cout<<"nca!=min case"<<std::endl;
      TreeNode currentNode = min;
      while (nca<currentNode /*|| (nca.isAncestor(currentNode) && nca!=currentNode)*/) {
        TreeNode parentOfCurrent = currentNode.getParent();
        // if (!rank) std::cout << "Rank:" << rank << " Parent Node:" << parentOfCurrent << std::endl;
        std::vector<ot::TreeNode> myBros;
        parentOfCurrent.addChildren(myBros);
        for (unsigned int i = 0; i < myBros.size(); i++) {
#ifdef __DEBUG_OCT__
          assert(areComparable(myBros[i], max));
#endif
          if ( (myBros[i] > min) && (myBros[i] < max) && (!(myBros[i].isAncestor(max)))  ) {
            //Bottom-up here

            // if (!rank) std::cout << rank << " adding to new nodes" << myBros[i] << std::endl;
            newNodes.push_back(myBros[i]);
          } else if (myBros[i].isAncestor(max)) {

            // if (!rank) std::cout << rank << " Found ancesstor of max:" << myBros[i] << std::endl;
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
              for (unsigned int j = 0; j < tmpChildList.size(); j++) {
#ifdef __DEBUG_OCT__
                assert(areComparable(tmpChildList[j], max));
#endif
                if ((tmpChildList[j] < max) &&
                    (!(tmpChildList[j].isAncestor(max)))) {
                  newNodes.push_back(tmpChildList[j]);
                  // if (!rank) std::cout << "Adding child " << j << " to newnode " << tmpChildList[j] << std::endl;

                } else if (tmpChildList[j].isAncestor(max)) {
                  tmpAncestor = tmpChildList[j];
                  repeatLoop = true;
                  // if (!rank) std::cout << rank << " Repeating Loop " << tmpAncestor.getLevel() << std::endl;

                  break;
                } else {

//                 if (tmpChildList[j] != max) {
//                   std::cout << "min " << min << " " << min.getMaxDepth() << std::endl;
//                   std::cout << "max " << max << " " << max.getMaxDepth() << std::endl;
//                   std::cout << "tmp " << tmpChildList[j] << " " << tmpChildList[j].getMaxDepth() << std::endl;
//                   bool status=(min>tmpChildList[j]);
//                   std::cout << " min > tmp " << status << std::endl;
// 		  status=(tmpChildList[j]<max);
//                   std::cout << " tmp < max " << status << std::endl;
//                 }

                  if (tmpChildList[j] != max) {
                    std::cout << "  - Invalid Case" << std::endl;

                    std::cout << "   - tmp:" << tmpChildList[j] << std::endl;
                    std::cout << "   - max:" << max << std::endl;
                    TreeNode dld = tmpChildList[j].getDLD();
                    std::cout << "   - DLD of the tmpChild:" << dld << std::endl;
                  }
                  // std::cout << (tmpChildList[j] < max) << std::endl;
                  // std::cout << (max < dld) << std::endl;

                  //max.printTreeNode();

                  assert(tmpChildList[j] == max);
                  break;
                }
              } //end for j
            }
            while (repeatLoop);
            break;
          } //end select between bottom-up/ top-down
        } //end for i
        currentNode = parentOfCurrent;
      } //end while
    } //end if special case

    if (includeMax) {
      newNodes.push_back(max);
    }

    // std::cout << "leaving " << __func__ << std::endl;
    PROF_COMPLETE_REGION_END
  } //end function

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

} //end namespace

