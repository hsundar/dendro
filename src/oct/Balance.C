
/**
  @file Balance.C
  @brief A set of functions for 2:1 balancing octrees. 
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  @author Hari Sundar, hsundar@gmail.com
  */

#include "parUtils.h"
#include "seqUtils.h"
#include "TreeNode.h"
#include "TreeNodePointer.h"
#include "testUtils.h"
#include "dendro.h"

#ifdef __DEBUG__
#ifndef __DEBUG_OCT__
#define __DEBUG_OCT__
#endif
#endif

#ifdef __DEBUG_OCT__
#ifndef __MEASURE_BAL_COMM__
#define __MEASURE_BAL_COMM__
#endif
#endif

namespace ot {

  //Assumption: in is globally sorted, linear and complete.
  int balanceOctree (std::vector<TreeNode > &in, std::vector<TreeNode > &out,
      unsigned int dim, unsigned int maxDepth, bool incCorner,
      MPI_Comm comm, MPI_Comm* newCommPtr, bool* iAmActive) {
#ifdef __PROF_WITH_BARRIER__
    MPI_Barrier(comm);
#endif
    PROF_BAL_BEGIN

      int rank, size;
    MPI_Comm_size(comm,&size);
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
      std::cout<<"Balance Octree. inpSize: "<<in.size()<<" activeNpes: 1"<<std::endl; 
#endif

      out = in;
      in.clear();
      comboRipple(out, incCorner);
      PROF_BAL_END
    }

    MPI_Comm_rank(comm,&rank);

    //Check if the input is too small...
    DendroIntL locSize = in.size();
    DendroIntL globSize;
    PROF_BAL_COMM_BEGIN

      par::Mpi_Allreduce<DendroIntL>(&locSize, &globSize, 1, MPI_SUM, comm);

    PROF_BAL_COMM_END

      //min grain size = 1000
      const DendroIntL THOUSAND = 1;
    if (globSize < (THOUSAND*size)) {
      int splittingSize = (globSize/THOUSAND); 
      if(splittingSize == 0) {
        splittingSize = 1; 
      }

      unsigned int avgLoad = (globSize/splittingSize);
      int leftOvers = (globSize - (splittingSize*avgLoad));

      PROF_BAL_SCATTER_BEGIN

        std::vector<TreeNode> tmpIn;
      if(rank >= splittingSize) {
        par::scatterValues<ot::TreeNode>(in, tmpIn, 0, comm);
      }else if(rank < leftOvers) {
        par::scatterValues<ot::TreeNode>(in, tmpIn, (avgLoad+1), comm);
      }else {
        par::scatterValues<ot::TreeNode>(in, tmpIn, avgLoad, comm);
      }

      PROF_BAL_SCATTER_END

        in.clear();

      PROF_BAL_SPLIT_COMM_BEGIN

        MPI_Comm newComm;
      par::splitCommUsingSplittingRank(splittingSize, &newComm, comm);

      PROF_BAL_SPLIT_COMM_END

#ifndef __SILENT_MODE__
        if(!rank) {
          std::cout<<"Input to Balance is small ("<<globSize
            <<"). npes = "<<size<<" Splitting Comm. "<<std::endl;
        }
#endif

      if(rank < splittingSize) {
        balanceOctree (tmpIn, out, dim, maxDepth, incCorner, newComm, NULL, NULL);
      } else {
        if(iAmActive) {
          *iAmActive = false;
        }
      }
      tmpIn.clear();

      if(newCommPtr) {
        *newCommPtr = newComm;
      }

      PROF_BAL_END
    }//end if reduce procs

#ifndef __SILENT_MODE__
    if(!rank) {
      std::cout<<"Balance Octree. inpSize: "<<globSize <<" activeNpes: "<<size<<std::endl; 
    }
#endif

    //////////////////////////////////////////////////////////////////////////////////////////////////
    //Partition in and create blocks (blocks must be globally sorted).
    //Prepare for blockPart.

    std::vector<ot::TreeNode> blocks;    
    std::vector<ot::TreeNode> minsAllBlocks;

#ifdef __PROF_WITH_BARRIER__
    MPI_Barrier(comm);
#endif
    PROF_BAL_BPART1_BEGIN

      blockPartStage1(in, blocks, dim, maxDepth, comm);            

    PROF_BAL_BPART1_END

#ifdef __PROF_WITH_BARRIER__
      MPI_Barrier(comm);
#endif

    PROF_BAL_BPART2_BEGIN

      blockPartStage2(in, blocks, minsAllBlocks, dim, maxDepth, comm);

    PROF_BAL_BPART2_END

      //blocks will be sorted.

      assert(!blocks.empty());
    TreeNode myFirstBlock = blocks[0];
    TreeNode myLastBlock = blocks[blocks.size()-1];

#ifdef __DEBUG_OCT__
    assert(areComparable(myFirstBlock, in[0]));
    assert(myFirstBlock.isAncestor(in[0]) || myFirstBlock == in[0]);
#endif

#ifdef __USE_AGV_FOR_BAL__
    // 1. locally assign wgts to all blocks. The wgts are set to the rank of
    // the processor so that the owner of the block can be easily identified.
    //We don't need to reset the weights of blocks later since weights are not
    //used anywhere else. All comparisons are based only on the Morton
    //ordering.
    for(unsigned int i = 0; i < blocks.size(); i++) {
      blocks[i].setWeight(rank);
    }

    // 2. Communicate the blocks to all processors.
    int totalSize=0;
    int *blockSizes = new int [size];
    int *blockDisps = new int [size];

    int localBlockSize = blocks.size();
    PROF_BAL_COMM_BEGIN

      par::Mpi_Allgather<int>(&localBlockSize, blockSizes, 1, comm);

    PROF_BAL_COMM_END

      totalSize += blockSizes[0];
    blockDisps[0] = 0;
    for (unsigned int i = 1; i < size; i++) {
      totalSize += blockSizes[i];
      blockDisps[i] = blockDisps[i-1] + blockSizes[i-1];
    }//end for i

    // allocate for the blocks ...
    std::vector<TreeNode> allBlocks(totalSize);

    ot::TreeNode* blocksPtr = NULL;
    ot::TreeNode* allBlocksPtr = NULL;
    if(!blocks.empty()) {
      blocksPtr = &(*(blocks.begin()));
    }
    if(!allBlocks.empty()) {
      allBlocksPtr = &(*(allBlocks.begin()));
    }
    PROF_BAL_COMM_BEGIN

      par::Mpi_Allgatherv<ot::TreeNode>( blocksPtr, blockSizes[rank],
          allBlocksPtr, blockSizes, blockDisps, comm );

    PROF_BAL_COMM_END

#ifndef __SILENT_MODE__
      if(!rank) {
        std::cout<<"# AllBlocks: "<<allBlocks.size()<<std::endl;
        for(int i = 0; i < size; i++) {
          std::cout<<"blockSizes["<<i<<"] = "<<blockSizes[i]<<std::endl;
        }
        for(int i = 0; i < allBlocks.size(); i++) {
          std::cout<<"allBlocks["<<i<<"].Wt() = "<<(allBlocks[i].getWeight())<<std::endl;
        }
      }
#endif

#endif

    //////////////////////////////////////////////////////////////////////////////////////////////////
    //Local Balancing

    //Block Balance
    std::vector<ot::TreeNode> allBoundaryLeaves;
    std::vector<unsigned int> maxBlockBndVec;

    balanceBlocks(in, blocks, out, allBoundaryLeaves, incCorner, &maxBlockBndVec);
    in.clear();

#ifdef __USE_AGV_FOR_BAL__
    std::vector<TreeNode> myNhBlocks;
    selectNeighboringBlocks(allBlocks, blocks, maxBlockBndVec, rank, myNhBlocks);

    allBlocks.clear();

#ifdef __MEASURE_BAL_COMM__
    MPI_Barrier(comm);
    unsigned int localMyNhBlocksSize = myNhBlocks.size();
    unsigned int* globalMyNhBlocksSize = new unsigned int[size];
    par::Mpi_Gather<unsigned int>(&localMyNhBlocksSize, globalMyNhBlocksSize, 1, 0, comm);
    if(!rank) {
      for(int i = 0; i < size; i++) {
        std::cout<<"globalMyNhBlocksSize["<<i<<"] = "<<globalMyNhBlocksSize[i]<<std::endl;
      }//end for i
    }
    delete [] globalMyNhBlocksSize;
    MPI_Barrier(comm);
#endif
#endif

    blocks.clear();
    maxBlockBndVec.clear();

#ifdef __USE_AGV_FOR_BAL__
    if(!myNhBlocks.empty()) {
      //Merge (result must remain sorted and linear) our local
      //boundary leaves and myNhBlocks. We need to perform the
      //intra-processor local balance on this
      //combined list. Note, this step is slightly different from 
      //what was described in the 2:1 balancing paper. In the algorithm, descibed
      //in the paper the blocks from other processors are not taken into account
      //for this step. By taking these blocks into account as well, we expect
      //to identify some very large octants early on and hence reduce the
      //insulation layer sizes before the actual 2-stage inter-processor step
      std::vector<ot::TreeNode> boundaryAndOtherBlocks;

      //Three Mutually Exclusive Cases:
      if(myNhBlocks[myNhBlocks.size()-1] < myFirstBlock) {
        //1) All myNhBlocks are < myFirstBlock

        boundaryAndOtherBlocks.insert(boundaryAndOtherBlocks.end(),
            myNhBlocks.begin(), myNhBlocks.end());

        boundaryAndOtherBlocks.insert(boundaryAndOtherBlocks.end(),
            allBoundaryLeaves.begin(), allBoundaryLeaves.end());

      } else if(myNhBlocks[0] > myLastBlock) {
        //2) All myNhBlocks are > myLastBlock

        boundaryAndOtherBlocks.insert(boundaryAndOtherBlocks.end(),
            allBoundaryLeaves.begin(), allBoundaryLeaves.end());

        boundaryAndOtherBlocks.insert(boundaryAndOtherBlocks.end(),
            myNhBlocks.begin(), myNhBlocks.end());

      } else {
        //3 some myNhBlocks are < myFirstBlock and some are > myLastBlock

        unsigned int myStIdxInNhBlocks = 0;
        unsigned int myEndIdxInNhBlocks = 0;

        for(unsigned int i = 0; i < myNhBlocks.size(); i++) {
          if(myNhBlocks[i] > myFirstBlock) {
            myStIdxInNhBlocks = i;
            break;
          }
        }

        //We don't need to use DLD here for comparison since there can not be
        //any overlap between our blocks/octants and myNhBlocks 
        for(unsigned int i = myNhBlocks.size(); i > 0; i--) {
          if(myNhBlocks[i-1] < myLastBlock) {
            myEndIdxInNhBlocks = i;
            break;
          }
        }

        boundaryAndOtherBlocks.insert(boundaryAndOtherBlocks.begin(),
            myNhBlocks.begin(), (myNhBlocks.begin() + myStIdxInNhBlocks));

        boundaryAndOtherBlocks.insert(boundaryAndOtherBlocks.end(),
            allBoundaryLeaves.begin(), allBoundaryLeaves.end());

        boundaryAndOtherBlocks.insert(boundaryAndOtherBlocks.end(),
            (myNhBlocks.begin() + myEndIdxInNhBlocks), myNhBlocks.end());

      }

      allBoundaryLeaves.clear();

      //Intra-processor Balance
      ripple(boundaryAndOtherBlocks, incCorner);

      //Extract the local boundary elements
      unsigned int myStartIdx = 0;
      unsigned int myEndIdx = 0;

      bool foundFirst = false;
      bool foundLast = false;
      for(unsigned int i = 0; i < boundaryAndOtherBlocks.size(); i++) {
        if(boundaryAndOtherBlocks[i] >= myFirstBlock) {
          myStartIdx = i;
          foundFirst = true;
          break;
        }
      }

      //We need DLD here since we need to account for descendants of
      //myLastBlock as well
      for(unsigned int i = boundaryAndOtherBlocks.size(); i > 0; i--) {
        if(boundaryAndOtherBlocks[i-1] <= myLastBlock.getDLD()) {
          myEndIdx = i;
          foundLast = true;
          break;
        }
      }

#ifdef __DEBUG_OCT__
      //There must be at least one element belonging to our processor 
      assert(foundFirst && foundLast);
#endif

      allBoundaryLeaves.insert(allBoundaryLeaves.begin(),
          boundaryAndOtherBlocks.begin() + myStartIdx, 
          boundaryAndOtherBlocks.begin() + myEndIdx);

      boundaryAndOtherBlocks.clear();

    } else {
      //Intra-processor Balance
      ripple(allBoundaryLeaves, incCorner);
    }//end if myNhBlocks is empty
    myNhBlocks.clear();
#else
    //Intra-processor Balance
    ripple(allBoundaryLeaves, incCorner);
#endif

    mergeComboBalAndPickBoundary(out, allBoundaryLeaves, myFirstBlock, myLastBlock);

#ifdef __MEASURE_BAL_COMM__
    MPI_Barrier(comm);
    unsigned int localAllBndSize = allBoundaryLeaves.size();
    unsigned int* globalAllBndSize = new unsigned int[size];
    par::Mpi_Gather<unsigned int>(&localAllBndSize, globalAllBndSize, 1, 0, comm);
    if(!rank) {
      for(int i = 0; i < size; i++) {
        std::cout<<" globalAllBndSize["<<i<<"]= "<<globalAllBndSize[i]<<std::endl;
      }//end for i
    }
    delete [] globalAllBndSize;
    MPI_Barrier(comm);
#endif

    //There is no need to sort, the lists are already sorted.

    //////////////////////////////////////////////////////////////////////////////////////////////////
    //Preparation for inter-processor balance....
    //First Stage of Communication...
    std::vector<TreeNode> * sendNodes = new std::vector<TreeNode> [size];
    std::vector<unsigned int> * sentToPid = NULL;

    std::vector<ot::TreeNode> sendK;
    int *sendCnt = new int[size];
    int *recvCnt = new int[size];
    int *sendOffsets = new int[size];

    // 3. Loop through all local boundary nodes and determine which processors
    // need to be aware of those nodes. Create lists that shall be sent.

    for (int i = 0; i < size; i++) {
      sendCnt[i] = 0;
      sendNodes[i].clear();
    }

    if(!allBoundaryLeaves.empty()) {
      sentToPid = new std::vector<unsigned int> [allBoundaryLeaves.size()];
    } else {
      sentToPid = NULL;
    }

    //create sendNodes
    prepareBalComm1MessagesType2(allBoundaryLeaves, minsAllBlocks, rank, dim, maxDepth,
        sendNodes, sentToPid, sendCnt);

    // 4. Actual send/recv. to exchange nodes.
    //
    // 4a.
    // Now do an All2All to get numKeysRecv
    PROF_BAL_COMM_BEGIN

      par::Mpi_Alltoall<int>(sendCnt, recvCnt, 1, comm);

    PROF_BAL_COMM_END

#ifdef __MEASURE_BAL_COMM__
      unsigned int numProcsSend1 = 0;
    unsigned int numProcsRecv1 = 0;
#endif

    // 4b. Concatenate all nodes into one single Carray ...
    unsigned int totalSend = 0;
    unsigned int totalRecv = 0;
    for (unsigned int i = 0; i < size; i++) {
      totalSend += sendCnt[i];
      totalRecv += recvCnt[i];
#ifdef __MEASURE_BAL_COMM__
      if(sendCnt[i]) {
        numProcsSend1++;
      }
      if(recvCnt[i]) {
        numProcsRecv1++;
      }
#endif
    }//end for i

    // create the send and recv buffers ...

    int *recvOffsets1 = new int[size];

    // Now create sendK
    sendOffsets[0] = 0;
    recvOffsets1[0] = 0;

    // compute offsets ...
    for (int i = 1; i < size; i++) {
      sendOffsets[i] = sendOffsets[i-1] + sendCnt[i-1];
      recvOffsets1[i] = recvOffsets1[i-1] + recvCnt[i-1];
    }//end for i

    sendK.resize(totalSend);
    for (int i = 0; i < size; i++) {
#ifdef __DEBUG_OCT__
      assert(seq::test::isUniqueAndSorted(sendNodes[i]));
#endif
      for (unsigned int j = 0; j < sendCnt[i]; j++) {
        sendK[sendOffsets[i] + j] = sendNodes[i][j];
      }//end for j
    }//end for i

    for (unsigned int i = 0; i < size; i++) {
      sendNodes[i].clear();
    }

    // 4c. Perform SendRecv to send and receive all keys ...
    std::vector<ot::TreeNode> recvK1(totalRecv);

    ot::TreeNode* sendKptr = NULL;
    ot::TreeNode* recvK1ptr = NULL;
    if(!sendK.empty()) {
      sendKptr = &(*(sendK.begin()));
    }
    if(!recvK1.empty()) {
      recvK1ptr = &(*(recvK1.begin()));
    }

    PROF_BAL_COMM_BEGIN

      par::Mpi_Alltoallv_sparse<ot::TreeNode>(sendKptr, sendCnt, sendOffsets,        
          recvK1ptr, recvCnt, recvOffsets1, comm);

    PROF_BAL_COMM_END

      sendK.clear();

#ifdef __MEASURE_BAL_COMM__
    MPI_Barrier(comm);
    unsigned int* allTotalSend1 = new unsigned int[size];
    unsigned int* allTotalRecv1 = new unsigned int[size];
    unsigned int* allNumProcsSend1 = new unsigned int[size];
    unsigned int* allNumProcsRecv1 = new unsigned int[size];
    par::Mpi_Gather<unsigned int>(&totalSend, allTotalSend1, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&totalRecv, allTotalRecv1, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&numProcsSend1, allNumProcsSend1, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&numProcsRecv1, allNumProcsRecv1, 1, 0, comm);
    if(!rank) {
      for(int i = 0; i < size; i++) {
        std::cout<<"allTotalSend1["<<i<<"] in bal: "<<allTotalSend1[i]<<std::endl; 
        std::cout<<"allTotalRecv1["<<i<<"] in bal: "<<allTotalRecv1[i]<<std::endl; 
        std::cout<<"allNumProcsSend1["<<i<<"] in bal: "<<allNumProcsSend1[i]<<std::endl; 
        std::cout<<"allNumProcsRecv1["<<i<<"] in bal: "<<allNumProcsRecv1[i]<<std::endl; 
      }
    }
    delete [] allTotalSend1;
    delete [] allTotalRecv1;
    delete [] allNumProcsSend1;
    delete [] allNumProcsRecv1;
    MPI_Barrier(comm);
#endif

    //////////////////////////////////////////////////////////////////////////////////////////////////
    //Second Stage of Communication...
    // 5. Loop through all local boundary nodes and determine which processors
    // need to be aware of those nodes. Create lists that shall be sent.

    std::vector<ot::TreeNode> wList;
    std::vector<std::vector<unsigned int> > wListRanks; 

    prepareWlistInBal(recvK1, recvCnt, size, myFirstBlock, myLastBlock, wList, wListRanks);

    for (int i = 0; i < size; i++) {
      sendCnt[i] = 0;
      sendNodes[i].clear();
    }

    prepareBalComm2Messages(allBoundaryLeaves, wList, wListRanks, 
        sendNodes, sentToPid, sendCnt);

    for (int i = 0; i < allBoundaryLeaves.size(); i++) {
      sentToPid[i].clear();
    }

    if(sentToPid) {
      delete [] sentToPid;
      sentToPid = NULL;
    }

    // 6. Actual send/recv. to exchange nodes.
    //
    // 6a.
    // Now do an All2All to get numKeysRecv
    PROF_BAL_COMM_BEGIN

      par::Mpi_Alltoall<int>(sendCnt, recvCnt, 1, comm);

    PROF_BAL_COMM_END

#ifdef __MEASURE_BAL_COMM__
      unsigned int numProcsSend2 = 0;
    unsigned int numProcsRecv2 = 0;
#endif
    // 6b. Concatenate all nodes into one single Carray ...
    totalSend = 0;
    totalRecv = 0;
    for (unsigned int i = 0; i < size; i++) {
      totalSend+= sendCnt[i];
      totalRecv+= recvCnt[i];
#ifdef __MEASURE_BAL_COMM__
      if(sendCnt[i]) {
        numProcsSend2++;
      }
      if(recvCnt[i]) {
        numProcsRecv2++;
      }
#endif
    }//end for i

    // create the send and recv buffers ...

    int *recvOffsets2 = new int[size];
    // Now create sendK
    sendOffsets[0] = 0;
    recvOffsets2[0] = 0;

    // compute offsets ...
    for (int i=1; i<size; i++) {
      sendOffsets[i] = sendOffsets[i-1] + sendCnt[i-1];
      recvOffsets2[i] = recvOffsets2[i-1] + recvCnt[i-1];
    }//end for i

    sendK.resize(totalSend);
    for (int i=0; i<size; i++) {
#ifdef __DEBUG_OCT__
      assert(seq::test::isUniqueAndSorted(sendNodes[i]));
#endif
      for (unsigned int j = 0; j < sendCnt[i]; j++) {
        sendK[sendOffsets[i] + j] = sendNodes[i][j];
      }//end for j
    }//end for i

    for (unsigned int i = 0; i < size; i++) {
      sendNodes[i].clear();
    }

    delete [] sendNodes;
    sendNodes = NULL;

    std::vector<ot::TreeNode> recvK2(totalRecv);

    sendKptr = NULL;
    ot::TreeNode* recvK2ptr = NULL;
    if(!sendK.empty()) {
      sendKptr = &(*(sendK.begin()));
    }
    if(!recvK2.empty()) {
      recvK2ptr = &(*(recvK2.begin()));
    }

    PROF_BAL_COMM_BEGIN

      par::Mpi_Alltoallv_sparse<ot::TreeNode>(sendKptr, sendCnt, sendOffsets,        
          recvK2ptr, recvCnt, recvOffsets2, comm);

    PROF_BAL_COMM_END

      sendK.clear();

    delete [] sendCnt;
    sendCnt = NULL;

    delete [] recvCnt;
    recvCnt = NULL;

    delete [] sendOffsets;
    sendOffsets = NULL;


#ifdef __MEASURE_BAL_COMM__
    MPI_Barrier(comm);
    unsigned int* allTotalSend2 = new unsigned int[size];
    unsigned int* allTotalRecv2 = new unsigned int[size];
    unsigned int* allNumProcsSend2 = new unsigned int[size];
    unsigned int* allNumProcsRecv2 = new unsigned int[size];
    par::Mpi_Gather<unsigned int>(&totalSend, allTotalSend2, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&totalRecv, allTotalRecv2, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&numProcsSend2, allNumProcsSend2, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&numProcsRecv2, allNumProcsRecv2, 1, 0, comm);
    if(!rank) {
      for(int i = 0; i < size; i++) {
        std::cout<<"allTotalSend2["<<i<<"] in bal: "<<allTotalSend2[i]<<std::endl; 
        std::cout<<"allTotalRecv2["<<i<<"] in bal: "<<allTotalRecv2[i]<<std::endl; 
        std::cout<<"allNumProcsSend2["<<i<<"] in bal: "<<allNumProcsSend2[i]<<std::endl; 
        std::cout<<"allNumProcsRecv2["<<i<<"] in bal: "<<allNumProcsRecv2[i]<<std::endl; 
      }
    }
    delete [] allTotalSend2;
    delete [] allTotalRecv2;
    delete [] allNumProcsSend2;
    delete [] allNumProcsRecv2;
    MPI_Barrier(comm);
#endif

    //////////////////////////////////////////////////////////////////////////////////////////////////

    //7.a.Merge (In-place) recieved keys from both stages of communication....
    std::vector<ot::TreeNode> recvK;

    unsigned int myOff1 = recvOffsets1[rank];
    unsigned int myOff2 = recvOffsets2[rank];
    unsigned int myOff = myOff1 + myOff2;

    mergeRecvKeysInBal(recvK1, recvOffsets1, recvK2, recvOffsets2, size, recvK);

    delete [] recvOffsets1;
    recvOffsets1 = NULL;

    recvK1.clear();

    delete [] recvOffsets2;
    recvOffsets2 = NULL;

    recvK2.clear();

    //7.b.Merge (In-place) recieved octants with local octants....
    unsigned int n = static_cast<unsigned int>(allBoundaryLeaves.size());
    std::vector <TreeNode > nodes( recvK.size() + n );

    //Note this is an inplace insertion. This is done to preserve the sorted
    //order. Each recvK[i] is independently sorted. Also, since the intial
    //blocks were globally sorted and since the elements on each processor are
    //only decendants of these blocks, hence recvK[i][j] < recvK[i+][k] for all i,j
    //and k.
    for (int i = 0; i < myOff; i++) {
      nodes[i] = recvK[i];
    }
    for (int i = myOff; i < myOff + n; i++) {
      nodes[i] = allBoundaryLeaves[i - myOff];
    }
    for (int i = myOff + n; i < nodes.size(); i++) {
      nodes[i] = recvK[i - n];
    }
    recvK.clear();
    allBoundaryLeaves.clear();

    //////////////////////////////////////////////////////////////////////////////////////////////////

    // 8. Now perform seq. balance on the new list of local nodes.
    ripple(nodes, incCorner);

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // 9. Discard all nodes that do not belong to the original domain. This can
    // easily be obtained from the mins of blocks. 
    unsigned int myStIdxInNodes = 0;
    unsigned int myEndIdxInNodes = 0;
    bool foundStart = false;
    bool foundEnd = false;
    for(unsigned int i = 0; i < nodes.size(); i++) {
      if(nodes[i] >= myFirstBlock) {
        myStIdxInNodes = i;
        foundStart = true;
        break;
      }
    }

    //We need DLD here since we need to account for descendants of
    //myLastBlock as well
    for(unsigned int i = nodes.size(); i > 0; i--) {
      if(nodes[i-1] <= myLastBlock.getDLD()) {
        myEndIdxInNodes = i;
        foundEnd = true;
        break;
      }
    }

#ifdef __DEBUG_OCT__
    assert(foundStart && foundEnd);
#endif

    allBoundaryLeaves.insert(allBoundaryLeaves.begin(),
        (nodes.begin() + myStIdxInNodes), (nodes.begin() + myEndIdxInNodes));

    nodes.clear();

    //10.Merge (In-place) results from all three stages....
    finalMergeInBal(out, allBoundaryLeaves);

    PROF_BAL_END
  }//end function

  //Assumption: in is sorted, linear and the morton space between the first and
  //last elements is complete.
  int comboRipple(std::vector<TreeNode> & in, bool incCorner, const unsigned int maxNum) {
    PROF_COMBO_RIPPLE_BEGIN

      if (in.size() < 2) {
        PROF_COMBO_RIPPLE_END
      }

    std::vector<ot::TreeNode> blocks;
    std::vector<ot::TreeNode> out;
    std::vector<ot::TreeNode> allBoundaryLeaves;

    TreeNode nca = getNCA(in[0],in[in.size()-1]);
    nca.addChildren(blocks);

    assert(maxNum > 0);

    if (in.size() > maxNum) {
      std::vector<TreeNode> *splitInp = NULL;

      if(!blocks.empty()) {
        splitInp = new std::vector<TreeNode> [blocks.size()];
      }

      unsigned int nextPt = 0;
      unsigned int nextNode = 0;
      //All elements of inp are inside some element in blocks.
      while (nextPt < in.size()) {
        //The first pt must be inside some block.
#ifdef __DEBUG_OCT__
        assert(nextNode < blocks.size());
        assert(areComparable(blocks[nextNode], in[nextPt]));
#endif
        if ((blocks[nextNode].isAncestor(in[nextPt])) ||
            (blocks[nextNode] == in[nextPt])) {
          splitInp[nextNode].push_back(in[nextPt]);
          nextPt++;
        } else {
          nextNode++;
          if (nextNode == blocks.size()) {
            assert(false);
          }
        }
      }//end while
      for (unsigned int i = 0; i < blocks.size(); i++) {
        comboRipple(splitInp[i], incCorner, maxNum);
      }//end for i

      unsigned int allBoundarySz = 0;
      unsigned int blockOutSize = 0;

      std::vector<TreeNode> *blockBoundaries = NULL;

      if(!blocks.empty()) {
        blockBoundaries = new std::vector<TreeNode>[blocks.size()];
      }

      unsigned int numEmptyBlocks=0;
      for (unsigned int bi = 0; bi < blocks.size(); bi++) {
        if (splitInp[bi].empty()) {
          continue; 
        }
        blockOutSize += splitInp[bi].size();
        if (splitInp[bi].size() > 1) {
          blocks[bi].pickInternalBoundaryCells(splitInp[bi], blockBoundaries[bi]);
        }
        if (blockBoundaries[bi].empty()) {
          //the case where splitInp[bi] = blocks[bi].
          blockBoundaries[bi].push_back(blocks[bi]);  
          numEmptyBlocks++;
        }
        allBoundarySz += blockBoundaries[bi].size();
      }//end for bi

      //Concatenate the lists into out and allBoundaryLeaves
      out.resize(blockOutSize);
      allBoundaryLeaves.resize(allBoundarySz);
      unsigned int nodeCtr = 0;
      unsigned int boundaryCtr = 0;
      for (unsigned int bi = 0; bi < blocks.size(); bi++) {
        for (unsigned int bj = 0; bj < splitInp[bi].size(); bj++) {
          out[nodeCtr] = splitInp[bi][bj];
          nodeCtr++;
        }
        for (unsigned int bj = 0; bj < blockBoundaries[bi].size(); bj++) {
          allBoundaryLeaves[boundaryCtr] = blockBoundaries[bi][bj];
          boundaryCtr++;
        }
      }//end for bi
      delete [] splitInp;
      splitInp = NULL;

      delete [] blockBoundaries;
      blockBoundaries = NULL;
    } else {
      balanceBlocks (in, blocks, out, allBoundaryLeaves, incCorner); 
    }

    in.clear();

    //Inter-Block Balance
    ripple(allBoundaryLeaves, incCorner);

    //Merge (in-place) results from the two stages....
    in.resize(out.size() + allBoundaryLeaves.size());
    unsigned int tmpLsz=0;
    unsigned int bndCnt=0;
    for (unsigned int i =0;i<out.size();i++) {
      if ( bndCnt < allBoundaryLeaves.size() ) {
        if ( out[i] == allBoundaryLeaves[bndCnt] ) {
          in[tmpLsz++] = out[i];
          bndCnt++;
        } else if (out[i] < allBoundaryLeaves[bndCnt] ) {
#ifdef __DEBUG_OCT__
          assert(areComparable(out[i], allBoundaryLeaves[bndCnt]));
#endif
          if (out[i].isAncestor(allBoundaryLeaves[bndCnt]) ) {
            while ( (bndCnt < allBoundaryLeaves.size()) && 
                out[i].isAncestor(allBoundaryLeaves[bndCnt]) ) {
              in[tmpLsz++] = allBoundaryLeaves[bndCnt++];
#ifdef __DEBUG_OCT__
              if(bndCnt < allBoundaryLeaves.size()) {
                assert(areComparable(out[i], allBoundaryLeaves[bndCnt]));
              }
#endif
            }
          } else {
            in[tmpLsz++] = out[i];
          }
        } else {
          // nodes[i] > allBdy .. so insert 
          in[tmpLsz++] = allBoundaryLeaves[bndCnt++];
        }
      } else {
        in[tmpLsz++] = out[i];
      }
    }//end for i

    in.resize(tmpLsz);
    out.clear();

    //There is no need to sort, the lists are already sorted.
    PROF_COMBO_RIPPLE_END
  }//end function

  //Assumption: All elements of inp are inside some unique element in blocks.
  int balanceBlocks (const std::vector<TreeNode> &inp,
      const std::vector<TreeNode> &blocks, std::vector<TreeNode> &nodes,
      std::vector<TreeNode> &allBoundaryLeaves, bool incCorner, 
      std::vector<unsigned int> *maxBlockBndVec ) {
    PROF_CON_BAL_BEGIN

      if (inp.empty()) {
        nodes.clear();
        allBoundaryLeaves.clear();
        if (maxBlockBndVec != NULL) {
          maxBlockBndVec->clear();
        }
        PROF_CON_BAL_END
      }

    unsigned int nextPt;
    unsigned int nextNode;

    assert(!blocks.empty());
    std::vector<TreeNode> *blockOut = new std::vector<TreeNode> [blocks.size()];
    std::vector<TreeNode> *splitInp = new std::vector<TreeNode> [blocks.size()];

    nextPt =0;
    nextNode=0;
    //All elements of inp are inside some element in blocks.
    while (nextPt<inp.size()) {
      //The first pt must be inside some block.
#ifdef __DEBUG_OCT__
      assert(areComparable(blocks[nextNode], inp[nextPt]));
#endif
      if ((blocks[nextNode].isAncestor(inp[nextPt])) ||
          (blocks[nextNode] == inp[nextPt])) {
        splitInp[nextNode].push_back(inp[nextPt]);
        nextPt++;
      } else {
        nextNode++;
        if (nextNode == blocks.size()) {
          assert(false);
        }
      }
    }//end while

    //Create Local Trees
    for (unsigned int bi=0;bi<blocks.size();bi++) {
      //This also sorts and makes the vector unique inside.
      blocks[bi].balanceSubtree(splitInp[bi],blockOut[bi],incCorner,true);

      splitInp[bi].clear();
      //This tackles the case where blocks[bi] has no decendants.
      if (blockOut[bi].empty()) {
        blockOut[bi].push_back(blocks[bi]);
      }
    }//end for bi

    delete [] splitInp; 
    splitInp = NULL;

    unsigned int allBoundarySz = 0;
    unsigned int blockOutSize = 0;

    std::vector<TreeNode> *blockBoundaries = 
      new std::vector<TreeNode>[blocks.size()];

    unsigned int numEmptyBlocks=0;
    for (unsigned int bi = 0; bi < blocks.size(); bi++) {
      blockOutSize += blockOut[bi].size();
      if (blockOut[bi].size() > 1) {
        blocks[bi].pickInternalBoundaryCells(blockOut[bi],blockBoundaries[bi]);
      }
      if (blockBoundaries[bi].empty()) {
        //the case where blockOut[bi] = blocks[bi].
        blockBoundaries[bi].push_back(blocks[bi]);  
        numEmptyBlocks++;
      }
      allBoundarySz += blockBoundaries[bi].size();
    }//end for bi

    //Concatenate the lists into nodes and allBoundaryLeaves
    nodes.resize(blockOutSize);
    allBoundaryLeaves.resize(allBoundarySz);
    unsigned int nodeCtr = 0;
    unsigned int boundaryCtr = 0;

    if (maxBlockBndVec != NULL) {
      maxBlockBndVec->resize(blocks.size());
      unsigned int maxDepth;
      if (!blocks.empty()) {
        maxDepth = blocks[0].getMaxDepth();  
      }
      for (unsigned int bi = 0; bi < blocks.size(); bi++) {
        (*maxBlockBndVec)[bi] = 0;
        for (unsigned int bj = 0; bj < blockOut[bi].size(); bj++) {
          nodes[nodeCtr] = blockOut[bi][bj];
          nodeCtr++;
        }//end for bj
        for (unsigned int bj = 0; bj < blockBoundaries[bi].size(); bj++) {
          unsigned int len = (1u<<(maxDepth - (blockBoundaries[bi][bj].getLevel())));
          if (len > ((*maxBlockBndVec)[bi])) {
            (*maxBlockBndVec)[bi] = len;
          }
          allBoundaryLeaves[boundaryCtr] = blockBoundaries[bi][bj];
          boundaryCtr++;
        }//end for bj
      }//end for bi
    } else {
      for (unsigned int bi = 0; bi < blocks.size(); bi++) {
        for (unsigned int bj = 0; bj < blockOut[bi].size(); bj++) {
          nodes[nodeCtr] = blockOut[bi][bj];
          nodeCtr++;
        }//end for bj
        for (unsigned int bj = 0; bj < blockBoundaries[bi].size(); bj++) {
          allBoundaryLeaves[boundaryCtr] = blockBoundaries[bi][bj];
          boundaryCtr++;
        }//end for bj
      }//end for bi
    }

    delete [] blockOut;
    blockOut = NULL;

    delete [] blockBoundaries;
    blockBoundaries = NULL;

    PROF_CON_BAL_END
  }//end function

  /*
  //New implementation of the ripple algorithm. Implemented on Dec 22, 2007.
  int ripple(std::vector<TreeNode> & nodes, bool incCorners) {
  PROF_RIPPLE_BAL_BEGIN 

  if (!nodes.size()) {
  PROF_RIPPLE_BAL_END 
  }

  unsigned int dim = nodes[0].getDim();
  unsigned int maxDepth = nodes[0].getMaxDepth();
  TreeNode root(dim,maxDepth);

  unsigned int maxLev = 1;
  for (unsigned int i=0;i<nodes.size();i++) {
  if (nodes[i].getLevel() > maxLev) {
  maxLev = nodes[i].getLevel();
  }
  }//end for 

  for (unsigned int lev = maxLev; lev > 2; lev--) {
  unsigned int mLev = maxDepth;
  for (unsigned int i=0;i<nodes.size();i++) {
  if (nodes[i].getLevel() < mLev) {
  mLev = nodes[i].getLevel();
  }
  }//end for 

  if (mLev >= (lev-1)) {
//Difference between min and max levels is less than 2 .
break;
}

std::vector<TreeNode> wList;      
unsigned int wLen = 0;
wList.resize(nodes.size());
//Pick leaves at a particular level
for (int i =0;i<nodes.size();i++) {
if (nodes[i].getLevel() == lev) {
wList[wLen++] = nodes[i];
}
}//end for i
wList.resize(wLen);

std::vector<std::vector<TreeNode> > tList(wLen);
for (unsigned int i=0; i < wLen; i++) {
tList[i] = wList[i].getSearchKeys(incCorners);        
}//end for i
wList.clear();

std::vector<TreeNode> allKeys;
for (unsigned int i=0; i < wLen; i++) {
for(unsigned int j=0; j < tList[i].size(); j++) {
if(tList[i][j] > root) {
allKeys.push_back(tList[i][j]);
}
}
}    
tList.clear();

par::makeVectorUnique<ot::TreeNode>(allKeys, false);
unsigned int keyLen = allKeys.size();

std::vector<std::vector<TreeNode> > seedList(nodes.size());
unsigned int lastIdx = 0;
for (int i = 0; i < keyLen; i++) {
unsigned int idx;
//assumes nodes is sorted and unique.          
bool flag = par::maxLowerBound<TreeNode>(nodes, allKeys[i], idx, &lastIdx, NULL);
if (flag) {
lastIdx = idx;
if (nodes[idx].isAncestor(allKeys[i])) {
  if (lev > (nodes[idx].getLevel() + 1 )) {
    seedList[idx].push_back(allKeys[i].getAncestor(lev-1));                
  }//end if correct result
}
} //end if flag
}//end for i   
allKeys.clear();

//Include the new octants in the octree
//while preserving linearity and sorted order

//Seedlist may have duplicates. Seedlist 
// will be sorted.

std::vector<TreeNode> tmpList;
std::vector<std::vector<TreeNode> >allInternalLeaves(nodes.size());
unsigned int tmpSz = 0;
for (unsigned int i=0;i<nodes.size();i++) {                
  par::makeVectorUnique<TreeNode>(seedList[i],true);
  if (!seedList[i].empty()) {
    nodes[i].completeSubtree(seedList[i], allInternalLeaves[i]);
  }
  tmpSz += allInternalLeaves[i].size();
  if (allInternalLeaves[i].empty()) {
    tmpSz++;
  }
  seedList[i].clear();
}//end for i

seedList.clear();
tmpList.resize(tmpSz);
tmpSz=0;
for (unsigned int i=0;i<nodes.size();i++) {
  if (!allInternalLeaves[i].empty()) {
    for (int k=0;k<allInternalLeaves[i].size();k++) {
      tmpList[tmpSz++] = allInternalLeaves[i][k];
    }//end for k
    allInternalLeaves[i].clear();
  } else {
    tmpList[tmpSz++] = nodes[i];
  }
}//end for i
allInternalLeaves.clear();
nodes = tmpList;
tmpList.clear();
}//end for lev

PROF_RIPPLE_BAL_END 
}//end function
*/

//Original implementation of the ripple algorithm
int ripple(std::vector<TreeNode> & nodes, bool incCorners) {
  PROF_RIPPLE_BAL_BEGIN 

    if (!nodes.size()) {
      PROF_RIPPLE_BAL_END 
    }

  unsigned int dim = nodes[0].getDim();
  unsigned int maxDepth = nodes[0].getMaxDepth();
  TreeNode root(dim,maxDepth);

  unsigned int maxLev = 1;
  for (unsigned int i=0;i<nodes.size();i++) {
    if (nodes[i].getLevel() > maxLev) {
      maxLev = nodes[i].getLevel();
    }
  }//end for 

  for (unsigned int lev = maxLev; lev > 2; lev--) {
    unsigned int mLev = maxDepth;
    for (unsigned int i=0;i<nodes.size();i++) {
      if (nodes[i].getLevel() < mLev) {
        mLev = nodes[i].getLevel();
      }
    }//end for 

    if (mLev >= (lev-1)) {
      //Difference between min and max levels is less than 2 .
      break;
    }

    std::vector<TreeNode> wList;
    std::vector<std::vector<TreeNode> > seedList(nodes.size());
    unsigned int wLen = 0;
    wList.resize(nodes.size());

    //Pick leaves at a particular level
    for (int i =0;i<nodes.size();i++) {
      if (nodes[i].getLevel() == lev) {
        wList[wLen++] = nodes[i];
      }
    }//end for i
    wList.resize(wLen);

    for (unsigned int i=0; i < wList.size(); i++) {
      std::vector<TreeNode> tList = wList[i].getSearchKeys(incCorners);
      for (int j=0; j<tList.size(); j++) {
        unsigned int idx;
        //assumes nodes is sorted and unique.          
        bool flag = seq::maxLowerBound<TreeNode>(nodes, tList[j], idx,NULL,NULL);
        if (flag) {
#ifdef __DEBUG_OCT__
          assert(areComparable(nodes[idx], tList[j]));
#endif
          if (nodes[idx].isAncestor(tList[j])) {
            if (wList[i].getLevel() > (nodes[idx].getLevel() + 1 )) {
              nodes[idx].addBalancingDescendants(wList[i],
                  seedList[idx], incCorners);
            }//end if failing balance
          }//end if valid result
        } //end if flag
      }//end for j
      tList.clear();
    }//end for i

    wList.clear();
    std::vector<TreeNode> tmpList;
    std::vector<std::vector<TreeNode> >allInternalLeaves(nodes.size());
    unsigned int tmpSz = 0;
    for (unsigned int i=0;i<nodes.size();i++) {
      seq::makeVectorUnique<TreeNode>(seedList[i],false);
      if (!seedList[i].empty()) {
        if (seedList[i][0] == root) {
          seedList[i].erase(seedList[i].begin());
        }
      }
      if (!seedList[i].empty()) {
        nodes[i].completeSubtree(seedList[i], allInternalLeaves[i]);
      }
      tmpSz += allInternalLeaves[i].size();
      if (allInternalLeaves[i].empty()) {
        tmpSz++;
      }
      seedList[i].clear();
    }//end for i
    seedList.clear();
    tmpList.resize(tmpSz);
    tmpSz=0;
    for (unsigned int i=0;i<nodes.size();i++) {
      if (!allInternalLeaves[i].empty()) {
        for (int k=0;k<allInternalLeaves[i].size();k++) {
          tmpList[tmpSz++] = allInternalLeaves[i][k];
        }//end for k
        allInternalLeaves[i].clear();
      } else {
        tmpList[tmpSz++] = nodes[i];
      }
    }//end for i
    allInternalLeaves.clear();
    nodes = tmpList;
    tmpList.clear();
  }//end for lev

  PROF_RIPPLE_BAL_END 
}//end function


int pointerBasedRipple(std::vector<ot::TreeNode> & nodes, bool incCorners) {
  PROF_PTR_RIPPLE_BAL_BEGIN

    if (!nodes.size()) {
      PROF_PTR_RIPPLE_BAL_END 
    }

  unsigned int dim = nodes[0].getDim();
  unsigned int maxDepth = nodes[0].getMaxDepth();
  TreeNode root(dim,maxDepth);

  unsigned int maxLev = 1;
  for (unsigned int i = 0; i < nodes.size(); i++) {
    if (nodes[i].getLevel() > maxLev) {
      maxLev = nodes[i].getLevel();
    }
  }//end for i

  TreeNodePointer ptrOctree;

  convertLinearToPointer(nodes, ptrOctree);

  for(unsigned int lev = maxLev; lev > 2; lev--) {

    //Traverse the tree, select octants at level = lev and process them

    std::vector<ot::TreeNode> wList;

    appendOctantsAtLevel(ptrOctree, wList, lev);

    for(unsigned int i = 0; i < wList.size(); i++) {
      std::vector<TreeNode> tList = wList[i].getSearchKeys(incCorners);
      for(int j = 0; j < tList.size(); j++) {
        TreeNodePointer* searchResult = NULL;
        findOctantOrFinestAncestor(ptrOctree, tList[j], searchResult);
#ifdef __DEBUG_OCT__
        assert(searchResult);
        assert( ((searchResult->m_tnMe).isAncestor(tList[j]))
            || ((searchResult->m_tnMe) == tList[j]) );
        assert( (searchResult->m_tnpMyChildren) == NULL );
#endif
        //Check balance constraint
        if( ((searchResult->m_tnMe).getLevel()) < (lev - 1) ) {
          addOctantToTreeNodePointer((*searchResult),
              (tList[j].getAncestor((lev - 1))));
        }
      }//end for j
      tList.clear();
    }//end for i

    wList.clear();

  }//end for lev

  std::vector<ot::TreeNode> linOct;

  convertPointerToLinear(linOct, ptrOctree);

  //Select only those octants from linOct, which are either present in nodes or
  //are decendants of octants in nodes.
  //assumes nodes is sorted, linear and unique.          
  std::vector<ot::TreeNode> finalOct;
  unsigned int lCtr = 0;
  for(unsigned int i = 0; i < nodes.size(); i++) {
    while( (lCtr < linOct.size()) && (linOct[lCtr] < nodes[i]) ) {
      lCtr++;
    }
    if(lCtr < linOct.size()) {
      if(linOct[lCtr] == nodes[i]) {
        finalOct.push_back(nodes[i]);
      } else {
#ifdef __DEBUG_OCT__
        assert(nodes[i].isAncestor(linOct[lCtr]));
#endif
        while( (lCtr < linOct.size()) &&
            (nodes[i].isAncestor(linOct[lCtr])) ) {
          finalOct.push_back(linOct[lCtr]);
          lCtr++;
        }
      }
    } else {
      //Can't exhaust linOct before exhausting nodes
      assert(false);
    }
  }

  nodes = finalOct;

  finalOct.clear();
  linOct.clear();

  deleteTreeNodePointer(ptrOctree);

  PROF_PTR_RIPPLE_BAL_END
}//end function

int parallelRippleType3(std::vector<TreeNode> & nodes,
    bool incCorners, bool checkBailOut, bool rePart,
    unsigned int dim, unsigned int maxDepth, MPI_Comm comm)
{
  PROF_PAR_RIPPLE_TYPE3_BEGIN 

    TreeNode root(dim,maxDepth);

  unsigned int maxLev = 1;
  for(unsigned int i=0; i<nodes.size(); i++) {
    if (nodes[i].getLevel() > maxLev) {
      maxLev = nodes[i].getLevel();
    }
  }//end for 

  int rank,npes;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&npes);

#ifdef __DEBUG_OCT__
  MPI_Barrier(comm);	
  if(!rank) {
    std::cout<<"Computed local maxLev."<<std::endl;
  }
  MPI_Barrier(comm);	
#endif

  unsigned int globalMaxLev = maxLev;
  par::Mpi_Allreduce<unsigned int>(&maxLev, &globalMaxLev, 1, MPI_MAX, comm);

#ifdef __DEBUG_OCT__
  MPI_Barrier(comm);	
  if(!rank) {
    std::cout<<"Global maxLev: "<<globalMaxLev<<std::endl;
  }
  MPI_Barrier(comm);	
#endif

  for (unsigned int lev = globalMaxLev; lev > 2; lev--)
  {
#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Lev: "<<lev<<std::endl;
    }
    MPI_Barrier(comm);	
#endif
    if(checkBailOut) {
      unsigned int minLev = maxDepth;
      for (unsigned int i=0; i<nodes.size(); i++) {
        if (nodes[i].getLevel() < minLev) {
          minLev = nodes[i].getLevel();
        }
      }//end for i

      unsigned int globalMinLev = minLev;
      par::Mpi_Allreduce<unsigned int>(&minLev, &globalMinLev, 1, MPI_MIN, comm);

      if (globalMinLev >= (lev-1)) {
        //Difference between min and 
        //max levels is less than 2 .
        break;
      }
    }//end if check to bail out

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 1."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    std::vector<TreeNode> wList;
    unsigned int wLen = 0;
    wList.resize(nodes.size());
    //Pick leaves at a particular level
    for (int i=0; i<nodes.size(); i++) {
      if (nodes[i].getLevel() == lev) {
        wList[wLen] = nodes[i];
        wLen++;
      }
    }//end for i
    wList.resize(wLen);

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 2."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    std::vector<std::vector<TreeNode> >tList(wLen);
    for (unsigned int i=0; i < wLen; i++) {
      tList[i] = wList[i].getSearchKeys(incCorners);        
    }//end for i

    wList.clear();

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 3."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    std::vector<TreeNode> allKeys;
    for (unsigned int i=0; i < wLen; i++) {
      for(unsigned int j=0; j < tList[i].size(); j++) {
        if(tList[i][j] > root) {
          allKeys.push_back(tList[i][j]);
        }
      }
    }    
    tList.clear();

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 4."<<std::endl;
    }
    std::cout<<rank<<": allKeys.size(): "
      <<allKeys.size()<<std::endl;
    MPI_Barrier(comm);	
#endif

    //
    //Parallel searches and subsequent processing      
    //1. Compute the ranges controlled by each processor
    //2. Compare the keys with the ranges and 
    //send the keys to the appropriate processors.
    //3. Perform local searches using the keys
    //recieved. 
    //4. Compute the balancing descendants for the results
    //and create the 'seedList' vector. Note,
    //we process one level at a time. So all the keys
    //were generated by octants at the same level. The 
    //corresponding balancing descendants are simply the
    //ancestors of the keys at one level lower than this
    //level.
    //

    seq::makeVectorUnique<TreeNode>(allKeys,false);	

    unsigned int keyLen = allKeys.size();


    //First Get the mins from each processor.

    // allocate memory for the mins array
    std::vector<ot::TreeNode> mins (npes); 

    ot::TreeNode sendMin;
    if(!nodes.empty()) {
      sendMin = nodes[0]; //local min
    }else {
      sendMin = root;
    }

    par::Mpi_Allgather<ot::TreeNode>(&sendMin, &(*mins.begin()), 1, comm);    

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 5."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    //Now determine the processors which own these keys.

    unsigned int *partKeys = NULL;

    if(keyLen) {
      partKeys = new unsigned int[keyLen];    
    }

    for (unsigned int i=0; i<keyLen; i++) {
      unsigned int idx;
      //maxLB returns the last index in a 
      //sorted array such that a[ind] <= key 
      //and  a[index +1] > key
      bool found = seq::maxLowerBound<TreeNode >(mins,
          allKeys[i], idx, NULL, NULL);
      if (!found ) {
        //Can happen on incomplete domains
        partKeys[i] = rank;
      } else {
        partKeys[i] = idx;
      }
    }//end for i

    mins.clear();

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 6."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    int *numKeysSend = new int[npes];
    int *numKeysRecv = new int[npes];    

    for (int i=0; i<npes; i++) {
      numKeysSend[i] = 0;
    }

    // calculate the number of keys to send ...
    //and create the send buffer.
    std::vector<std::vector<unsigned int> >sendKtmp(npes);

    for (unsigned int i=0; i<keyLen; i++) {
      numKeysSend[partKeys[i]]++;
      sendKtmp[partKeys[i]].push_back(i);      
    }

    // Now do an All2All to get numKeysRecv

    par::Mpi_Alltoall<int>(numKeysSend, numKeysRecv, 1, comm);

    // Pre-processing for sending

    int *sendOffsets = new int[npes];
    sendOffsets[0] = 0;
    int *recvOffsets = new int[npes];
    recvOffsets[0] = 0;

    // compute offsets ...

    for (int i = 1; i < npes; i++) {
      sendOffsets[i] = sendOffsets[i-1] + 
        numKeysSend[i-1];
      recvOffsets[i] = recvOffsets[i-1] +
        numKeysRecv[i-1];
    }

    std::vector<ot::TreeNode> sendK(sendOffsets[npes-1] + numKeysSend[npes-1]);

    for(unsigned int i = 0; i < npes; i++) {
#ifdef __DEBUG_OCT__
      assert(sendKtmp[i].size() == numKeysSend[i]);
#endif
      for(unsigned int j = 0; j < numKeysSend[i]; j++) {
        sendK[sendOffsets[i] + j] = allKeys[sendKtmp[i][j]];
      }//end for j
    }//end for i

    if(partKeys) {
      delete [] partKeys;
      partKeys = NULL;
    }

    allKeys.clear();
    sendKtmp.clear();

    std::vector<ot::TreeNode> recvK(recvOffsets[npes-1] + numKeysRecv[npes-1]);

    ot::TreeNode* sendKptr = NULL;
    ot::TreeNode* recvKptr = NULL;
    if(!sendK.empty()) {
      sendKptr = &(*(sendK.begin()));
    }
    if(!recvK.empty()) {
      recvKptr = &(*(recvK.begin()));
    }
    par::Mpi_Alltoallv_sparse<ot::TreeNode>( sendKptr, numKeysSend, sendOffsets,        
        recvKptr, numKeysRecv, recvOffsets, comm);


    sendK.clear();

    delete [] sendOffsets;
    sendOffsets = NULL;

    delete [] recvOffsets;
    recvOffsets = NULL;

    delete [] numKeysSend;
    numKeysSend = NULL;

    delete [] numKeysRecv;
    numKeysRecv = NULL;

    seq::makeVectorUnique<TreeNode>(recvK,false);	

    keyLen = recvK.size();

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 7."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    //Local searches and creating seedList      

    std::vector<std::vector<TreeNode> > seedList(nodes.size());
    unsigned int lastIdx = 0;
    for(unsigned int i=0; i < keyLen; i++) {
      unsigned int idx;
      //assumes nodes is sorted and unique.          
      bool flag = seq::maxLowerBound<TreeNode>(nodes, recvK[i], idx, &lastIdx,NULL);
      if (flag) {
        lastIdx = idx;
#ifdef __DEBUG_OCT__
        assert(areComparable(nodes[idx], recvK[i]));
#endif
        if (nodes[idx].isAncestor(recvK[i])) {
          if (lev > (nodes[idx].getLevel() + 1)) {
            seedList[idx].push_back(recvK[i].getAncestor(lev-1));                
          }
        }//end if correct result
      }//end if flag        
    }//end for i

    recvK.clear();

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 8."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    //Include the new octants in the octree
    //while preserving linearity and sorted order

    //Seedlist may have duplicates. Seedlist will be sorted.

    std::vector<TreeNode> tmpList;
    std::vector<std::vector<TreeNode> >allInternalLeaves(nodes.size());
    unsigned int tmpSz = 0;
    for (unsigned int i=0;i<nodes.size();i++) {	
      seq::makeVectorUnique<TreeNode>(seedList[i],true);	
      if (!seedList[i].empty()) {
        nodes[i].completeSubtree(seedList[i], allInternalLeaves[i]);
      }
      tmpSz += allInternalLeaves[i].size();
      if (allInternalLeaves[i].empty()) {
        tmpSz++;
      }
      seedList[i].clear();
    }//end for i

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 9."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    seedList.clear();
    tmpList.resize(tmpSz);
    tmpSz=0;
    for (unsigned int i=0;i<nodes.size();i++) {
      if (!allInternalLeaves[i].empty()) {
        for (int k=0;k<allInternalLeaves[i].size();k++) {
          tmpList[tmpSz++] = allInternalLeaves[i][k];
        }//end for k
        allInternalLeaves[i].clear();
      } else {
        tmpList[tmpSz++] = nodes[i];
      }
    }//end for i
    allInternalLeaves.clear();
    nodes = tmpList;
    tmpList.clear();

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 10."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    if(rePart) {
      par::partitionW<ot::TreeNode>(nodes, NULL,comm);
    }

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 11."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

  }//end for lev

  PROF_PAR_RIPPLE_TYPE3_END 
}//end function


int parallelRippleType2(std::vector<TreeNode> & nodes,
    bool incCorners, bool checkBailOut, bool rePart,
    unsigned int dim, unsigned int maxDepth, MPI_Comm comm)
{
  PROF_PAR_RIPPLE_TYPE2_BEGIN 

    TreeNode root(dim,maxDepth);

  unsigned int maxLev = 1;
  for(unsigned int i=0; i<nodes.size(); i++) {
    if (nodes[i].getLevel() > maxLev) {
      maxLev = nodes[i].getLevel();
    }
  }//end for 

  int rank,npes;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&npes);

#ifdef __DEBUG_OCT__
  MPI_Barrier(comm);	
  if(!rank) {
    std::cout<<"Computed local maxLev."<<std::endl;
  }
  MPI_Barrier(comm);	
#endif

  unsigned int globalMaxLev = maxLev;
  par::Mpi_Allreduce<unsigned int>(&maxLev, &globalMaxLev, 1, MPI_MAX, comm);

#ifdef __DEBUG_OCT__
  MPI_Barrier(comm);	
  if(!rank) {
    std::cout<<"Global maxLev: "<<globalMaxLev<<std::endl;
  }
  MPI_Barrier(comm);	
#endif

  for (unsigned int lev = globalMaxLev; lev > 2; lev--)
  {
#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Lev: "<<lev<<std::endl;
    }
    MPI_Barrier(comm);	
#endif
    if(checkBailOut) {
      unsigned int minLev = maxDepth;
      for (unsigned int i=0; i<nodes.size(); i++) {
        if (nodes[i].getLevel() < minLev) {
          minLev = nodes[i].getLevel();
        }
      }//end for i

      unsigned int globalMinLev = minLev;
      par::Mpi_Allreduce<unsigned int>(&minLev, &globalMinLev, 1, MPI_MIN, comm);

      if (globalMinLev >= (lev-1)) {
        //Difference between min and 
        //max levels is less than 2 .
        break;
      }
    }//end if check to bail out

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 1."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    std::vector<TreeNode> wList;
    unsigned int wLen = 0;
    wList.resize(nodes.size());
    //Pick leaves at a particular level
    for (int i=0; i<nodes.size(); i++) {
      if (nodes[i].getLevel() == lev) {
        wList[wLen] = nodes[i];
        wLen++;
      }
    }//end for i
    wList.resize(wLen);

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 2."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    std::vector<std::vector<TreeNode> >tList(wLen);
    for (unsigned int i=0; i < wLen; i++) {
      tList[i] = wList[i].getSearchKeys(incCorners);        
    }//end for i

    wList.clear();

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 3."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    std::vector<TreeNode> allKeys;
    for (unsigned int i=0; i < wLen; i++) {
      for(unsigned int j=0; j < tList[i].size(); j++) {
        if(tList[i][j] > root) {
          allKeys.push_back(tList[i][j]);
        }
      }
    }    
    tList.clear();

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 4."<<std::endl;
    }
    std::cout<<rank<<": allKeys.size(): "
      <<allKeys.size()<<std::endl;
    MPI_Barrier(comm);	
#endif

    //
    //Parallel searches and subsequent processing      
    //1. Compute the ranges controlled by each processor
    //2. Compare the keys with the ranges and 
    //send the keys to the appropriate processors.
    //3. Perform local searches using the keys
    //recieved. 
    //4. Compute the balancing descendants for the results
    //and create the 'seedList' vector. Note,
    //we process one level at a time. So all the keys
    //were generated by octants at the same level. The 
    //corresponding balancing descendants are simply the
    //ancestors of the keys at one level lower than this
    //level.
    //

    unsigned int keyLen = allKeys.size();


    //First Get the mins from each processor.

    // allocate memory for the mins array
    std::vector<ot::TreeNode> mins (npes); 

    ot::TreeNode sendMin;
    if(!nodes.empty()) {
      sendMin = nodes[0]; //local min
    }else {
      sendMin = root;
    }

    par::Mpi_Allgather<ot::TreeNode>(&sendMin, &(*mins.begin()), 1, comm);    

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 5."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    //Now determine the processors which own these keys.

    unsigned int *partKeys = NULL;

    if(keyLen) {
      partKeys = new unsigned int[keyLen];    
    }

    for (unsigned int i=0; i<keyLen; i++) {
      unsigned int idx;
      //maxLB returns the last index in a 
      //sorted array such that a[ind] <= key 
      //and  a[index +1] > key
      bool found = seq::maxLowerBound<TreeNode >(mins,
          allKeys[i], idx, NULL, NULL);
      if (!found ) {
        //Can happen on incomplete domains
        partKeys[i] = rank;
      } else {
        partKeys[i] = idx;
      }
    }//end for i

    mins.clear();

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 6."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    int *numKeysSend = new int[npes];
    int *numKeysRecv = new int[npes];    

    for (int i=0; i<npes; i++) {
      numKeysSend[i] = 0;
    }

    // calculate the number of keys to send ...
    //and create the send buffer.
    std::vector<std::vector<unsigned int> >sendKtmp(npes);

    for (unsigned int i=0; i<keyLen; i++) {
      numKeysSend[partKeys[i]]++;
      sendKtmp[partKeys[i]].push_back(i);      
    }

    // Now do an All2All to get numKeysRecv

    par::Mpi_Alltoall<int>(numKeysSend, numKeysRecv, 1, comm);

    // Pre-processing for sending

    int *sendOffsets = new int[npes];
    sendOffsets[0] = 0;
    int *recvOffsets = new int[npes];
    recvOffsets[0] = 0;

    // compute offsets ...

    for (int i = 1; i < npes; i++) {
      sendOffsets[i] = sendOffsets[i-1] + 
        numKeysSend[i-1];
      recvOffsets[i] = recvOffsets[i-1] +
        numKeysRecv[i-1];
    }

    std::vector<ot::TreeNode> sendK(sendOffsets[npes-1] + numKeysSend[npes-1]);

    for(unsigned int i = 0; i < npes; i++) {
#ifdef __DEBUG_OCT__
      assert(sendKtmp[i].size() == numKeysSend[i]);
#endif
      for(unsigned int j = 0; j < numKeysSend[i]; j++) {
        sendK[sendOffsets[i] + j] = allKeys[sendKtmp[i][j]];
      }//end for j
    }//end for i

    if(partKeys) {
      delete [] partKeys;
      partKeys = NULL;
    }

    allKeys.clear();
    sendKtmp.clear();

    std::vector<ot::TreeNode> recvK(recvOffsets[npes-1] + numKeysRecv[npes-1]);

    ot::TreeNode* sendKptr = NULL;
    ot::TreeNode* recvKptr = NULL;
    if(!sendK.empty()) {
      sendKptr = &(*(sendK.begin()));
    }
    if(!recvK.empty()) {
      recvKptr = &(*(recvK.begin()));
    }
    par::Mpi_Alltoallv_sparse<ot::TreeNode>(sendKptr, numKeysSend, sendOffsets,        
        recvKptr, numKeysRecv,recvOffsets, comm);

    sendK.clear();

    delete [] sendOffsets;
    sendOffsets = NULL;

    delete [] recvOffsets;
    recvOffsets = NULL;

    delete [] numKeysSend;
    numKeysSend = NULL;

    delete [] numKeysRecv;
    numKeysRecv = NULL;

    keyLen = recvK.size();

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 7."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    //Local searches and creating seedList      

    std::vector<std::vector<TreeNode> > seedList(nodes.size());
    for(unsigned int i=0; i < keyLen; i++) {
      unsigned int idx;
      //assumes nodes is sorted and unique.          
      bool flag = seq::maxLowerBound<TreeNode>(nodes, recvK[i], idx,NULL,NULL);
      if (flag) {
#ifdef __DEBUG_OCT__
        assert(areComparable(nodes[idx], recvK[i]));
#endif
        if (nodes[idx].isAncestor(recvK[i])) {
          if (lev > (nodes[idx].getLevel() + 1)) {
            seedList[idx].push_back(recvK[i].getAncestor(lev-1));                
          }
        }//end if correct result
      }//end if flag        
    }//end for i

    recvK.clear();

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 8."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    //Include the new octants in the octree
    //while preserving linearity and sorted order

    //Seedlist may have duplicates. Seedlist 
    // may not be sorted.

    std::vector<TreeNode> tmpList;
    std::vector<std::vector<TreeNode> >allInternalLeaves(nodes.size());
    unsigned int tmpSz = 0;
    for (unsigned int i=0;i<nodes.size();i++) {	
      seq::makeVectorUnique<TreeNode>(seedList[i],false);	
      if (!seedList[i].empty()) {
        nodes[i].completeSubtree(seedList[i], allInternalLeaves[i]);
      }
      tmpSz += allInternalLeaves[i].size();
      if (allInternalLeaves[i].empty()) {
        tmpSz++;
      }
      seedList[i].clear();
    }//end for i

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 9."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    seedList.clear();
    tmpList.resize(tmpSz);
    tmpSz=0;
    for (unsigned int i=0;i<nodes.size();i++) {
      if (!allInternalLeaves[i].empty()) {
        for (int k=0;k<allInternalLeaves[i].size();k++) {
          tmpList[tmpSz++] = allInternalLeaves[i][k];
        }//end for k
        allInternalLeaves[i].clear();
      } else {
        tmpList[tmpSz++] = nodes[i];
      }
    }//end for i
    allInternalLeaves.clear();
    nodes = tmpList;
    tmpList.clear();

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 10."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    if(rePart) {
      par::partitionW<ot::TreeNode>(nodes, NULL,comm);
    }

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 11."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

  }//end for lev

  PROF_PAR_RIPPLE_TYPE2_END 
}//end function

int parallelRippleType1(std::vector<TreeNode> & nodes,
    bool incCorners, bool checkBailOut, bool rePart,
    unsigned int dim, unsigned int maxDepth, MPI_Comm comm)
{
  PROF_PAR_RIPPLE_TYPE1_BEGIN 

    TreeNode root(dim,maxDepth);

  unsigned int maxLev = 1;
  for(unsigned int i=0; i<nodes.size(); i++) {
    if (nodes[i].getLevel() > maxLev) {
      maxLev = nodes[i].getLevel();
    }
  }//end for 

  int rank,npes;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&npes);

#ifdef __DEBUG_OCT__
  MPI_Barrier(comm);	
  if(!rank) {
    std::cout<<"Computed local maxLev."<<std::endl;
  }
  MPI_Barrier(comm);	
#endif

  unsigned int globalMaxLev = maxLev;
  par::Mpi_Allreduce<unsigned int>(&maxLev,&globalMaxLev, 1, MPI_MAX, comm);

#ifdef __DEBUG_OCT__
  MPI_Barrier(comm);	
  if(!rank) {
    std::cout<<"Global maxLev: "<<globalMaxLev<<std::endl;
  }
  MPI_Barrier(comm);	
#endif

  for (unsigned int lev = globalMaxLev; lev > 2; lev--)
  {
#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Lev: "<<lev<<std::endl;
    }
    MPI_Barrier(comm);	
#endif
    if(checkBailOut) {
      unsigned int minLev = maxDepth;
      for (unsigned int i=0; i<nodes.size(); i++) {
        if (nodes[i].getLevel() < minLev) {
          minLev = nodes[i].getLevel();
        }
      }//end for i

      unsigned int globalMinLev = minLev;
      par::Mpi_Allreduce<unsigned int>(&minLev, &globalMinLev, 1, MPI_MIN, comm);

      if (globalMinLev >= (lev-1)) {
        //Difference between min and 
        //max levels is less than 2 .
        break;
      }
    }//end if check to bail out

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 1."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    std::vector<TreeNode> wList;
    unsigned int wLen = 0;
    wList.resize(nodes.size());
    //Pick leaves at a particular level
    for (int i=0; i<nodes.size(); i++) {
      if (nodes[i].getLevel() == lev) {
        wList[wLen] = nodes[i];
        wLen++;
      }
    }//end for i
    wList.resize(wLen);

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 2."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    std::vector<std::vector<TreeNode> >tList(wLen);
    for (unsigned int i=0; i < wLen; i++) {
      tList[i] = wList[i].getSearchKeys(incCorners);        
    }//end for i

    wList.clear();

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 3."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    std::vector<TreeNode> allKeys;
    for (unsigned int i=0; i < wLen; i++) {
      for(unsigned int j=0; j < tList[i].size(); j++) {
        if(tList[i][j] > root) {
          allKeys.push_back(tList[i][j]);
        }
      }
    }    
    tList.clear();

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 4."<<std::endl;
    }
    std::cout<<rank<<": allKeys.size(): "
      <<allKeys.size()<<std::endl;
    MPI_Barrier(comm);	
#endif

    //
    //Parallel searches and subsequent processing
    //1. Make the list of keys unique globally.
    //2. Compute the ranges controlled by each processor
    //3. Compare the keys with the ranges and 
    //send the keys to the appropriate processors.
    //4. Perform local searches using the keys
    //recieved. Use the fact that the keys are
    // sorted to do this efficiently.
    //5. Compute the balancing descendants for the results
    //and create the 'seedList' vector. Note,
    //we process one level at a time. So all the keys
    //were generated by octants at the same level. The 
    //corresponding balancing descendants are simply the
    //ancestors of the keys at one level lower than this
    //level.
    //

    par::removeDuplicates<ot::TreeNode>(allKeys,
        false,comm);
    unsigned int keyLen = allKeys.size();

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 5."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    //First Get the mins from each processor.

    // allocate memory for the mins array
    std::vector<ot::TreeNode> mins (npes); 

    ot::TreeNode sendMin;
    if(!nodes.empty()) {
      sendMin = nodes[0]; //local min
    }else {
      sendMin = root;
    }

    par::Mpi_Allgather<ot::TreeNode>(&sendMin, &(*mins.begin()), 1, comm);    

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 6."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    //Now determine the processors which own these keys.

    unsigned int *partKeys = NULL;

    if(keyLen) {
      partKeys = new unsigned int[keyLen];    
    }

    for (unsigned int i=0; i<keyLen; i++) {
      unsigned int idx;
      //maxLB returns the last index in a 
      //sorted array such that a[ind] <= key 
      //and  a[index +1] > key
      bool found = seq::maxLowerBound<TreeNode >(mins,
          allKeys[i], idx, NULL, NULL);
      if (!found ) {
        //Can happen on incomplete domains
        partKeys[i] = rank;
      } else {
        partKeys[i] = idx;
      }
    }//end for i

    mins.clear();

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 7."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    int *numKeysSend = new int[npes];
    int *numKeysRecv = new int[npes];    

    for (int i=0; i<npes; i++) {
      numKeysSend[i] = 0;
    }

    // calculate the number of keys to send ...
    for (unsigned int i=0; i<keyLen; i++) {
      numKeysSend[partKeys[i]]++;      
    }

    // Now do an All2All to get numKeysRecv

    par::Mpi_Alltoall<int>(numKeysSend, numKeysRecv, 1, comm);

    // Pre-processing for sending

    int *sendOffsets = new int[npes];
    sendOffsets[0] = 0;
    int *recvOffsets = new int[npes];
    recvOffsets[0] = 0;


    // compute offsets ...

    for (int i = 1; i < npes; i++) {
      sendOffsets[i] = sendOffsets[i-1] + 
        numKeysSend[i-1];
      recvOffsets[i] = recvOffsets[i-1] +
        numKeysRecv[i-1];
    }

    //Since, allKeys is sorted globally in this case
    // allKeys must be the same as sendK. So,
    //there is no need to create 
    //a separate send buffer.

    std::vector<ot::TreeNode> recvK(recvOffsets[npes-1] 
        + numKeysRecv[npes-1]);

    if(partKeys) {
      delete [] partKeys;
      partKeys = NULL;
    }

    ot::TreeNode* allKeysPtr = NULL;
    ot::TreeNode* recvKptr = NULL;
    if(!allKeys.empty()) {
      allKeysPtr = &(*(allKeys.begin()));
    }
    if(!recvK.empty()) {
      recvKptr = &(*(recvK.begin()));
    }
    par::Mpi_Alltoallv_sparse<ot::TreeNode>(allKeysPtr, numKeysSend, sendOffsets, 
        recvKptr, numKeysRecv, recvOffsets, comm);

    allKeys.clear();

    delete [] sendOffsets;
    sendOffsets = NULL;

    delete [] recvOffsets;
    recvOffsets = NULL;

    delete [] numKeysSend;
    numKeysSend = NULL;

    delete [] numKeysRecv;
    numKeysRecv = NULL;

    keyLen = recvK.size();

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 8."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    //Local searches and creating seedList
    //recvK will be sorted. 
    //since, recvK is sorted, we can pass the previous
    // result as a lower bound for subsequent 
    //searches instead of NULL for everything. This
    //should reduce the constant a bit.

    std::vector<std::vector<TreeNode> > seedList(nodes.size());
    unsigned int lastIdx = 0;	
    for(unsigned int i=0; i < keyLen; i++) {
      unsigned int idx;
      //assumes nodes is sorted and unique.          
      bool flag = seq::maxLowerBound<TreeNode>(nodes,
          recvK[i], idx,&lastIdx,NULL);
      if (flag) {
        lastIdx = idx;	
#ifdef __DEBUG_OCT__
        assert(areComparable(nodes[idx], recvK[i]));
#endif
        if (nodes[idx].isAncestor(recvK[i])) {
          if (lev > (nodes[idx].getLevel() + 1)) {
            seedList[idx].push_back(
                recvK[i].getAncestor(lev-1));                
          }
        }//end if correct result
      }//end if flag        
    }//end for i

    recvK.clear();

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 9."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    //Include the new octants in the octree
    //while preserving linearity and sorted order

    //Seedlist may have duplicates. Seedlist 
    // will be sorted.

    std::vector<TreeNode> tmpList;
    std::vector<std::vector<TreeNode> >allInternalLeaves(nodes.size());
    unsigned int tmpSz = 0;
    for (unsigned int i=0;i<nodes.size();i++) {	
      seq::makeVectorUnique<TreeNode>(seedList[i],true);	
      if (!seedList[i].empty()) {
        nodes[i].completeSubtree(seedList[i], allInternalLeaves[i]);
      }
      tmpSz += allInternalLeaves[i].size();
      if (allInternalLeaves[i].empty()) {
        tmpSz++;
      }
      seedList[i].clear();
    }//end for i

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 10."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    seedList.clear();
    tmpList.resize(tmpSz);
    tmpSz=0;
    for (unsigned int i=0;i<nodes.size();i++) {
      if (!allInternalLeaves[i].empty()) {
        for (int k=0;k<allInternalLeaves[i].size();k++) {
          tmpList[tmpSz++] = allInternalLeaves[i][k];
        }//end for k
        allInternalLeaves[i].clear();
      } else {
        tmpList[tmpSz++] = nodes[i];
      }
    }//end for i
    allInternalLeaves.clear();
    nodes = tmpList;
    tmpList.clear();

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 11."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

    if(rePart) {
      par::partitionW<ot::TreeNode>(nodes, NULL,comm);
    }

#ifdef __DEBUG_OCT__
    MPI_Barrier(comm);	
    if(!rank) {
      std::cout<<"Passed Step 12."<<std::endl;
    }
    MPI_Barrier(comm);	
#endif

  }//end for lev

  PROF_PAR_RIPPLE_TYPE1_END 
}//end function

//myNhBlocks is automatically sorted and unique at the end of the loop.
int selectNeighboringBlocks(const std::vector<TreeNode>& allBlocks, 
    const std::vector<TreeNode>& blocks, const std::vector<unsigned int>& maxBlockBndVec,
    int myRank, std::vector<TreeNode>& myNhBlocks) {
  PROF_PICK_NH_BLOCKS_BEGIN

    for (int k = 0; k < allBlocks.size(); k++) {
      if (allBlocks[k].getWeight() == myRank) {
        //Only need to send to other processors.
        continue;
      }

      unsigned int othMinX = allBlocks[k].minX();
      unsigned int othMinY = allBlocks[k].minY();
      unsigned int othMinZ = allBlocks[k].minZ();
      unsigned int othMaxX = allBlocks[k].maxX();
      unsigned int othMaxY = allBlocks[k].maxY();
      unsigned int othMaxZ = allBlocks[k].maxZ();

      //Even if allBlocks[k] interescts the region of influence of any of the
      //local blocks, it is included in the neighbour list.
      for (int j = 0; j < blocks.size(); j++) {
        unsigned int myLen = maxBlockBndVec[j];
        unsigned int myMinX = blocks[j].minX();
        unsigned int myMinY = blocks[j].minY();
        unsigned int myMinZ = blocks[j].minZ();
        unsigned int myMaxX = blocks[j].maxX();
        unsigned int myMaxY = blocks[j].maxY();
        unsigned int myMaxZ = blocks[j].maxZ();
        unsigned int xlow = ( (myMinX >= myLen) ? (myMinX - myLen) : myMinX );
        unsigned int ylow = ( (myMinY >= myLen) ? (myMinY - myLen) : myMinY );
        unsigned int zlow = ( (myMinZ >= myLen) ? (myMinZ - myLen) : myMinZ );
        unsigned int xhigh = ( myMaxX + myLen );
        unsigned int yhigh = ( myMaxY + myLen );
        unsigned int zhigh = ( myMaxZ + myLen );

        if ( (othMinX < xhigh) && (othMinY < yhigh) && (othMinZ < zhigh)
            && (othMaxX > xlow) && (othMaxY > ylow) && (othMaxZ > zlow) ) {
          myNhBlocks.push_back(allBlocks[k]);
          break;
        }//end if to be sent
      }//end for j
    }//end for k   

  PROF_PICK_NH_BLOCKS_END
}//end function

int mergeComboBalAndPickBoundary(std::vector<ot::TreeNode>& out, 
    std::vector<ot::TreeNode>& allBoundaryLeaves, 
    const ot::TreeNode& firstBlock, const ot::TreeNode& lastBlock) {
  PROF_MERGE_COMBO_BAL_BEGIN

    //Merge (in-place) results from the two stages, i.e. blockBalance and
    //intra-processor rippleBalance..
    std::vector<TreeNode> tmpNodeList(out.size() + allBoundaryLeaves.size());

  unsigned int tmpLsz = 0;
  unsigned int bndCnt = 0;
  unsigned int bndSz = allBoundaryLeaves.size();

  for (unsigned int i = 0;i < out.size();i++) {
    if ( bndCnt < allBoundaryLeaves.size() ) {
      if ( out[i] == allBoundaryLeaves[bndCnt] ) {
        tmpNodeList[tmpLsz] = out[i];
        tmpLsz++;
        bndCnt++;
      } else if (out[i] < allBoundaryLeaves[bndCnt] ) {
#ifdef __DEBUG_OCT__
        assert(areComparable(out[i], allBoundaryLeaves[bndCnt]));
#endif
        if (out[i].isAncestor(allBoundaryLeaves[bndCnt]) ) {
          while ( (bndCnt < bndSz) && 
              out[i].isAncestor(allBoundaryLeaves[bndCnt]) ) {
            tmpNodeList[tmpLsz] = allBoundaryLeaves[bndCnt];
            tmpLsz++;
            bndCnt++;
#ifdef __DEBUG_OCT__
            if(bndCnt < bndSz) {
              assert(areComparable(out[i], allBoundaryLeaves[bndCnt]));
            }
#endif
          }
        } else {
          tmpNodeList[tmpLsz] = out[i];
          tmpLsz++;
        }
      } else {
        // nodes[i] > allBdy .. so insert
        tmpNodeList[tmpLsz] = allBoundaryLeaves[bndCnt];
        tmpLsz++;
        bndCnt++;
      }
    } else {
      tmpNodeList[tmpLsz] = out[i];
      tmpLsz++;
    }
  }//end for i

  tmpNodeList.resize(tmpLsz);
  out = tmpNodeList;

  tmpNodeList = allBoundaryLeaves;
  pickInterProcessorBoundaryNodes(tmpNodeList, allBoundaryLeaves,
      firstBlock, lastBlock);
  tmpNodeList.clear();

  PROF_MERGE_COMBO_BAL_END 
}//end function

int finalMergeInBal(std::vector<ot::TreeNode>& out, std::vector<ot::TreeNode>& allBoundaryLeaves) {
  PROF_FINAL_MERGE_IN_BAL_BEGIN

    std::vector<ot::TreeNode> tmpNodeList(out.size() + allBoundaryLeaves.size());

  unsigned int tmpLsz = 0;
  unsigned int bndSz = allBoundaryLeaves.size();
  unsigned int bndCnt = 0;
  for (unsigned int i = 0; i < out.size(); i++) {
    if ( bndCnt < bndSz ) {
      if ( out[i] == allBoundaryLeaves[bndCnt] ) {
        tmpNodeList[tmpLsz++] = out[i];
        bndCnt++;
      } else if (out[i] < allBoundaryLeaves[bndCnt] ) {
#ifdef __DEBUG_OCT__
        assert(areComparable(out[i], allBoundaryLeaves[bndCnt]));
#endif
        if (out[i].isAncestor(allBoundaryLeaves[bndCnt]) ) {
          while ( (bndCnt < bndSz ) && out[i].isAncestor(allBoundaryLeaves[bndCnt]) ) {
            tmpNodeList[tmpLsz++] = allBoundaryLeaves[bndCnt++];
#ifdef __DEBUG_OCT__
            if(bndCnt < bndSz) {
              assert(areComparable(out[i], allBoundaryLeaves[bndCnt]));
            }
#endif
          }
        } else {
          tmpNodeList[tmpLsz++] = out[i];
        }
      } else {
        // nodes[i] > allBdy .. so insert
        tmpNodeList[tmpLsz++] = allBoundaryLeaves[bndCnt++];
      }
    } else {
      tmpNodeList[tmpLsz++] = out[i];
    }
  }//end for i

  tmpNodeList.resize(tmpLsz);
  out = tmpNodeList;
  tmpNodeList.clear();
  allBoundaryLeaves.clear();

  PROF_FINAL_MERGE_IN_BAL_END
}//end function

int prepareBalComm1MessagesType1(const std::vector<ot::TreeNode>& allBoundaryLeaves, 
    const std::vector<ot::TreeNode>& myNhBlocks, int npes, unsigned int maxDepth, 
    std::vector<TreeNode>* sendNodes, std::vector<unsigned int>* sentToPid, int* sendCnt) {
  PROF_PREP_BAL_COMM1_MSSG_BEGIN 

    unsigned int allBndSz = allBoundaryLeaves.size();
  for (int j = 0; j < allBndSz; j++) {
    unsigned int myLen = (1u << (maxDepth - (allBoundaryLeaves[j].getLevel())));
    unsigned int myMinX = allBoundaryLeaves[j].minX();
    unsigned int myMinY = allBoundaryLeaves[j].minY();
    unsigned int myMinZ = allBoundaryLeaves[j].minZ();
    unsigned int xlow = ( (myMinX >= myLen) ? (myMinX - myLen) : myMinX );
    unsigned int ylow = ( (myMinY >= myLen) ? (myMinY - myLen) : myMinY );
    unsigned int zlow = ( (myMinZ >= myLen) ? (myMinZ - myLen) : myMinZ );
    unsigned int xhigh = ( myMinX + (2*myLen) );
    unsigned int yhigh = ( myMinY + (2*myLen) );
    unsigned int zhigh = ( myMinZ + (2*myLen) );

    unsigned int lastP = npes;
    //Each node could be sent to multiple processors. So must check with
    //different blocks. All the blocks from the same processor are
    //contiguous.
    for (int k = 0; k < myNhBlocks.size(); k++) {
      if (myNhBlocks[k].getWeight() == lastP) {
        //Already sent to this processor. It doesn't matter if there are
        //multiple blocks on the same processor in the region of influence of
        //the same node.
        continue;
      }//end if sent already.

      unsigned int othMinX = myNhBlocks[k].minX();
      unsigned int othMinY = myNhBlocks[k].minY();
      unsigned int othMinZ = myNhBlocks[k].minZ();
      unsigned int othMaxX = myNhBlocks[k].maxX();
      unsigned int othMaxY = myNhBlocks[k].maxY();
      unsigned int othMaxZ = myNhBlocks[k].maxZ();

      if ( (othMinX < xhigh) && (othMinY < yhigh) && (othMinZ < zhigh)
          && (othMaxX > xlow) && (othMaxY > ylow) && (othMaxZ > zlow) ) {
        sendNodes[myNhBlocks[k].getWeight()].push_back(allBoundaryLeaves[j]);
        sentToPid[j].push_back(myNhBlocks[k].getWeight());
        sendCnt[myNhBlocks[k].getWeight()]++;
        lastP = myNhBlocks[k].getWeight();
      }//end if to be sent
    }//end for k   
  }//end for j

  PROF_PREP_BAL_COMM1_MSSG_END 
}//end function

int prepareWlistInBal(const std::vector<ot::TreeNode>& recvK1, 
    const int* recvCnt, int npes, const ot::TreeNode& myFirstBlock,
    const ot::TreeNode& myLastBlock, std::vector<TreeNode>& wList,
    std::vector<std::vector<unsigned int> >& wListRanks) {
  PROF_PREP_BAL_WLIST_BEGIN

    std::vector<ot::TreeNode> recvKtemp = recvK1;

  int counter = 0;
  for(int i = 0; i < npes; i++) {
    for(int j = counter; j < (counter + recvCnt[i]); j++) {
      recvKtemp[j].setWeight(i);
    }//end for j
    counter += recvCnt[i];
  }//end for i

  std::vector<std::vector<ot::TreeNode> > wListTmp(recvKtemp.size());
  for (int i = 0; i < recvKtemp.size(); i++) {
    std::vector<TreeNode> myNh =  recvKtemp[i].getAllNeighbours();
    unsigned int myNhSz = myNh.size();
    for (int k = 0; k < myNhSz; k++) {
      //Those that overlap the local domain...
      //There are only 3 types of overlaps between octants A and B.
      //a) A is an ancestor of B
      //b) B is an ancestor of A
      //c) A and B are the same.
      //Any octant >= myFirstBlock and <= myLastBlock.DLD if and only if it
      //lies on my processor and it is either a decendant of one of my blocks or equal to
      //one of my blocks. 
      //Fact: If A is an ancestor of B, then A must also be the ancestor of
      //any octant > A and < B. So, if myNh[k] is not an ancestor of
      //myFirstBlock, then it can be an ancestor of some octant > myFirstBlock only
      //if myNh[k] >= myFirstBlock
#ifdef __DEBUG_OCT__
      assert(areComparable(myNh[k], myFirstBlock));
#endif
      if ( (myNh[k].isAncestor(myFirstBlock)) ||
          ( (myNh[k] >= myFirstBlock) && (myNh[k] <= myLastBlock.getDLD()) ) ) {
        myNh[k].setWeight(recvKtemp[i].getWeight());
        wListTmp[i].push_back(myNh[k]);
      }
    }//end for k
    myNh.clear();
  }//end for i

  recvKtemp.clear();

  for (int i = 0;i < wListTmp.size(); i++) {
    wList.insert(wList.end(), wListTmp[i].begin(), wListTmp[i].end());
    wListTmp[i].clear();
  }
  wListTmp.clear();

  std::sort(wList.begin(),wList.end());

  unsigned int wCtr=0;
  unsigned int wListSz = wList.size();
  if (wCtr < wListSz) {
    std::vector<unsigned int> tmpWr;
    tmpWr.push_back(wList[0].getWeight());
    wList[0].setWeight(1);
    wListRanks.push_back (tmpWr);
    tmpWr.clear();
    wCtr++;
  }

  unsigned int rowCtr = 0;
  while (wCtr < wListSz) {
    if (wList[wCtr] == wList[wCtr-1]) {
      wListRanks[rowCtr].push_back(wList[wCtr].getWeight());
      wList[wCtr].setWeight(1);
    } else {
      seq::makeVectorUnique<unsigned int>(wListRanks[rowCtr],false);
      std::vector<unsigned int> tmpWr;
      tmpWr.push_back(wList[wCtr].getWeight());
      wList[wCtr].setWeight(1);
      wListRanks.push_back(tmpWr);
      tmpWr.clear();
      rowCtr++;
    }
    wCtr++;
  }//end while

  seq::makeVectorUnique<TreeNode>(wList, true);

  PROF_PREP_BAL_WLIST_END
}//end function

/*
   Some Facts: allBoundaryLeaves is linear, sorted,unique, incomplete
   wList is not linear (overlaps are allowed), it is sorted, unique and incomplete.
   All elements in wList are in the domain controlled by this processor.
   wList and wListRanks are in sync always.
   Each element in allBoundaryLeaves can overlap with multiple elements in wList and vice-versa.
   Important Property: If A < B < C and if A and B do not overlap then A and C also do not overlap.
   This can be proved by contradiction, Suppose A and C overlap,
   then either A is an ancestor of C or C is an ancestor of A or both are equal.
   The third is automatically ruled out. The second is not possible since C > A.
   If A was an ancestor of C then A has to be an ancestor of B as well and since this 
   is not the case, A and C can't overlap.
   */
int prepareBalComm2Messages(const std::vector<ot::TreeNode>& allBoundaryLeaves,
    const std::vector<ot::TreeNode>& wList,
    const std::vector<std::vector<unsigned int> >& wListRanks,
    std::vector<TreeNode>* sendNodes, std::vector<unsigned int>* sentToPid, int* sendCnt) {
  PROF_PREP_BAL_COMM2_MSSG_BEGIN 

    /*
       Why is sendNodes[i] sorted at the end of this loop?
       1. For each processor, the blocks (W) corresponding to that
       processor are visited in a sorted order.

       2. For any given processor say W1 is visited before W2 => W1 < W2.

       3. Either W1 is an ancestor of W2 or W1 and W2 do not overlap

       4. All bnd sent to W1 are sent in a sorted order.

       5. Let B1 be the list of elements that overlap with W1 and 
       let B2 be the ones that overlap with W2.

       6. There are only two possibilites: either B2 is a subset of B1 or
       B1 and B2 have no common elements and min(B2) > max(B1).

       7. IF W2 is a decendant of W1, then clearly any ancestor of W2 or
       decendant of W2 also overlaps W1, thus B2 is a subset of B1. Since, all 
       of these were sent already they will not be sent to the same processor 
       again and thus they will not disturb the order.

       8. If W2 > W1 and W2 and W1 do not overlap, then the only octants that overlap both
       W1 and W2 are the common ancestors.

       9. If there is a common ancestor of both in allBnd, it will be sent during the comparison
       with W1 and no other elements in allBnd will overlap either W1 or W2. In this case B1=B2.

       10. If there are no common ancestors then B1 and B2 have no common elements.

       11. Since W2 >W1 any octant that overlaps W2 but not W1 has to be > all decendants of W1.
       Hence, min(B2) > max(B1).
       */

    int lastStart = 0;
  unsigned int wListSz = wList.size();
  unsigned int allBndSz = allBoundaryLeaves.size();
  for (int ii = 0; ii < wListSz; ii++) {
    for (int j = lastStart; j < allBndSz; j++) {
      //Each node could be sent to multiple processors.
#ifdef __DEBUG_OCT__
      assert(areComparable(wList[ii], allBoundaryLeaves[j]));
#endif
      if ((wList[ii] == allBoundaryLeaves[j]) || (wList[ii].isAncestor(allBoundaryLeaves[j]))
          || (allBoundaryLeaves[j].isAncestor(wList[ii]))) {
        //Overlap...
        for (int jj = 0; jj < wListRanks[ii].size(); jj++) {
          bool sentAlready = false;
          for (int kk = 0; kk < sentToPid[j].size(); kk++) {
            //loopCtr++;
            if (sentToPid[j][kk] == wListRanks[ii][jj]) {
              sentAlready = true;
              break;
            }
          }//end for kk
          if (!sentAlready) {
            sendNodes[wListRanks[ii][jj]].push_back(allBoundaryLeaves[j]);
            sendCnt[wListRanks[ii][jj]]++;
            sentToPid[j].push_back(wListRanks[ii][jj]);
          }
        }//end for jj
      } else if (wList[ii] < allBoundaryLeaves[j]) {
        //No overlap w < allBnd
        //Since, W does not intersect this it will not intersect any element that follows.
        //This is the justification for early termination.
        break;
      } else {
        //No overlap allBnd < w
        //Since, allBnd does not intersect this w, it will not intersect any
        //element in w that follows. So subsequent w need not be compared
        //against this.
        lastStart++;
      }//end if-else overlaps
    }//end for j
  }//end for ii

  PROF_PREP_BAL_COMM2_MSSG_END 
}//end function

int mergeRecvKeysInBal(const std::vector<ot::TreeNode>& recvK1, const int* recvOffsets1,
    const std::vector<ot::TreeNode>& recvK2, const int* recvOffsets2, 
    int npes, std::vector<ot::TreeNode>& recvK) {
  PROF_MERGE_RECV_KEYS_BAL_BEGIN 

    recvK.resize(recvK1.size() + recvK2.size());

  //Note, you only recieve from other processors and not from yourself.
  //Merge recvK1 and recvK2 inplace....

  //Basic idea... All elements from processor i are less than those from i+1.
  //The elements in recvK1 and recvK2 are independently sorted.
  unsigned int recvKcnt = 0;
  for (int i = 0; i < npes; i++) {
    unsigned int nextFrom1 = recvOffsets1[i];
    unsigned int nextFrom2 = recvOffsets2[i];
    unsigned int end1 = ((i <(npes - 1)) ? recvOffsets1[i+1] : recvK1.size());
    unsigned int end2 = ((i <(npes - 1)) ? recvOffsets2[i+1] : recvK2.size());
    while ( (nextFrom1 < end1) || (nextFrom2 < end2) ) {

      if (nextFrom1 < end1 && nextFrom2 < end2) {
        while (recvK1[nextFrom1] <= recvK2[nextFrom2]) {
          recvK[recvKcnt++] = recvK1[nextFrom1];
          nextFrom1++;
          if (nextFrom1 >= end1) {
            break;
          }
        }
      }

      if (nextFrom1 < end1 && nextFrom2 < end2) {
        while (recvK1[nextFrom1] > recvK2[nextFrom2]) {
          recvK[recvKcnt++] = recvK2[nextFrom2];
          nextFrom2++;
          if (nextFrom2 >= end2) {
            break;
          }
        }
      }

      if (nextFrom2 >= end2) {
        while (nextFrom1 < end1) {
          recvK[recvKcnt++] = recvK1[nextFrom1];
          nextFrom1++;
        }
      }

      if (nextFrom1 >= end1) {
        while (nextFrom2 < end2) {
          recvK[recvKcnt++] = recvK2[nextFrom2];
          nextFrom2++;
        }
      }

    }
  }//end for i

  PROF_MERGE_RECV_KEYS_BAL_END 
}//end function

int prepareBalComm1MessagesType2(const std::vector<ot::TreeNode>& allBoundaryLeaves, 
    const std::vector<ot::TreeNode>& minsAllBlocks, int rank, unsigned int dim, 
    unsigned int maxDepth, std::vector<TreeNode>* sendNodes,
    std::vector<unsigned int>* sentToPid, int* sendCnt) {
  PROF_PREP_BAL_COMM1_MSSG_BEGIN 

    //Each octant must be sent to all processors which overlap its insulation
    //layer. So we generate all the neighbours (nh) of this octant at this level and
    //find all processors that lie in between maxLowerBound(DFD(nh(i))) and
    //maxLowerBound(DLD(nh(i))) for all i

    unsigned int allBndSz = allBoundaryLeaves.size();
  ot::TreeNode rootNode(dim, maxDepth);
  for (int j = 0; j < allBndSz; j++) {
    std::vector<ot::TreeNode> myNh = allBoundaryLeaves[j].getAllNeighbours();

    seq::makeVectorUnique<ot::TreeNode>(myNh, false);
    unsigned int stIdx = ( (myNh[0] == rootNode) ? 1 : 0 );

    std::vector<unsigned int> pIds;
    for(unsigned int i = stIdx; i < myNh.size(); i++) {
      unsigned int idx1;
      unsigned int idx2;
      seq::maxLowerBound<ot::TreeNode>(minsAllBlocks, myNh[i].getDFD(), idx1, NULL, NULL);
      seq::maxLowerBound<ot::TreeNode>(minsAllBlocks, myNh[i].getDLD(), idx2, NULL, NULL);
      for(int k = idx1; k <= idx2; k++) {
        pIds.push_back(k);
      }//end for k
    }//end for i

#ifdef __DEBUG_OCT__
    //myNh is explicitly sorted. Moreover, no two elements of myNh overlap and since
    //myNh(i) < myNh(i+1) this implies that
    //myNh(i).getDFD() <= myNh(i).getDLD() < myNh(i+1).getDFD() <=
    //myNh(i+1).getDLD() 
    //If a < b then MLB(a) <= MLB(b). So uniqueness if not guaranteed.
    assert(seq::test::isSorted<unsigned int>(pIds));
#endif

    seq::makeVectorUnique<unsigned int>(pIds, true);

    for(int i = 0; i < pIds.size(); i++) {
      if(pIds[i] != rank) {
        sendNodes[pIds[i]].push_back(allBoundaryLeaves[j]);
        sentToPid[j].push_back(pIds[i]);
        sendCnt[pIds[i]]++;
      }
    }//end for i
  }//end for j

  PROF_PREP_BAL_COMM1_MSSG_END 
}//end function

}//end namespace



