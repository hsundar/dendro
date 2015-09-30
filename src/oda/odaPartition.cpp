
/**
  @file odaPartition.C
  @brief Partition related functions for octree meshing
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  @author Hari Sundar, hsundar@gmail.com
  */

#include "oda.h"
#include "parUtils.h"
#include "seqUtils.h"
#include "testUtils.h"

#ifdef __DEBUG__
#ifndef __DEBUG_DA__
#define __DEBUG_DA__
#endif
#endif

namespace ot {

  //If block Part equals Morton part, then we only correct singular blocks with
  //hanging anchors. Unlike in the general Block Part, the blocks and the
  //nodes are already aligned in Morton Part
#ifdef __BLOCK_PART_EQUALS_MORTON_PART__

  int DA_blockPartStage2(std::vector<TreeNode> &nodes, std::vector<TreeNode> &globalCoarse,
      unsigned int dim, unsigned int maxDepth, MPI_Comm commActive) {
#ifdef __PROF_WITH_BARRIER__
    MPI_Barrier(commActive);
#endif
    PROF_BLKPART2_BEGIN

      assert(!globalCoarse.empty());
    assert(!nodes.empty());
    assert(globalCoarse[0] == nodes[0]);
    assert(globalCoarse[globalCoarse.size() - 1] == nodes[nodes.size() - 1]);

    int npesActive, rankActive;

    MPI_Comm_rank(commActive, &rankActive);
    MPI_Comm_size(commActive, &npesActive);

    //Reset weights
    std::vector<bool> isSingular(globalCoarse.size());
    unsigned int nodeCnt = 0;
    for (int i = 0; i < globalCoarse.size(); i++) {
      globalCoarse[i].setWeight(1);
      isSingular[i] = false;
      while(nodes[nodeCnt] < globalCoarse[i]) {
        nodeCnt++;
      }
      if( (nodes[nodeCnt] == globalCoarse[i]) && 
          (!(nodes[nodeCnt].getFlag() & ot::TreeNode::NODE))  ) {
        isSingular[i] = true;
      }
    }//end for i 

    if(npesActive > 1) {
      //For DA only.....
      //Pick singular blocks on this processor...
      std::vector<ot::TreeNode> singularBlocks;
      for(unsigned int i = 0; i < globalCoarse.size(); i++) {
        if(isSingular[i]) {
          singularBlocks.push_back(globalCoarse[i]);
        }
      }//end for i

      //Gather all singular blocks on all processors.
      std::vector<int> numSingular(npesActive);
      std::vector<int> singularDisps(npesActive);

      int singularSz = singularBlocks.size();

      par::Mpi_Allgather<int>(&singularSz, &(*numSingular.begin()), 1, commActive);

      unsigned int totSingular = 0;
      for(int i = 0; i < npesActive; i++) {
        totSingular += numSingular[i];
      }

      std::vector<TreeNode> allSingular(totSingular);

      singularDisps[0] = 0;
      for (unsigned int i=1; i < npesActive; i++) {
        singularDisps[i] = singularDisps[i-1] + numSingular[i-1];
      }

      ot::TreeNode* singularBlocksPtr = NULL;
      ot::TreeNode* allSingularPtr = NULL;
      if(!singularBlocks.empty()) {
        singularBlocksPtr = &(*(singularBlocks.begin()));
      }
      if(!allSingular.empty()) {
        allSingularPtr = &(*(allSingular.begin()));
      }
      par::Mpi_Allgatherv<ot::TreeNode>(singularBlocksPtr, singularSz,
          allSingularPtr, &(*numSingular.begin()), &(*singularDisps.begin()), commActive);

      singularBlocks.clear();
      numSingular.clear();
      singularDisps.clear();

#ifdef __DEBUG_DA__
      MPI_Barrier(commActive);
      assert(seq::test::isUniqueAndSorted(allSingular));
      assert(par::test::isUniqueAndSorted(globalCoarse, commActive));
      MPI_Barrier(commActive);
#endif

      //Loop through globalCoarse and set wts of all elements in between some
      //singular Block's parent and the singular Block to 0. So that the global
      //scan of all these elements in partW is the same and hence they will be
      //sent to the same processor...
      unsigned int lastIdxFound = (globalCoarse.size() -1);
      for(int singCnt = (allSingular.size()-1); singCnt >= 0; singCnt--) {
        unsigned int idxMLB;          
        bool foundMLB = seq::maxLowerBound<ot::TreeNode>(globalCoarse, 
            allSingular[singCnt], idxMLB, NULL, &lastIdxFound);
        if(foundMLB) {
          ot::TreeNode requiredOct = allSingular[singCnt].getParent().getDFD().
            getAncestor(allSingular[singCnt].getLevel());
          while(globalCoarse[idxMLB] > requiredOct) {
            globalCoarse[idxMLB].setWeight(0);
            if(idxMLB > 0) {
              idxMLB--;
            }else {
              break;
            }              
          }
          lastIdxFound = idxMLB;
          while( (singCnt >= 0) && (allSingular[singCnt] > requiredOct) ) {
            singCnt--;
          }
          singCnt++;
        }else {
          break;
        }//end if found
      }//end for i

      allSingular.clear();
    }//end if npes > 1

    isSingular.clear();

    par::partitionW<ot::TreeNode>(globalCoarse, getNodeWeight, commActive);

    //Reset weights
    for (unsigned int i = 0; i < globalCoarse.size(); i++) {
      globalCoarse[i].setWeight(1);
    }

    PROF_BLKPART2_END
  } // end blockPart

#else

  int DA_blockPartStage2(std::vector<TreeNode> &nodes, std::vector<TreeNode> &globalCoarse,
      unsigned int dim, unsigned int maxDepth, MPI_Comm commActive) {
#ifdef __PROF_WITH_BARRIER__
    MPI_Barrier(commActive);
#endif
    PROF_BLKPART2_BEGIN

      int npesActive, rankActive;

    MPI_Comm_rank(commActive, &rankActive);
    MPI_Comm_size(commActive, &npesActive);

    int *sendCnt = new int[npesActive];
    int *recvCnt = new int[npesActive];
    int *sendOffsets = new int[npesActive];
    int *recvOffsets = new int[npesActive];

    //1. Now compute the wts of these cells ...
    //    a. Get the min and max nodes at each processor.
    std::vector<TreeNode> _mins_maxs(2*npesActive);

    // communicate ...
    TreeNode sendMinMax[2];
    TreeNode rootNode (dim,maxDepth);

    if (!nodes.empty()) {
      sendMinMax[0] =  nodes[0];
      sendMinMax[1] =  nodes[nodes.size()-1];
    } else {
      sendMinMax[0] = rootNode;
      sendMinMax[1] = rootNode;
    }

    par::Mpi_Allgather<ot::TreeNode>(sendMinMax, &(*_mins_maxs.begin()), 2, commActive);

    std::vector<std::vector<TreeNode> > sendNodes(npesActive);
    std::vector<std::vector<unsigned int> > keymap(npesActive);

    for (int i = 0; i < npesActive; i++) {
      sendCnt[i] = 0;
      sendNodes[i].clear();
      keymap[i].clear();	
    }

    //    b. Now compute which cells go to which cells ...
    //       logic is that if the coarse cell is between the min and max at a
    //       processor or if it is an ancestor of min, then it is sent to that
    //       processor.
    //Naive Logic:
    for (unsigned int i = 0; i < globalCoarse.size(); i++) {
      for (int p = 0; p < npesActive; p++) {
        if ( (globalCoarse[i].isAncestor(_mins_maxs[2*p])) ||
            ( (globalCoarse[i] >= _mins_maxs[2*p]) &&
              (globalCoarse[i] <=_mins_maxs[(2*p)+1]) ) ) {
          sendNodes[p].push_back(globalCoarse[i]);
          // save keymap so that we can assign weights back to globalCoarse.
          keymap[p].push_back(i);    
          sendCnt[p]++;
        }//end if
      }//end for
    }//end for

    _mins_maxs.clear();

    //2. Send nodes to all cells to compute the wts ... locally ...

    //    a. Communicate how many you'll be sending and how many will be
    //       received.

    // Now do an All2All to get numKeysRecv
    par::Mpi_Alltoall<int>( sendCnt, recvCnt, 1, commActive);

    //    b. Concatenate all nodes into one single Carray ...
    unsigned int totalSend = 0;
    unsigned int totalRecv = 0;
    for (unsigned int i = 0; i < npesActive; i++) {
      totalSend+= sendCnt[i];
      totalRecv+= recvCnt[i];
    }

    // create the send and recv buffers ...
    std::vector<ot::TreeNode> sendK (totalSend);
    std::vector<ot::TreeNode> recvK (totalRecv);

    // Now create sendK
    sendOffsets[0] = 0;
    recvOffsets[0] = 0;

    // compute offsets ...
    for (int i = 1; i < npesActive; i++) {
      sendOffsets[i] = sendOffsets[i-1] + sendCnt[i-1];
      recvOffsets[i] = recvOffsets[i-1] + recvCnt[i-1];
    }

    for (int i = 0; i < npesActive; i++) {
#ifdef __DEBUG_DA__
      assert( sendCnt[i]  == sendNodes[i].size() );
#endif
      for (unsigned int j=0; j<sendCnt[i]; j++) {
#ifdef __DEBUG_DA__
        assert( (sendOffsets[i] + j) < totalSend);
#endif
        sendK[sendOffsets[i] + j] = sendNodes[i][j];
      }//end for j
    }//end for i

    //3. send and receive all keys ...

    ot::TreeNode* sendKptr = NULL;
    ot::TreeNode* recvKptr = NULL;
    if(!sendK.empty()) {
      sendKptr = &(*(sendK.begin()));
    }
    if(!recvK.empty()) {
      recvKptr = &(*(recvK.begin()));
    }
    par::Mpi_Alltoallv_sparse<ot::TreeNode>( sendKptr, sendCnt, sendOffsets,
        recvKptr, recvCnt, recvOffsets, commActive);

    sendK.clear();

    //4. Now compute the wts of the locally received nodes ...
    //    a. loop through nodes and update the wts of the local chunks ...
    unsigned int *wts = NULL;
    char * isAnchorHanging = NULL;
    if(totalRecv) {
      wts = new unsigned int [totalRecv];
      isAnchorHanging = new char [totalRecv];
    }

    for (unsigned int i = 0; i < totalRecv; i++) {
      wts[i] = 0;
      isAnchorHanging[i] = 0;
    }

    //decendants and chunks are both sorted at this point.
    unsigned int nextPt = 0;
    unsigned int nextNode = 0;
    //Every element in nodes is inside some element in recvK.
    while (nextPt < nodes.size()) {
      //The first pt. lies in some block.
#ifdef __DEBUG_DA__
      assert(nextNode < recvK.size());
#endif
      if ((recvK[nextNode].isAncestor(nodes[nextPt])) ||
          (recvK[nextNode] == nodes[nextPt])) {
        wts[nextNode]++;
        if( (recvK[nextNode].getAnchor() == nodes[nextPt].getAnchor()) &&
            (!(nodes[nextPt].getFlag() & ot::TreeNode::NODE)) ) {
          isAnchorHanging[nextNode] = 1;
#ifdef __DEBUG_DA__
          //Only singular blocks can have hanging anchors
          assert(recvK[nextNode] == nodes[nextPt]);
#endif
        }
        nextPt++;
      } else {
        nextNode++;
        if (nextNode >= totalRecv) {
          //If this fails then either recvK and nodes are not sorted or
          //Some pt in nodes is not in any of recvK
          assert(false);
        }
      }//end if-else
    }//end while

    recvK.clear();

    //5. Now communicate the wts back to the procs ...
    unsigned int *recvWts = NULL;
    char *recvChars = NULL;
    if(totalSend) {
      recvWts = new unsigned int[totalSend];
      recvChars = new char[totalSend];
    }

    par::Mpi_Alltoallv_sparse<unsigned int>( wts, recvCnt, recvOffsets, 
        recvWts, sendCnt, sendOffsets,  commActive);

    par::Mpi_Alltoallv_sparse<char>( isAnchorHanging, recvCnt, recvOffsets, 
        recvChars, sendCnt, sendOffsets, commActive);

    //6. Now map them back to the blocks ...
    std::vector<bool> isSingular(globalCoarse.size());
    for (int i = 0; i < globalCoarse.size(); i++) {
      globalCoarse[i].setWeight(0);
      isSingular[i] = false;
    }

    for (int i = 0; i < npesActive; i++) {
      for (int j = 0; j < sendCnt[i]; j++) {
#ifdef __DEBUG_DA__
        assert(j < keymap[i].size());
        assert(keymap[i][j] < globalCoarse.size());
        assert( (sendOffsets[i] + j) < totalSend );
#endif
        globalCoarse[keymap[i][j]].addWeight(recvWts[sendOffsets[i] + j]);
        isSingular[keymap[i][j]] = ( isSingular[keymap[i][j]] ||
            recvChars[sendOffsets[i] + j]  );
      }//end for j
    }//end for i

    for (unsigned int i = 0; i < npesActive; i++) {
      keymap[i].clear();
      sendNodes[i].clear();
    }//end for i

    sendNodes.clear();
    keymap.clear();

    if(recvWts) {
      delete [] recvWts;
      recvWts = NULL;
    }

    if(recvChars) {
      delete [] recvChars;
      recvChars = NULL;
    }

    if(wts) {
      delete [] wts;
      wts = NULL;
    }

    if(isAnchorHanging) {
      delete [] isAnchorHanging;
      isAnchorHanging = NULL;
    }

    if(npesActive > 1) {
      //For DA only.....
      //Pick singular blocks on this processor...
      std::vector<ot::TreeNode> singularBlocks;
      for(unsigned int i=0;i<globalCoarse.size(); i++) {
        if(isSingular[i]) {
          singularBlocks.push_back(globalCoarse[i]);
        }
      }//end for i

      //Gather all singular blocks on all processors.
      std::vector<int> numSingular(npesActive);
      std::vector<int> singularDisps(npesActive);

      int singularSz = singularBlocks.size();

      par::Mpi_Allgather<int>(&singularSz, &(*numSingular.begin()), 1, commActive);

      unsigned int totSingular = 0;
      for(int i = 0; i < npesActive; i++) {
        totSingular += numSingular[i];
      }

      std::vector<TreeNode> allSingular(totSingular);

      singularDisps[0] = 0;
      for (unsigned int i=1; i < npesActive; i++) {
        singularDisps[i] = singularDisps[i-1] + numSingular[i-1];
      }

      ot::TreeNode* singularBlocksPtr = NULL;
      ot::TreeNode* allSingularPtr = NULL;
      if(!singularBlocks.empty()) {
        singularBlocksPtr = &(*(singularBlocks.begin()));
      }
      if(!allSingular.empty()) {
        allSingularPtr = &(*(allSingular.begin()));
      }
      par::Mpi_Allgatherv<ot::TreeNode>(singularBlocksPtr, singularSz,
          allSingularPtr, &(*numSingular.begin()), &(*singularDisps.begin()), commActive);

      singularBlocks.clear();
      numSingular.clear();
      singularDisps.clear();

#ifdef __DEBUG_DA__
      MPI_Barrier(commActive);
      assert(seq::test::isUniqueAndSorted(allSingular));
      assert(par::test::isUniqueAndSorted(globalCoarse, commActive));
      MPI_Barrier(commActive);
#endif

      //Loop through globalCoarse and set wts of all elements in between some
      //singular Block's parent and the singular Block to 0. So that the global
      //scan of all these elements in partW is the same and hence they will be
      //sent to the same processor...
      unsigned int lastIdxFound = (globalCoarse.size() -1);
      for(int singCnt = (allSingular.size()-1); singCnt >= 0; singCnt--) {
        unsigned int idxMLB;          
        bool foundMLB = seq::maxLowerBound<ot::TreeNode>(globalCoarse, 
            allSingular[singCnt], idxMLB, NULL, &lastIdxFound);
        if(foundMLB) {
          ot::TreeNode requiredOct = allSingular[singCnt].getParent().getDFD().
            getAncestor(allSingular[singCnt].getLevel());
          while(globalCoarse[idxMLB] > requiredOct) {
            globalCoarse[idxMLB].setWeight(0);
            if(idxMLB > 0) {
              idxMLB--;
            }else {
              break;
            }              
          }
          lastIdxFound = idxMLB;
          while( (singCnt >= 0) && (allSingular[singCnt] > requiredOct) ) {
            singCnt--;
          }
          singCnt++;
        }else {
          break;
        }//end if found
      }//end for i

      allSingular.clear();
    }//end if npes > 1

    isSingular.clear();

    par::partitionW<ot::TreeNode>(globalCoarse, getNodeWeight, commActive);

    //Reset weights
    for (unsigned int i=0;i<globalCoarse.size(); i++) {
      globalCoarse[i].setWeight(1);
    }

    // clean up ...
    delete [] sendCnt;
    sendCnt = NULL;

    delete [] recvCnt;
    recvCnt = NULL;

    delete [] sendOffsets;
    sendOffsets = NULL;

    delete [] recvOffsets;
    recvOffsets = NULL;

    PROF_BLKPART2_END
  } // end blockPart

#endif

  int DA_blockPartStage3(std::vector<TreeNode> &nodes, std::vector<TreeNode>& globalCoarse,
      std::vector<ot::TreeNode>& minsAllBlocks, unsigned int dim,
      unsigned int maxDepth, MPI_Comm commActive) {
#ifdef __PROF_WITH_BARRIER__
    MPI_Barrier(commActive);
#endif
    PROF_BLKPART3_BEGIN

      int npesActive, rankActive;

    MPI_Comm_rank(commActive, &rankActive);
    MPI_Comm_size(commActive, &npesActive);

    int *sendCnt = new int[npesActive];
    int *recvCnt = new int[npesActive];
    int *sendOffsets = new int[npesActive];
    int *recvOffsets = new int[npesActive];

    TreeNode rootNode (dim,maxDepth);

    // Now communicate the nodes ...

    //7. Determine locally which keys to send to which proc ...
    //Compute Dist on globalCoarse....

    TreeNode *sendMin;
    std::vector<ot::TreeNode> vtkDist(npesActive);
    if (!globalCoarse.empty()) {
      sendMin = (TreeNode *)&(*(globalCoarse.begin()));
    } else {
      sendMin = &(rootNode);
    }

    par::Mpi_Allgather<ot::TreeNode>(sendMin, &(* vtkDist.begin()), 1, commActive);

    minsAllBlocks.clear();
    for(int j = 0; j < npesActive; j++) {
      if(vtkDist[j] != rootNode) {
        minsAllBlocks.push_back(vtkDist[j]);
      }
    }//end for j

    for (unsigned int j = 1; j < npesActive ; j++) {
      if (vtkDist[j] == rootNode) {
        vtkDist[j] = vtkDist[j-1];
      }
    }//end for j

    // correct dist ...
    if (npesActive > 1) {
      if (vtkDist[npesActive - 1] == vtkDist[npesActive - 2]) {
        vtkDist[npesActive - 1] = rootNode;
      }//end if

      for (int i = npesActive - 2; i > 0; i--) {
        if (vtkDist[i] == vtkDist[i-1]) {
          vtkDist[i] = vtkDist[i+1];
        }//end if
      }//end for
    }//end if npes > 1

    unsigned int *part = NULL;
    if(!nodes.empty()) {
      part = new unsigned int[nodes.size()];
    }

    if (npesActive > 1) {
      unsigned int pCnt=0;
      for (unsigned int i=0; i< nodes.size(); i++) {
#ifdef __DEBUG_DA__
        assert(pCnt < npesActive);
#endif
        if ( (nodes[i] >= vtkDist[pCnt]) && ( (pCnt == (npesActive - 1)) ||
              ( nodes[i] < vtkDist[pCnt+1] ) || (vtkDist[pCnt+1] == rootNode) ) ) {
          part[i] = pCnt;
        } else {
          while ( (pCnt < (npesActive -1)) && (nodes[i] >= vtkDist[pCnt+1])
              && (vtkDist[pCnt+1] != rootNode)  ) {
            pCnt++;
          }//end while
          part[i] = pCnt;
        }//end if-else
      }//end for i
    }//end if np>1

    vtkDist.clear();
    //_________________________________________________________________________
    // Now the partitions should be contiguous since the two lists are globally
    // sorted ... and it's simply a shift between the two.
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //compute the total number of nodes being sent to each proc ...
    for (int i = 0; i < npesActive; i++) {
      sendCnt[i]=0;
      recvCnt[i]=0;
    }

    if (npesActive > 1) {
      for (unsigned int i=0; i<nodes.size(); i++) {
#ifdef __DEBUG_DA__
        assert(part[i] < npesActive);
#endif
        sendCnt[part[i]]++;
      }//end for i
    } else {
      sendCnt[0] += (nodes.size());
    }//end if-else

    if(part) {
      delete [] part;
      part = NULL;
    }

    // communicate with other procs how many you shall be sending and get how
    // many to recieve from whom.

    par::Mpi_Alltoall<int>( sendCnt, recvCnt, 1, commActive);

    unsigned int totalRecv = 0;
    for (unsigned int i = 0; i < npesActive; i++) {
      totalRecv += recvCnt[i];
    }//end for i

    sendOffsets[0] = 0;
    recvOffsets[0] = 0;

    // compute offsets ...
    for (int i=1; i < npesActive; i++) {
      sendOffsets[i] = sendOffsets[i-1] + sendCnt[i-1];
      recvOffsets[i] = recvOffsets[i-1] + recvCnt[i-1];
    }//end for i

    // Allocate for new array ...
    std::vector<ot::TreeNode > newNodes(totalRecv);

    // perform All2Allv
    ot::TreeNode* nodesPtr = NULL;
    ot::TreeNode* newNodesPtr = NULL;
    if(!nodes.empty()) {
      nodesPtr = &(*(nodes.begin()));
    }
    if(!newNodes.empty()) {
      newNodesPtr = &(*(newNodes.begin()));
    }
    par::Mpi_Alltoallv_sparse<ot::TreeNode>( nodesPtr, sendCnt, sendOffsets,
        newNodesPtr, recvCnt, recvOffsets, commActive);

    // reset the pointer ...
    nodes = newNodes;

    // clean up ...
    delete [] sendCnt;
    sendCnt = NULL;

    delete [] recvCnt;
    recvCnt = NULL;

    delete [] sendOffsets;
    sendOffsets = NULL;

    delete [] recvOffsets;
    recvOffsets = NULL;

    newNodes.clear();

    PROF_BLKPART3_END
  } // end blockPart

  //Assumptions: Blocks are sorted and unique and linear, nodes (Volume) is sorted and unique and linear.
  //Every element in nodes is a decendant of or equal to some block.
  //Some Blocks may be empty.
  //Assumes nodes is balanced.
  void pickGhostCandidates(const std::vector<ot::TreeNode> & blocks, const std::vector<ot::TreeNode> &nodes,
      std::vector<ot::TreeNode>& res, unsigned int dim, unsigned int maxDepth) {	

    PROF_PICK_GHOSTS_BEGIN 
      std::vector<unsigned int> firstLayer;
    //1. First Pick true inter-processor boundary octants 
    //(any octant that touches the inter-processor bdy.)

    pickInterProcessorBoundaryNodes(nodes, firstLayer, blocks[0], blocks[blocks.size() - 1]);

    std::vector<TreeNode> secondLayerBlocks(26*firstLayer.size());
    unsigned int secondLayerCount = 0;
    for(unsigned int i = 0; i < firstLayer.size(); i++) {
#ifdef __DEBUG_DA__
      assert( firstLayer[i] < nodes.size() );
#endif
      std::vector<TreeNode> myNh = nodes[firstLayer[i]].getAllNeighbours();
      for(unsigned int j = 0; j < myNh.size(); j++) {
        if( (myNh[j] >= blocks[0]) && (myNh[j] <= blocks[blocks.size()-1].getDLD()) ) {
          //myNh[j] lies in the domain controlled by my processor so add it.
          secondLayerBlocks[secondLayerCount] = myNh[j];
          secondLayerCount++;
        }
      }//end for j
    }//end for i

    secondLayerBlocks.resize(secondLayerCount);

    seq::makeVectorUnique<ot::TreeNode>(secondLayerBlocks, false);

    //2 Generate Second Ring Candidates depending on the way the first layer touches the block boundaries.
    //secondRing is sorted, unique and exists in nodes. (linearity follows directly)
    std::vector<unsigned int> secondRing;

    addSecondRing(nodes,firstLayer,blocks,secondRing);

    //2.c. Compare with the Second Ring Candidates and eliminate unwanted octants.
    //secondLayer is sorted, unique and not linear.
    //secondRing is sorted, unique and linear.
    unsigned int layerCnt = 0;
    unsigned int ringCnt = 0;
    res.clear();     
    while( (layerCnt < secondLayerBlocks.size()) && (ringCnt < secondRing.size()) ) {
      if( secondLayerBlocks[layerCnt] < nodes[secondRing[ringCnt]] ) {		
        if(secondLayerBlocks[layerCnt].isAncestor(nodes[secondRing[ringCnt]])) {
          unsigned int tmpCnt = ringCnt;
          while ( (tmpCnt < secondRing.size()) && 
              (secondLayerBlocks[layerCnt].isAncestor(nodes[secondRing[tmpCnt]]) ) ) {
            if( secondLayerBlocks[layerCnt] == nodes[secondRing[tmpCnt]].getParent() ) {
              res.push_back(nodes[secondRing[tmpCnt]]);
            }
            tmpCnt++;
          }
        }
        layerCnt++;
      }else {
        ringCnt++;
      }
    }

    //3. Add the first layer and makeVectorUnique.
    for(unsigned int i=0; i < firstLayer.size(); i++) {
      res.push_back(nodes[firstLayer[i]]);
    }
    seq::makeVectorUnique<ot::TreeNode>(res,false);

    firstLayer.clear();
    secondRing.clear();
    secondLayerBlocks.clear();

    PROF_PICK_GHOSTS_END 
  }// end fn.

  //secondRing is sorted, unique and exists in nodes. (linearity follows directly)
  //blocks is sorted, unique and linear.
  //firstLayer is sorted, unique and linear.
  int addSecondRing(const std::vector<ot::TreeNode> & nodes, const std::vector<unsigned int> & firstLayer,
      const std::vector<ot::TreeNode> & blocks , std::vector<unsigned int> & secondRing) {
    unsigned int blkCnt =0;
    unsigned int firstCnt = 0;
    std::vector<ot::TreeNode> tmpLayer;

    while( (blkCnt < blocks.size()) && (firstCnt < firstLayer.size()) ) {
      if( blocks[blkCnt] <= nodes[firstLayer[firstCnt]] ) {
        if( (blocks[blkCnt].isAncestor(nodes[firstLayer[firstCnt]])) ||
            (blocks[blkCnt] == nodes[firstLayer[firstCnt]]) ) {
          unsigned char flags;
          bool isBnd = nodes[firstLayer[firstCnt]].isBoundaryOctant(blocks[blkCnt],
              ((ot::TreeNode::NEGATIVE) | (ot::TreeNode::POSITIVE)),&flags);						

          if(isBnd) {
            bool negX = (flags & 1);
            bool negY = (flags & 2);
            bool negZ = (flags & 4);
            bool posX = (flags & 32);
            bool posY = (flags & 64);
            bool posZ = (flags & 128);
            bool isX = (negX || posX);
            bool isY = (negY || posY);
            bool isZ = (negZ || posZ);

            //Faces...
            if(isX) {
              if(negX) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getRight());					
              }
              if(posX) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getLeft());
              }		
            }
            if(isY) {
              if(negY) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getBack()); 
              }
              if(posY) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getFront());
              }
            }
            if(isZ) {
              if(negZ) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getTop());
              }
              if(posZ) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getBottom());
              }
            }

            //Edges...						
            if(isX && isY) {						
              if(negX && negY) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getRightBack());						
              }
              if(negX && posY) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getRightFront());
              }
              if(posX && negY) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getLeftBack());
              }			
              if(posX && posY) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getLeftFront());
              }
            }

            if(isZ && isY) {									
              if(negZ && negY) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getTopBack());
              }
              if(negZ && posY) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getTopFront());					
              }
              if(posZ && negY) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getBottomBack());
              }			
              if(posZ && posY) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getBottomFront());
              }
            }

            if(isX && isZ) {
              if(negX && negZ) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getTopRight());
              }
              if(negX && posZ) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getBottomRight());
              }
              if(posX && negZ) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getTopLeft());
              }			
              if(posX && posZ) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getBottomLeft()) ; 
              }
            }

            //Corners...				 				 
            if(isX && isY && isZ) {
              if( negX && negY && negZ ) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getTopRightBack());
              }
              if( posX && negY && negZ ) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getTopLeftBack());
              }
              if( posX && negY && posZ ) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getBottomLeftBack());
              }
              if( negX && negY && posZ ) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getBottomRightBack());
              }
              if( posX && posY && negZ ) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getTopLeftFront());
              }
              if( negX && posY && negZ ) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getTopRightFront());
              }
              if( posX && posY && posZ ) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getBottomLeftFront());
              }
              if( negX && posY && posZ ) {
                tmpLayer.push_back(nodes[firstLayer[firstCnt]].getBottomRightFront());
              }
            }
          }
          //The same block must be compared with other possible decendants.
          firstCnt++;
        }else {
          blkCnt++;
        }
      }else {
        firstCnt++;
      }
    }	

    seq::makeVectorUnique<ot::TreeNode>(tmpLayer,false);

    //tmpLayer is sorted and unique, but not linear.
    unsigned int tmpCnt = 0;
    unsigned int nodeCnt = 0;
    while( (tmpCnt < tmpLayer.size()) && (nodeCnt < nodes.size()) ) {
      if( nodes[nodeCnt] <= tmpLayer[tmpCnt] ) {
        if(nodes[nodeCnt] == tmpLayer[tmpCnt]) {
          secondRing.push_back(nodeCnt);
        }
        nodeCnt++;					
      }else {
        tmpCnt++;
      }
    }
    tmpLayer.clear();

    return 1;
  }//end fn.


} // end namespace ot


