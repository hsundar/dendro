
/**
  @file odaFactory.C
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  */

#include "oda.h"
#include "parUtils.h"
#include "colors.h"
#include "testUtils.h"
#include "dendro.h"

#ifdef __DEBUG__
#ifndef __DEBUG_DA__
#define __DEBUG_DA__
#endif
#endif

#ifdef __DEBUG_DA__
#ifndef __DEBUG_DA_PUBLIC__
#define __DEBUG_DA_PUBLIC__
#endif
#endif

#ifdef __DEBUG_DA_PUBLIC__
#ifndef __MEASURE_DA__
#define __MEASURE_DA__
#endif
#endif

namespace ot {

#define RESET_DA_BLOCK {\
  /*Initialize Member Data*/\
  m_ucpOctLevels = NULL;\
  m_bComputedLocalToGlobal = false;\
  m_bComputedLocalToGlobalElems = false;\
  m_tnBlocks.clear();\
  m_tnMinAllBlocks.clear();\
  m_uiNlist.clear();\
  m_uiNlistPtr = NULL;\
  m_dilpLocalToGlobal = NULL;\
  m_dilpLocalToGlobalElems = NULL;\
  m_ucpLutRemainders.clear();\
  m_ucpLutRemaindersPtr = NULL;\
  m_ucpSortOrders.clear();\
  m_ucpSortOrdersPtr = NULL;\
  m_uspLutQuotients.clear();\
  m_uspLutQuotientsPtr = NULL;\
  m_ucpLutMasks.clear();\
  m_ucpLutMasksPtr = NULL;\
  m_ucpPreGhostConnectivity.clear();\
  m_ucpPreGhostConnectivityPtr = NULL;\
  m_ptsPreGhostOffsets.clear();\
  m_ptsPreGhostOffsetsPtr = NULL;\
  PreGhostAnchors.clear();\
  m_uiQuotientCounter = 0;\
  m_uiPreGhostQuotientCnt = 0;\
  m_uiElementQuotient = 0;\
  m_uiIndependentElementQuotient = 0;\
  m_uiNodeSize = 0;\
  m_uiBoundaryNodeSize = 0;\
  m_uiElementSize = 0;\
  m_uiIndependentElementSize = 0;\
  m_uiPreGhostElementSize = 0;\
  m_uiPreGhostNodeSize = 0;\
  m_uiPreGhostBoundaryNodeSize = 0;\
  m_uiPostGhostNodeSize = 0;\
  m_uiLocalBufferSize = 0;\
  m_uiElementBegin = 0;\
  m_uiElementEnd = 0;\
  m_uiPostGhostBegin = 0;\
  m_uiIndependentElementBegin = 0;\
  m_uiIndependentElementEnd = 0;\
  m_uiCurrent = 0;\
  m_uiDimension = 0;\
  m_uiMaxDepth = 0;\
  m_uipScatterMap.clear();\
  m_uipSendProcs.clear();\
  m_uipSendOffsets.clear();\
  m_uipSendCounts.clear();\
  m_uipElemScatterMap.clear();\
  m_uipElemSendProcs.clear();\
  m_uipElemSendOffsets.clear();\
  m_uipElemSendCounts.clear();\
  m_uipRecvProcs.clear();\
  m_uipRecvOffsets.clear();\
  m_uipRecvCounts.clear();\
  m_uipElemRecvProcs.clear();\
  m_uipElemRecvOffsets.clear();\
  m_uipElemRecvCounts.clear();\
  m_mpiContexts.clear();\
  m_bCompressLut = compressLut;\
  m_uiCommTag = 1;\
  m_mpiCommAll = comm;\
  MPI_Comm_size(m_mpiCommAll,&m_iNpesAll);\
  MPI_Comm_rank(m_mpiCommAll,&m_iRankAll);\
}

void DA::DA_FactoryPart0(std::vector<ot::TreeNode>& in, MPI_Comm comm,
    MPI_Comm activeInputComm, bool compressLut, bool* iAmActive) {
  RESET_DA_BLOCK
    m_uiInputSize = static_cast<unsigned int>(in.size());

  //The default is NULL. 
  if(iAmActive != NULL) {
    (*iAmActive) = (!(in.empty())); 
    m_bIamActive = (*iAmActive);    
  }else {
    //Don't care for active state.
    m_bIamActive = (!(in.empty()));
  }

#ifdef __DEBUG_DA_PUBLIC__
  MPI_Comm expectedActiveInputComm;
  int commCompareResult;
  par::splitComm2way(m_bIamActive, &expectedActiveInputComm, comm);
  if(m_bIamActive) {
    MPI_Comm_compare(activeInputComm, expectedActiveInputComm, &commCompareResult);
    assert( (commCompareResult == MPI_CONGRUENT) || (commCompareResult == MPI_IDENT) );
  }
#endif

  m_mpiCommActive = activeInputComm;

  if(m_bIamActive) {
    MPI_Comm_size(m_mpiCommActive, &m_iNpesActive);
    MPI_Comm_rank(m_mpiCommActive, &m_iRankActive);
  } else {
    m_iNpesActive = 0;
    m_iRankActive = 0;
  }

  if(!m_bIamActive) {
    return;
  } else {
    assert(!in.empty());
    m_uiDimension = in[0].getDim();
    m_uiMaxDepth  = in[0].getMaxDepth();
  }
}//end function

void DA::DA_FactoryPart1(std::vector<ot::TreeNode>& in) {
#ifdef __PROF_WITH_BARRIER__
  MPI_Barrier(m_mpiCommActive);
#endif
  PROF_BUILD_DA_STAGE1_BEGIN

    assert(!in.empty());

  // first generate the boundary nodes ...
  std::vector<ot::TreeNode> positiveBoundaryOctants;
  //Assumption: in is globally sorted to begin with. 
  //Guarantee: in remains globally sorted in the end.
  //positiveBoundaryOctants need not be globally sorted, in fact I don't think
  //it will be sorted even on 1 processor.
  addBoundaryNodesType1(in, positiveBoundaryOctants, m_uiDimension, m_uiMaxDepth);

  // Update the maxDepth ...
  m_uiMaxDepth = m_uiMaxDepth + 1;

  //Most processors will not add any positive boundaries. So, there is no need
  //to unnecessarily distribute positive boundaries on all procs just for
  //sorting. So only the few processors touching the positive boundary will
  //participate in the parallel sort
  MPI_Comm bdyComm;
  par::splitComm2way((positiveBoundaryOctants.empty()), &bdyComm, m_mpiCommActive);

  if(!(positiveBoundaryOctants.empty())) {
    //Call Sample Sort  
    std::vector<ot::TreeNode > tmpVecTN;
    par::sampleSort<ot::TreeNode>(positiveBoundaryOctants, tmpVecTN, bdyComm);
    positiveBoundaryOctants = tmpVecTN;
    tmpVecTN.clear();
  }

  par::concatenate<ot::TreeNode>(in, positiveBoundaryOctants, m_mpiCommActive);
  positiveBoundaryOctants.clear();

  PROF_BUILD_DA_STAGE1_END
}//end function

void DA::DA_FactoryPart2(std::vector<ot::TreeNode>& in) {
#ifdef __PROF_WITH_BARRIER__
  MPI_Barrier(m_mpiCommActive);
#endif
  PROF_BUILD_DA_STAGE2_BEGIN

    assert(!in.empty());
  //Marks regular nodes. MUST be called before DA_blockPartStage2
  //Assumption: in is globally sorted at this point, including positive
  //boundaries 

  flagNodesType3(in, m_mpiCommActive);

  PROF_BUILD_DA_STAGE2_END
}//end function

void DA::DA_FactoryPart3(std::vector<ot::TreeNode>& in, MPI_Comm comm, bool compressLut,
    const std::vector<ot::TreeNode>* blocksPtr, bool* iAmActive) {
#ifdef __PROF_WITH_BARRIER__
  MPI_Barrier(m_mpiCommActive);
#endif
  PROF_BUILD_DA_STAGE3_BEGIN

    assert(!in.empty());
  //  1. First repartition  

  DendroIntL localSizeBefore = in.size();
  DendroIntL globalSizeBefore; 
  par::Mpi_Allreduce<DendroIntL>(&localSizeBefore, &globalSizeBefore, 1, MPI_SUM, m_mpiCommActive);

  //Partition in and create blocks (blocks must be globally sorted).
  std::vector<ot::TreeNode> blocks;
  if(blocksPtr == NULL) {

    //min grain size = 1000
    const DendroIntL THOUSAND = 1000;
    if (globalSizeBefore < (THOUSAND*m_iNpesActive)) {
      int splittingSize = (globalSizeBefore/THOUSAND); 
      if(splittingSize == 0) {
        splittingSize = 1; 
      }

      unsigned int avgLoad = (globalSizeBefore / splittingSize);
      int leftOvers = (globalSizeBefore % splittingSize);

      std::vector<TreeNode> tmpIn;
      if(m_iRankActive >= splittingSize) {
        par::scatterValues<ot::TreeNode>(in, tmpIn, 0, m_mpiCommActive);
      }else if(m_iRankActive < leftOvers) {
        par::scatterValues<ot::TreeNode>(in, tmpIn, (avgLoad+1), m_mpiCommActive);
      }else {
        par::scatterValues<ot::TreeNode>(in, tmpIn, avgLoad, m_mpiCommActive);
      }
      in = tmpIn;
      tmpIn.clear();

      MPI_Comm newComm;
      par::splitCommUsingSplittingRank(splittingSize, &newComm, m_mpiCommActive);
      m_mpiCommActive = newComm;
      MPI_Comm_size(m_mpiCommActive,&m_iNpesActive);
      MPI_Comm_rank(m_mpiCommActive,&m_iRankActive);

#ifndef __SILENT_MODE__
      if(!m_iRankActive) {
        std::cout<<" input to DA constructor is small("<<globalSizeBefore
          <<") npes = "<<m_iNpesActive<<" splitting comm."<<std::endl;
      }
#endif

      m_bIamActive = (!in.empty());   
      if(iAmActive != NULL) {	
        //Want the active state returned
        (*iAmActive) = m_bIamActive;
      }
      if(!m_bIamActive) {
        RESET_DA_BLOCK
          PROF_BUILD_DA_STAGE3_END
          return;
      }
    }//end check if total size is too small

    PROF_DA_BPART1_BEGIN

      ot::blockPartStage1(in, blocks, m_uiDimension, m_uiMaxDepth, m_mpiCommActive);    

    PROF_DA_BPART1_END

      PROF_DA_BPART2_BEGIN

      DA_blockPartStage2(in, blocks, m_uiDimension, m_uiMaxDepth, m_mpiCommActive);

    PROF_DA_BPART2_END
  } else {
    blocks = *blocksPtr;
  }

  PROF_DA_BPART3_BEGIN

    DA_blockPartStage3(in, blocks, m_tnMinAllBlocks, m_uiDimension, m_uiMaxDepth, m_mpiCommActive);

  PROF_DA_BPART3_END

    //Store Blocks. 
    m_tnBlocks = blocks;

  //we must split comm if new empty procs were created in blockPartStage2
  //empty procs are created using partW and 0 wts and so all empty procs should
  //be contiguous and only the last few procs should be empty
  if(m_tnMinAllBlocks.size() != m_iNpesActive) {
    m_bIamActive = (!in.empty());   
    if(iAmActive != NULL) {	
      //Want the active state returned
      (*iAmActive) = m_bIamActive;
    }
    assert( m_bIamActive == (m_iRankActive < m_tnMinAllBlocks.size()) );

    MPI_Comm tmpComm;
    par::splitCommUsingSplittingRank(static_cast<int>(m_tnMinAllBlocks.size()),
        &tmpComm, m_mpiCommActive);
    m_mpiCommActive = tmpComm;

#ifndef __SILENT_MODE__
    if(!m_iRankActive) {
      std::cout<<" splitComm: "<<m_iNpesActive<<" -> "<<
        (m_tnMinAllBlocks.size())<<" after bPart in DA constructor. "<<std::endl;
    }
#endif

    if(m_bIamActive) {
      MPI_Comm_size(m_mpiCommActive,&m_iNpesActive);
      MPI_Comm_rank(m_mpiCommActive,&m_iRankActive);
    } else {
      m_iNpesActive = 0;
      m_iRankActive = 0;
    }

    if(!m_bIamActive) {
      RESET_DA_BLOCK
        PROF_BUILD_DA_STAGE3_END
        return;
    }
  }//end if new empty procs after bPart

  PROF_BUILD_DA_STAGE3_END
#ifdef __PROF_WITH_BARRIER__
    MPI_Barrier(m_mpiCommActive);
#endif
  PROF_BUILD_DA_STAGE4_BEGIN

    // Set the Local offset
    assert(!in.empty());
  assert(m_tnMinAllBlocks.size() == m_iNpesActive);

  m_ptOffset = in[0].getAnchor();

#ifdef __DEBUG_DA_PUBLIC__
  MPI_Barrier(m_mpiCommActive);
  //Check that block Part did the right job.
  //1. Ensure that in is globally sorted and unique.
  assert(par::test::isUniqueAndSorted(in, m_mpiCommActive));

  MPI_Barrier(m_mpiCommActive);
  //2. Ensure that the global size before and after block part has not
  //changed. So we don't miss out octants
  DendroIntL localSizeAfter = in.size();
  DendroIntL globalSizeAfter; 
  par::Mpi_Allreduce<DendroIntL>(&localSizeAfter, &globalSizeAfter, 1,
      MPI_SUM, m_mpiCommActive);

  assert(globalSizeAfter == globalSizeBefore);

  MPI_Barrier(m_mpiCommActive);

  //3. Locally check that if an anchor is hanging then the parent's first
  // child exists.
  for(unsigned int i=0; i< in.size(); i++) {
    if( !(in[i].getFlag() & ot::TreeNode::NODE) ) {
      ot::TreeNode anchorOfParent = in[i].getParent().getDFD();
      unsigned int idxOfParent;
      bool foundAnchorOfParent = seq::maxLowerBound<ot::TreeNode>(in, anchorOfParent, idxOfParent, NULL, NULL);
      if(!foundAnchorOfParent) {
        std::cout<<m_iRankActive<<" failing for: "<<in[i]<<std::endl<<
          "Expected to find: "<<
          anchorOfParent.getAncestor(in[i].getLevel())<<std::endl;
        assert(false);
      }
      if( anchorOfParent.getAnchor() != in[idxOfParent].getAnchor() ) {
        std::cout<<m_iRankActive<<" failing for: "<<in[i]<<std::endl<<
          "Expected to find: "<<
          anchorOfParent.getAncestor(in[i].getLevel())<<std::endl;
        assert(false);
      }
      if( in[idxOfParent].getParent() != in[i].getParent() ) {
        std::cout<<m_iRankActive<<" failing for: "<<in[i]<<std::endl<<
          "Expected to find: "<<
          anchorOfParent.getAncestor(in[i].getLevel())<<std::endl;
        assert(false);
      }
    }
  }

  MPI_Barrier(m_mpiCommActive);

  if(!m_iRankActive) {
    std::cout<<"Partition is correct."<<std::endl;
  }

  MPI_Barrier(m_mpiCommActive);
#endif

  //  2. Obtain all ghost octants.

#ifdef __DEBUG_DA_PUBLIC__
  MPI_Barrier(m_mpiCommActive);
  if(!m_iRankActive) {
    std::cout<<"Pick Inter-Processor Boundary Nodes "<<std::endl;
  }
  MPI_Barrier(m_mpiCommActive);
#endif

  // first let's pick boundary nodes on each proc.
  std::vector<ot::TreeNode> allBoundaryLeaves;

  pickInterProcessorBoundaryNodes(in, allBoundaryLeaves, blocks[0], blocks[blocks.size() - 1]);

  ot::TreeNode myFirstOctant = in[0];
  ot::TreeNode myLastOctant = in[in.size() - 1];

  includeSiblingsOfBoundary(allBoundaryLeaves, myFirstOctant, myLastOctant);

  PROF_BUILD_DA_STAGE4_END
#ifdef __PROF_WITH_BARRIER__
    MPI_Barrier(m_mpiCommActive);
#endif
  PROF_BUILD_DA_STAGE5_BEGIN

#ifdef __DEBUG_DA_PUBLIC__
    MPI_Barrier(m_mpiCommActive);
  if(!m_iRankActive) {
    std::cout<<"Selected Extra Ghost Candidates. Now Checking "<<std::endl;      
  }

  assert(seq::test::isUniqueAndSorted<ot::TreeNode>(allBoundaryLeaves));

  unsigned int allBoundaryLeavesSize = allBoundaryLeaves.size();
  unsigned int maxBndSz;
  par::Mpi_Reduce<unsigned int>(&allBoundaryLeavesSize,&maxBndSz,1,MPI_MAX,0,m_mpiCommActive);

  if(!m_iRankActive) {
    std::cout<<"Max Bnd Size:  "<<maxBndSz<<std::endl;      
  }

  MPI_Barrier(m_mpiCommActive);
#endif

  PROF_BUILD_DA_STAGE5_END
#ifdef __PROF_WITH_BARRIER__
    MPI_Barrier(m_mpiCommActive);
#endif
  PROF_BUILD_DA_STAGE6_BEGIN

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // 4. Loop through all local boundary nodes and determine which processors
    // need to be aware of those nodes. Create lists that shall be sent.
    //

    int *sendCnt = new int[m_iNpesActive];
  int *recvCnt = new int[m_iNpesActive];
  int *sendOffsets = new int[m_iNpesActive];
  int *recvOffsets = new int[m_iNpesActive];

  std::vector< std::vector<unsigned int> > sendNodes;

  prepareAprioriCommMessagesInDAtype2(in, allBoundaryLeaves, blocks, m_tnMinAllBlocks,
      m_iRankActive, m_iNpesActive, sendCnt, sendNodes);

  allBoundaryLeaves.clear();
  blocks.clear();

#ifdef __DEBUG_DA_PUBLIC__
  MPI_Barrier(m_mpiCommActive);
  if(!m_iRankActive) {
    std::cout<<"Picked Octants for Apriori Comm. "<<std::endl;
  }
  MPI_Barrier(m_mpiCommActive);
#endif

  // 5. Actual send/recv. to exchange nodes.
  //
  // 5a.
  // Now do an All2All to get numKeysRecv
  par::Mpi_Alltoall<int>( sendCnt, recvCnt, 1, m_mpiCommActive);

  // 5b. Concatenate all nodes into one single Carray ...
  unsigned int totalSend=0;
  unsigned int totalRecv=0;
  for (unsigned int i=0; i < m_iNpesActive; i++) {
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
  for (int i=1; i < m_iNpesActive; i++) {
    sendOffsets[i] = sendOffsets[i-1] + sendCnt[i-1];
    recvOffsets[i] = recvOffsets[i-1] + recvCnt[i-1];
  }

  int myOff = recvOffsets[m_iRankActive];

#ifdef __DEBUG_DA_PUBLIC__
  assert(sendCnt[m_iRankActive] == 0);
#endif

  for (int i=0; i < m_iNpesActive; i++) {
    for (unsigned int j=0; j<sendCnt[i]; j++) {
      sendK[sendOffsets[i] + j] = in[sendNodes[i][j]];
      m_uipScatterMap.push_back(sendNodes[i][j] + myOff);
    }
  }

#ifdef __DEBUG_DA_PUBLIC__
  MPI_Barrier(m_mpiCommActive);
  if(!m_iRankActive) {
    std::cout<<"Built Primary ScatterMap. "<<std::endl;
  }
  MPI_Barrier(m_mpiCommActive);
#endif

  // 5c. Perform SendRecv to send and receive all keys ...

  ot::TreeNode* sendKptr = NULL;
  ot::TreeNode* recvKptr = NULL;
  if(!sendK.empty()) {
    sendKptr = &(*(sendK.begin()));
  }
  if(!recvK.empty()) {
    recvKptr = &(*(recvK.begin()));
  }

  par::Mpi_Alltoallv_sparse<ot::TreeNode>( sendKptr, sendCnt, sendOffsets, 
      recvKptr, recvCnt, recvOffsets, m_mpiCommActive);

  sendK.clear();
  for (unsigned int i=0; i < m_iNpesActive; i++) {
    sendNodes[i].clear();
  }
  sendNodes.clear();

  // Let's store the counts for later synchronizations  ...
  //The offsets will be computed just before compression, After primary and
  //secondary scatterMaps are merged.
  for ( unsigned int i=0; i < m_iNpesActive; i++) {
    if ( sendCnt[i] ) {
      m_uipSendProcs.push_back(i);
      m_uipSendCounts.push_back(sendCnt[i]);
#ifdef __DEBUG_DA_PUBLIC__
      assert(i != m_iRankActive);
#endif
    }
    if ( recvCnt[i] ) {
      m_uipRecvProcs.push_back(i);
      m_uipRecvCounts.push_back(recvCnt[i]);
#ifdef __DEBUG_DA_PUBLIC__
      assert(i != m_iRankActive);
#endif
    }
  }

  //Free some memory...
  delete [] sendCnt;
  delete [] recvCnt;
  delete [] sendOffsets;
  delete [] recvOffsets;

  PROF_BUILD_DA_STAGE6_END
#ifdef __PROF_WITH_BARRIER__
    MPI_Barrier(m_mpiCommActive);
#endif
  PROF_BUILD_DA_STAGE7_BEGIN

    //Now that we have the ghost nodes (in recvK), we can merge this with the
    //input nodes, and all tag all of these as elements. While merging we shall
    //also set the right offsets into the global vector.

#ifdef __DEBUG_DA_PUBLIC__        
    MPI_Barrier(m_mpiCommActive);
  //std::cout << "Now Merging recv ghosts with own octants, recvK size is " << recvK.size() << std::endl;
  //std::cout << "In size is " << in.size() << std::endl;
  assert(seq::test::isSorted<ot::TreeNode> (recvK));
  MPI_Barrier(m_mpiCommActive);
#endif        

  // Lets store offsets demarcating the global domain spanned by this
  // processor.

  ot::TreeNode globalMin, globalMax;
  globalMin = in[0];


  // Since recvK is sorted, and so are the 'global' nodes, we can merge
  // effectively. We also mark them as being elements.  

  //Merge (In-place) recieved octants with local octants....
  std::vector < ot::TreeNode > localOcts(recvK.size() + in.size() );

#ifdef __DEBUG_DA_PUBLIC__
  MPI_Barrier(m_mpiCommActive);
  if(!m_iRankActive) {
    std::cout<<"Merging ghosts and own octants into the buffer. "<<std::endl;
  }
  MPI_Barrier(m_mpiCommActive);
#endif


  //Note this is an inplace insertion. This is done to preserve the sorted
  //order. Each recvK[i] is independently sorted. Also, since the intial
  //blocks were globally sorted and since the elements on each processor are
  //only decendants of these blocks, hence recvK[i][j] < recvK[i+][k] for all i,j
  //and k.
  // PreGhost ....
  for (int i=0; i<myOff; i++) {
    localOcts[i] = recvK[i];
    if (!(localOcts[i].getFlag() & ot::TreeNode::BOUNDARY) ) {
      m_uiPreGhostElementSize++;
    }
  }

  //Mine ...
  m_uiElementBegin = myOff;
  m_uiElementEnd = myOff;
  for (int i=myOff; i<myOff+in.size(); i++) {
    localOcts[i] = in[i-myOff];
    if (!(localOcts[i].getFlag() & ot::TreeNode::BOUNDARY) ) {
      m_uiElementEnd++;
    } 
  }
  m_uiElementSize = (m_uiElementEnd - m_uiElementBegin);

  //PostGhosts....
  m_uiPostGhostBegin = myOff + static_cast<unsigned int>(in.size());
  for (int i=myOff+in.size(); i<localOcts.size(); i++) {
    localOcts[i] = recvK[i-in.size()];
  }

  // all indices should be set at this stage ...
  // free up some memory.
  in.clear();
  recvK.clear();

  m_uiLocalBufferSize = static_cast<unsigned int>(localOcts.size());

  PROF_BUILD_DA_STAGE7_END
#ifdef __PROF_WITH_BARRIER__
    MPI_Barrier(m_mpiCommActive);
#endif
  PROF_BUILD_DA_STAGE8_BEGIN

    // Now for the big step. Perform searches to build the look-up table.
    // LocalOcts will be modified inside buildNodeList.
    // All counters will be reset inside the function.
    // The list will continue to remain sorted.


#ifdef __DEBUG_DA_PUBLIC__
    MPI_Barrier(m_mpiCommActive);
  assert(seq::test::isUniqueAndSorted(localOcts));
  if(!m_iRankActive) {
    std::cout<<" Entering BuildNodeList. "<<std::endl;
  }    
  MPI_Barrier(m_mpiCommActive);
#endif

  buildNodeList(localOcts);

#ifdef __DEBUG_DA_PUBLIC__
  MPI_Barrier(m_mpiCommActive);
  if(!m_iRankActive) {
    std::cout<<" Leaving BuildNodeList. "<<std::endl;
  }
  assert(seq::test::isUniqueAndSorted(localOcts));
  MPI_Barrier(m_mpiCommActive);
#endif

  PROF_BUILD_DA_STAGE8_END
#ifdef __PROF_WITH_BARRIER__
    MPI_Barrier(m_mpiCommActive);
#endif
  PROF_BUILD_DA_STAGE9_BEGIN

    // Set the global offset 
    m_ptGhostedOffset = localOcts[0].getAnchor();

  //Compute pre-ghost offsets
#ifdef __DEBUG_DA_PUBLIC__
  if ( m_uiElementBegin) {
    PreGhostAnchors.push_back( localOcts[0].getAnchor());
  }
#endif

  for (unsigned int i = 1; i < m_uiElementBegin; i++) {
    unsigned char touchConfig = getTouchConfig(localOcts[i-1],
        localOcts[i], m_uiMaxDepth);
    m_ucpPreGhostConnectivity.push_back(touchConfig);
    if (!touchConfig) {
      m_ptsPreGhostOffsets.push_back( localOcts[i].getAnchor() );
    }
#ifdef __DEBUG_DA_PUBLIC__
    PreGhostAnchors.push_back( localOcts[i].getAnchor());
#endif
  }

  // Finally compress the octree, clean up and we are all set.

  // we simply retain the oct levels ...
  if(!localOcts.empty()) {
    m_ucpOctLevels = new unsigned char [localOcts.size()];
  } else {
    m_ucpOctLevels = NULL;
  }

  for (unsigned int i = 0; i < localOcts.size(); i++) {
    m_ucpOctLevels[i] = localOcts[i].getFlag();
  }

  // clean up ...
  localOcts.clear();

  // Set pointers ....
  if(!m_ucpLutRemainders.empty()) {
    m_ucpLutRemaindersPtr = &(*m_ucpLutRemainders.begin());
  } else {
    m_ucpLutRemaindersPtr = NULL;
  }

  if(!m_uspLutQuotients.empty()) {
    m_uspLutQuotientsPtr = &(*m_uspLutQuotients.begin());
  } else {
    m_uspLutQuotientsPtr = NULL;
  }

  if(!m_ucpLutMasks.empty()) {
    m_ucpLutMasksPtr = &(*m_ucpLutMasks.begin());
  } else {
    m_ucpLutMasksPtr = NULL;
  }

  if(!m_ucpSortOrders.empty()) {
    m_ucpSortOrdersPtr = &(*m_ucpSortOrders.begin());
  } else {
    m_ucpSortOrdersPtr = NULL;
  }

  if(!m_uiNlist.empty()) {
    m_uiNlistPtr = &(*m_uiNlist.begin());
  } else {
    m_uiNlistPtr = NULL;
  }

  PROF_BUILD_DA_STAGE9_END

#ifdef __DEBUG_DA_PUBLIC__
    //Check Loops
    MPI_Barrier(m_mpiCommActive);
  for( init<ot::DA_FLAGS::WRITABLE>(); curr() < end<ot::DA_FLAGS::WRITABLE>();
      next<ot::DA_FLAGS::WRITABLE>() ) {
    assert(curr() >= m_uiElementBegin);
    assert(curr() < m_uiElementEnd);
    unsigned int indices[8];
    getNodeIndices(indices);
    for(unsigned int vtxId = 0; vtxId < 8; vtxId++) {
      assert(indices[vtxId] >= m_uiElementBegin);
      assert(indices[vtxId] < m_uiLocalBufferSize);		
    }
  }
  MPI_Barrier(m_mpiCommActive);
  if(!m_iRankActive) {
    std::cout<<" Finished checking WRITABLE."<<std::endl;
  }
  MPI_Barrier(m_mpiCommActive);

  for( init<ot::DA_FLAGS::INDEPENDENT>(); curr() < end<ot::DA_FLAGS::INDEPENDENT>();
      next<ot::DA_FLAGS::INDEPENDENT>() ) {
    assert(curr() >= m_uiElementBegin);
    assert(curr() < m_uiElementEnd);
    unsigned int indices[8];
    getNodeIndices(indices);
    for(unsigned int vtxId = 0; vtxId < 8; vtxId++) {
      assert(indices[vtxId] >= m_uiElementBegin);	
      assert(indices[vtxId] < m_uiPostGhostBegin);		
    }
  }

  MPI_Barrier(m_mpiCommActive);

  if(!m_iRankActive) {
    std::cout<<" Finished checking INDEPENDENT."<<std::endl;
  }

  MPI_Barrier(m_mpiCommActive);

  for( init<ot::DA_FLAGS::DEPENDENT>(); curr() < end<ot::DA_FLAGS::DEPENDENT>();
      next<ot::DA_FLAGS::DEPENDENT>() ) {
    assert(curr() < m_uiElementEnd);
    unsigned int indices[8];
    getNodeIndices(indices);
    unsigned int numLocal = 0;
    for(unsigned int vtxId = 0; vtxId < 8; vtxId++) {
      assert(indices[vtxId] < m_uiLocalBufferSize);
      if( (indices[vtxId] >= m_uiElementBegin) &&  (indices[vtxId] < m_uiPostGhostBegin) ) {
        numLocal++;
      }		
    }
    assert(numLocal > 0);
    assert(numLocal < 8);
  }

  MPI_Barrier(m_mpiCommActive);

  if(!m_iRankActive) {
    std::cout<<" Finished checking DEPENDENT."<<std::endl;
  }
  MPI_Barrier(m_mpiCommActive);

  for( init<ot::DA_FLAGS::W_DEPENDENT>(); curr() < end<ot::DA_FLAGS::W_DEPENDENT>();
      next<ot::DA_FLAGS::W_DEPENDENT>() ) {
    assert(curr() >= m_uiElementBegin);
    assert(curr() < m_uiElementEnd);
    unsigned int indices[8];
    getNodeIndices(indices);
    unsigned int numLocal = 0;
    for(unsigned int vtxId = 0; vtxId < 8; vtxId++) {
      assert(indices[vtxId] >= m_uiElementBegin);	
      assert(indices[vtxId] < m_uiLocalBufferSize);		
      if( (indices[vtxId] >= m_uiElementBegin) &&  (indices[vtxId] < m_uiPostGhostBegin) ) {
        numLocal++;
      }		
    }
    assert(numLocal > 0);
    assert(numLocal < 8);
  }

  MPI_Barrier(m_mpiCommActive);

  if(!m_iRankActive) {
    std::cout<<" Finished checking W_DEPENDENT."<<std::endl;
  }

  MPI_Barrier(m_mpiCommActive);

  for( init<ot::DA_FLAGS::ALL>(); 
      curr() < end<ot::DA_FLAGS::ALL>();
      next<ot::DA_FLAGS::ALL>() ) {
    assert(curr() < m_uiElementEnd);
    unsigned int indices[8];
    getNodeIndices(indices);
    for(unsigned int vtxId = 0; vtxId < 8; vtxId++) {
      assert(indices[vtxId] < m_uiLocalBufferSize);		
    }
  }

  MPI_Barrier(m_mpiCommActive);

  if(!m_iRankActive) {
    std::cout<<" Finished checking ALL."<<std::endl;
  }

  MPI_Barrier(m_mpiCommActive);

#endif

#ifdef __DEBUG_DA_PUBLIC__
  MPI_Barrier(m_mpiCommActive);
  if(!m_iRankActive) {
    std::cout<<" Just At End of DA...."<<std::endl<<std::endl;
  }
  //std::cout<<m_iRankActive<<" "<<m_uiPreGhostNodeSize<<" "<<m_uiPreGhostBoundaryNodeSize<<" "<<m_uiPostGhostNodeSize<<" "<<m_uiNodeSize<<" "<<m_uiBoundaryNodeSize<<std::endl<<std::endl;
  MPI_Barrier(m_mpiCommActive);
#endif

}//end function

#undef RESET_DA_BLOCK

}//end namespace



