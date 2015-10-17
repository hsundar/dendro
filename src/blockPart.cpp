
/**
  @file BlockPart.C
  @brief A Specialized Octree Partitioning Algorithm
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  @author Hari Sundar, hsundar@gmail.com
  */

#include "TreeNode.h"
#include "parUtils.h"
#include <omp.h>
#include <ompUtils.h>
#include "testUtils.h"
#ifdef __DEBUG__
#ifndef __DEBUG_OCT__
#define __DEBUG_OCT__
#endif
#endif

#ifdef __DEBUG_OCT__
#ifndef __MEASURE_BPART_COMM__
#define __MEASURE_BPART_COMM__
#endif
#endif

namespace ot {

  //Assumes that nodes are sorted.
  //No processor must call this with an empty input.
  int blockPartStage1_p2o(std::vector<TreeNode> &nodes, std::vector<TreeNode> &blocks,
      unsigned int dim, unsigned int maxDepth,  MPI_Comm comm) {
#ifdef __PROF_WITH_BARRIER__
    MPI_Barrier(comm);
#endif
    PROF_BLKPART1_BEGIN
      int npes, rank;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    par::partitionW<ot::TreeNode>(nodes, NULL, comm);

    //assert(nodes.size() > (1 << dim) ); 

    // 1. First compute the localCoarse octree 
    std::vector<TreeNode> localCoarse;

    //include both the min and max elements as well  
    //The output will be sorted, unique and linear
    appendCompleteRegion(nodes[0], nodes[nodes.size()-1], localCoarse, true, true);
    
    // std::cout << rank << ": in BlkPart_1 after appendCompleteRegion" << std::endl;
    
    // 2. Get local Blocks. These will be input to completeOctree that will
    // produce globalCoarse.
    std::vector<TreeNode> localBlocks;

    //Select all the coarsest blocks and insert them into localBlocks
    unsigned int minLevel = maxDepth;
    for(int i = 0; i < localCoarse.size(); i++) {
      if(localCoarse[i].getLevel() < minLevel) {
        minLevel = localCoarse[i].getLevel();
      }
    }

    for(int i = 0; i < localCoarse.size(); i++) {
      if(localCoarse[i].getLevel() == minLevel) {
        localBlocks.push_back(localCoarse[i]);
      }
    }

    localCoarse.clear();

    for (unsigned int i = 0; i < localBlocks.size(); i++) {
      localBlocks[i].setWeight(1);
    }

    // 3. Call nodes2Oct on these cells to generate the 
    //    globalCoarse octree ...
    //localBlocks will be sorted and linear 
    //There is a pathological case which prevents us from asserting that
    //localBlocks will be globally unique. For example, if the last element in the input on processor i
    //is the same as the first element on processor i+1 and if they are both
    //selected in localBlocks.

    completeOctree(localBlocks, blocks, dim, maxDepth, false, true, false, comm);

    localBlocks.clear();

    PROF_BLKPART1_END
  } // end blockPart

  int blockPartStage2_p2o(std::vector<TreeNode> &nodes, std::vector<TreeNode> &globalCoarse,
      std::vector<ot::TreeNode>& minsAllBlocks, unsigned int dim, unsigned int maxDepth,
      MPI_Comm comm) {

#ifdef __PROF_WITH_BARRIER__
    MPI_Barrier(comm);
#endif
    PROF_BLKPART2_BEGIN

      ot::TreeNode rootNode (dim,maxDepth);

    int npes, rank;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);

    int *sendCnt = new int[npes];
    int *recvCnt = new int[npes];
    int *sendOffsets = new int[npes];
    int *recvOffsets = new int[npes];

    //1. Now compute the wts of these cells ...
    //    a. Get the min and max nodes at each processor.
    std::vector<TreeNode> _mins_maxs(2*npes);

    // communicate ...
    ot::TreeNode sendMinMax[2];

    if (!nodes.empty()) {
      sendMinMax[0] =  nodes[0];
      sendMinMax[1] =  nodes[nodes.size()-1];
    } else {
      sendMinMax[0] = rootNode;
      sendMinMax[1] = rootNode;
    }

    par::Mpi_Allgather<ot::TreeNode>(sendMinMax, &(*_mins_maxs.begin()), 2, comm);

    std::vector<std::vector<TreeNode> > sendNodes(npes);
    std::vector<std::vector<unsigned int> > keymap(npes);

    for (int i=0; i<npes; i++) {
      sendCnt[i] = 0;
    }

    //    b. Now compute which cells go to which cells ...
    //       logic is that if the coarse cell is between the min and max at a
    //       processor or if it is an ancestor of min, then it is sent to that
    //       processor.
    //Naive Logic:
    for (unsigned int i=0; i<globalCoarse.size(); i++) {
      for (int p=0;p<npes;p++) {
#ifdef __DEBUG_OCT__
        assert(areComparable(globalCoarse[i], _mins_maxs[2*p]));
#endif
        if ( (globalCoarse[i].isAncestor(_mins_maxs[2*p])) ||
            ( (globalCoarse[i] >= _mins_maxs[2*p]) && (globalCoarse[i] <=_mins_maxs[(2*p)+1]) ) ) {
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
    par::Mpi_Alltoall<int>( sendCnt, recvCnt, 1, comm);


    //    b. Concatenate all nodes into one single Carray ...
    unsigned int totalSend = 0;
    unsigned int totalRecv = 0;
#ifdef __MEASURE_BPART_COMM__
    unsigned int numProcsSendI = 0;
    unsigned int numProcsRecvI = 0;
#endif
    for (unsigned int i = 0; i < npes; i++) {
      totalSend+= sendCnt[i];
      totalRecv+= recvCnt[i];
#ifdef __MEASURE_BPART_COMM__
      if(sendCnt[i]) {
        numProcsSendI++;
      }
      if(recvCnt[i]) {
        numProcsRecvI++;
      }
#endif
    }

    // create the send and recv buffers ...
    std::vector<ot::TreeNode> sendK (totalSend);
    std::vector<ot::TreeNode> recvK (totalRecv);

    // Now create sendK
    sendOffsets[0] = 0;
    recvOffsets[0] = 0;

    // compute offsets ...
    for (int i=1; i<npes; i++) {
      sendOffsets[i] = sendOffsets[i-1] + sendCnt[i-1];
      recvOffsets[i] = recvOffsets[i-1] + recvCnt[i-1];
    }

    for (int i=0; i<npes; i++) {
      for (unsigned int j=0; j<sendCnt[i]; j++) {
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

    par::Mpi_Alltoallv_sparse<ot::TreeNode>(sendKptr, sendCnt, sendOffsets,      
        recvKptr, recvCnt, recvOffsets, comm);

    sendK.clear();

#ifdef __MEASURE_BPART_COMM__
    MPI_Barrier(comm);
    unsigned int* allTotalSendI = new unsigned int[npes];
    unsigned int* allTotalRecvI = new unsigned int[npes];
    unsigned int* allNumProcsSendI = new unsigned int[npes];
    unsigned int* allNumProcsRecvI = new unsigned int[npes]; 
    par::Mpi_Gather<unsigned int>(&totalSend, allTotalSendI, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&totalRecv, allTotalRecvI, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&numProcsSendI, allNumProcsSendI, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&numProcsRecvI, allNumProcsRecvI, 1, 0, comm);
    if(!rank) {
      for(int i = 0; i < npes; i++) {
        std::cout<<" allTotalSendI["<<i<<"] in Bpart: "<<allTotalSendI[i]<<std::endl;
        std::cout<<" allTotalRecvI["<<i<<"] in Bpart: "<<allTotalRecvI[i]<<std::endl;
        std::cout<<" allNumProcsSendI["<<i<<"] in Bpart: "<<allNumProcsSendI[i]<<std::endl;
        std::cout<<" allNumProcsRecvI["<<i<<"] in Bpart: "<<allNumProcsRecvI[i]<<std::endl;
      }
    }
    delete [] allTotalSendI;
    delete [] allTotalRecvI;
    delete [] allNumProcsSendI;
    delete [] allNumProcsRecvI;
    MPI_Barrier(comm);
#endif

    //4. Now compute the wts of the locally received nodes ...
    //    a. loop through nodes and update the wts of the local chunks ...
    unsigned int *wts = NULL;
    if(totalRecv) {
      wts = new unsigned int [totalRecv];
    }
/*
    #pragma omp parallel for
    for (unsigned int i=0; i<totalRecv; i++) {
      wts[i] = 0;
    }

    //decendants and chunks are both sorted at this point.
    unsigned int nextPt = 0;
    unsigned int nextNode = 0;
    //Every element in nodes is inside some element in recvK.
    while (nextPt < nodes.size()) {
      //The first pt. lies in some block.
#ifdef __DEBUG_OCT__
      assert(areComparable(recvK[nextNode], nodes[nextPt]));
#endif
      if ((recvK[nextNode].isAncestor(nodes[nextPt])) ||
          (recvK[nextNode] == nodes[nextPt])) {
        wts[nextNode]++;
        nextPt++;
      } else {
        nextNode++;

        //If this fails then either recvK and nodes are not sorted or
        //Some pt in nodes is not in any of recvK
        assert(nextNode != totalRecv);
      }//end if-else
    }//end while
*/
    #pragma omp parallel for
    for(int i=0;i<totalRecv;i++){
      TreeNode* a=std::lower_bound(&nodes[0],&nodes[nodes.size()],recvK[i]);
      TreeNode* b=(i<totalRecv-1?std::lower_bound(&nodes[0],&nodes[nodes.size()],recvK[i+1]):&nodes[nodes.size()]);
//      int wts1=b-a;
//      wts1[i]=b-a;
//      if(wts1[i] > 0) assert(recvK[i].isAncestor(b[-1]) || recvK[i]==b[-1]);
//      assert(wts1[i]==wts[i]);
      wts[i]=b-a;
    }

    recvK.clear();

    //5. Now communicate the wts back to the procs ...
    unsigned int *recvWts = NULL;
    if(totalSend) {
      recvWts = new unsigned int[totalSend];
    }

    par::Mpi_Alltoallv_sparse<unsigned int>( wts, recvCnt, recvOffsets, 
        recvWts, sendCnt, sendOffsets, comm);

    //6. Now map them back to the blocks ...
    for (int i=0;i<globalCoarse.size();i++) {
      globalCoarse[i].setWeight(0);
    }

    for (int i=0; i<npes; i++) {
      for (int j=0; j<sendCnt[i]; j++) {
        globalCoarse[keymap[i][j]].addWeight(recvWts[sendOffsets[i] + j]);
      }
    }

    for (unsigned int i=0; i<npes; i++) {
      keymap[i].clear();
      sendNodes[i].clear();
    }
    sendNodes.clear();
    keymap.clear();

    if(recvWts) {
      delete [] recvWts;
      recvWts = NULL;
    }

    if(wts) {
      delete [] wts;
      wts = NULL;
    }

    par::partitionW<ot::TreeNode>(globalCoarse,getNodeWeight,comm);

    //Reset weights
    for (unsigned int i=0;i<globalCoarse.size(); i++) {
      globalCoarse[i].setWeight(1);
    }


    // Now communicate the nodes ...

    //7. Determine locally which keys to send to which proc ...
    //Compute Dist on globalCoarse....

    std::vector<TreeNode> vtkDist(npes);

    ot::TreeNode *sendMin = NULL;
    if (!globalCoarse.empty()) {
      sendMin = (ot::TreeNode *)&(*(globalCoarse.begin()));
    } else {
      sendMin = &(rootNode);
    }

    par::Mpi_Allgather<ot::TreeNode>(sendMin, &(* vtkDist.begin()), 1, comm);

    minsAllBlocks.clear();
    for(int j = 0; j < npes; j++) {
      if(vtkDist[j] != rootNode) {
        minsAllBlocks.push_back(vtkDist[j]);
      }
    }//end for j

    for (unsigned int j = 1; j < npes ; j++) {
      if (vtkDist[j] == rootNode) {
        vtkDist[j] = vtkDist[j-1];
      }
    }//end for j


    int max_non_root=npes;
    // correct dist ...
    if (npes>1) {
      if (vtkDist[npes-1] == vtkDist[npes-2]) {
        vtkDist[npes-1] = rootNode;
	max_non_root=npes-1;
      }//end if

      for (unsigned int i=npes-2; i>0; i--) {
        if (vtkDist[i] == vtkDist[i-1]) {
          vtkDist[i] = vtkDist[i+1];
          max_non_root=(vtkDist[i]==rootNode?i:max_non_root);
        }//end if
      }//end for
    }//end if npes > 1

    unsigned int *part = NULL;
    if(!nodes.empty()) {
      part = new unsigned int[nodes.size()];
    }

    if (npes > 1) {
/*      unsigned int pCnt=0;
      for (unsigned int i=0; i< nodes.size(); i++) {
        if ( (nodes[i] >= vtkDist[pCnt]) && ( (pCnt == (npes-1)) 
              || ( nodes[i] < vtkDist[pCnt+1] ) || (vtkDist[pCnt+1] == rootNode) ) ) {
          part[i] = pCnt;
        } else {
          while ( (pCnt < (npes -1)) && (nodes[i] >= vtkDist[pCnt+1])
              && (vtkDist[pCnt+1] != rootNode)  ) {
            pCnt++;
          }//end while
          part[i] = pCnt;
        }//end if-else
      }//end for i  */

      int omp_p=omp_get_max_threads();
      #pragma omp parallel for
      for(int j=0;j<omp_p;j++){
	int a=(nodes.size()*j)/omp_p;
	int b=(nodes.size()*(j+1))/omp_p;

        int pCnt=std::lower_bound(&vtkDist[0],&vtkDist[max_non_root],nodes[a])-&vtkDist[0]-1;
	pCnt=(pCnt>=0?pCnt:0);

        for (unsigned int i=a; i< b; i++) {
          if ( (nodes[i] >= vtkDist[pCnt]) && ( (pCnt == (npes-1)) 
                || ( nodes[i] < vtkDist[pCnt+1] ) || (vtkDist[pCnt+1] == rootNode) ) ) {
            part[i] = pCnt;
          } else {
            while ( (pCnt < (npes -1)) && (nodes[i] >= vtkDist[pCnt+1]) && (vtkDist[pCnt+1] != rootNode)  )
              pCnt++;
            part[i] = pCnt;
          }//end if-else
        }//end for i
      }
    }//end if np>1

    vtkDist.clear();
    //_________________________________________________________________________
    // Now the partitions should be contiguous since the two lists are globally
    // sorted ... and it's simply a shift between the two.
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //compute the total number of nodes being sent to each proc ...
    #pragma omp parallel for
    for (int i=0; i<npes; i++) {
      sendCnt[i]=0;
      recvCnt[i]=0;
    }

    if (npes > 1) {
/*      for (unsigned int i=0; i<nodes.size(); i++) {
        sendCnt[part[i]]++;
      }//end for i  */
      #pragma omp parallel for
      for(int i=part[0];i<=part[nodes.size()-1];i++){
	sendCnt[i]=std::lower_bound(&part[0],&part[nodes.size()],i+1) - 
	  std::lower_bound(&part[0],&part[nodes.size()],i);
      }

    } else {
      sendCnt[0] += (nodes.size());
    }//end if-else

    if(part) {
      delete [] part;
      part = NULL;
    }

    // communicate with other procs how many you shall be sending and get how
    // many to recieve from whom.
    par::Mpi_Alltoall<int>( sendCnt, recvCnt, 1, comm);

    sendOffsets[0] = 0;
    recvOffsets[0] = 0;
    // compute offsets ...
//    for (int i=1; i<npes; i++) {
//      sendOffsets[i] = sendOffsets[i-1] + sendCnt[i-1];
//      recvOffsets[i] = recvOffsets[i-1] + recvCnt[i-1];
//    }//end for i
    omp_par::scan(&sendCnt[0],&sendOffsets[0],npes);
    omp_par::scan(&recvCnt[0],&recvOffsets[0],npes);
    totalRecv=recvOffsets[npes-1]+recvCnt[npes-1];

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

#ifdef __MEASURE_BPART_COMM__
    unsigned int totalSendSize = nodes.size();
    unsigned int totalRecvSize = newNodes.size();
    unsigned int numProcsSendF = 0;
    unsigned int numProcsRecvF = 0;
    for(int i = 0; i < npes; i++) {
      if(sendCnt[i]) {
        numProcsSendF++;
      }
      if(recvCnt[i]) {
        numProcsRecvF++;
      }
    }
#endif

    par::Mpi_Alltoallv_sparse<ot::TreeNode>( nodesPtr, sendCnt, sendOffsets,      
        newNodesPtr, recvCnt, recvOffsets, comm);

    // reset the pointer ...
    nodes = newNodes;
    newNodes.clear();

    // clean up ...
    delete [] sendCnt;
    sendCnt = NULL;

    delete [] recvCnt;
    recvCnt = NULL;

    delete [] sendOffsets;
    sendOffsets = NULL;

    delete [] recvOffsets;
    recvOffsets = NULL;

#ifdef __MEASURE_BPART_COMM__
    MPI_Barrier(comm);
    unsigned int* allTotalSendF = new unsigned int[npes];
    unsigned int* allTotalRecvF = new unsigned int[npes];
    unsigned int* allNumProcsSendF = new unsigned int[npes];
    unsigned int* allNumProcsRecvF = new unsigned int[npes]; 
    par::Mpi_Gather<unsigned int>(&totalSendSize, allTotalSendF, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&totalRecvSize, allTotalRecvF, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&numProcsSendF, allNumProcsSendF, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&numProcsRecvF, allNumProcsRecvF, 1, 0, comm);
    if(!rank) {
      for(int i = 0; i < npes; i++) {
        std::cout<<" allTotalSendF["<<i<<"] in Bpart: "<<allTotalSendF[i]<<std::endl;
        std::cout<<" allTotalRecvF["<<i<<"] in Bpart: "<<allTotalRecvF[i]<<std::endl;
        std::cout<<" allNumProcsSendF["<<i<<"] in Bpart: "<<allNumProcsSendF[i]<<std::endl;
        std::cout<<" allNumProcsRecvF["<<i<<"] in Bpart: "<<allNumProcsRecvF[i]<<std::endl;
      }
    }
    delete [] allTotalSendF;
    delete [] allTotalRecvF;
    delete [] allNumProcsSendF;
    delete [] allNumProcsRecvF;
    MPI_Barrier(comm);
#endif

    PROF_BLKPART2_END
  } // end blockPart


  //Assumes that nodes are globally sorted, linear, complete, unique.
  //No processor must call this with an empty input.
  int blockPartStage1(std::vector<TreeNode> &nodes, std::vector<TreeNode> &blocks,
      unsigned int dim, unsigned int maxDepth,  MPI_Comm comm) {
#ifdef __PROF_WITH_BARRIER__
    MPI_Barrier(comm);
#endif
    PROF_BLKPART1_BEGIN
      int npes, rank;
    const double thFac = 0.5;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);

    par::partitionW<ot::TreeNode>(nodes, NULL,comm);

    //std::cout << rank << ": " << __func__ << ":Block Part Node Size:" <<nodes.size()<< std::endl;
    assert(nodes.size() > (1 << dim) ); 

#ifdef __BLOCK_PART_EQUALS_MORTON_PART__

    blocks.clear();

    //include both the min and max elements as well  
    //The output will be sorted, unique and linear
    appendCompleteRegion(nodes[0], nodes[nodes.size()-1], blocks, true, true);

#else

    // 1. First compute the localCoarse octree 
    std::vector<TreeNode> localCoarse;
    //include both the min and max elements as well  
    //The output will be sorted, unique and linear

    appendCompleteRegion(nodes[0], nodes[nodes.size()-1], localCoarse, true, true);
    //treeNodesTovtk(localCoarse,rank,"localCoarse");

    // 2. Get local Blocks. These will be input to completeOctree that will
    // produce globalCoarse.

    std::vector<TreeNode> localBlocks;

    //Logic: Use bPartComparator to assign priorities and pick the ones
    //with high priority. The number to be picked is determined by a percentage
    //of the total local wt.

    for (unsigned int i = 0; i < localCoarse.size(); i++) {
      localCoarse[i].setWeight(0);
    }

    //decendants and chunks are both sorted at this point.
    unsigned int nextPt = 0;
    unsigned int nextNode = 0;
    //Every element in nodes is inside some element in localCoarse.
    //@har: Why we need an addiational sort for the nodes ?
    std::vector<ot::TreeNode> sorted_nodes;

    while (nextPt < nodes.size()) {
      //The first pt. lies in some block.
#ifdef __DEBUG_OCT__
      assert(areComparable(localCoarse[nextNode], nodes[nextPt]));
#endif
      if ((localCoarse[nextNode].isAncestor(nodes[nextPt])) ||
          (localCoarse[nextNode] == nodes[nextPt])) {
        localCoarse[nextNode].addWeight(1);
        nextPt++;
      } else {
        nextNode++;
        if (nextNode == localCoarse.size()) {
          //If this fails then either the lists are not sorted
          //or there is some node which is not inside any block 
          assert(false);
        }
      }//end if-else
    }//end while

    sort(localCoarse.begin(), localCoarse.end(), ot::bPartComparator);

    long localWt = 0;
    unsigned int cnt = 0;
    while ( ( cnt < localCoarse.size() ) &&
        (localWt <= ( (long)( thFac*( (double)(nodes.size()) ) ) ) ) ) {
      localBlocks.push_back(localCoarse[cnt]);
      localWt += (long)(localCoarse[cnt].getWeight());
      cnt++;
    }//end while

    localCoarse.clear();

    for (unsigned int i = 0; i < localBlocks.size(); i++) {
      // std::cout << rank << ": block[" << i << "] = " << localBlocks[i].getWeight() << std::endl;
      localBlocks[i].setWeight(1);
    }

    //Sorting is necessary here since bPartcomparator is different from <.
    sort(localBlocks.begin(), localBlocks.end());

    // 3. Call nodes2Oct on these cells to generate the 
    //    globalCoarse octree ...
    //localBlocks will be sorted, linear and unique at this point 
    //localBlocks will not be empty on any processor

    completeOctree(localBlocks, blocks, dim, maxDepth, true, true, true, comm);

    localBlocks.clear();

#endif

    PROF_BLKPART1_END
  } // end blockPart

  int blockPartStage2(std::vector<TreeNode> &nodes, std::vector<TreeNode> &globalCoarse,
      std::vector<ot::TreeNode>& minsAllBlocks, unsigned int dim, unsigned int maxDepth,
      MPI_Comm comm) {
#ifdef __PROF_WITH_BARRIER__
    MPI_Barrier(comm);
#endif
    PROF_BLKPART2_BEGIN

      ot::TreeNode rootNode (dim,maxDepth);

    int npes, rank;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);

#ifdef __BLOCK_PART_EQUALS_MORTON_PART__
    //Just set minsAllBlocks. nodes and globalCoarse are already aligned.
    //We will not do weighted partitioning here.

    std::vector<TreeNode> vtkDist(npes);

    ot::TreeNode *sendMin = NULL;
    if (!globalCoarse.empty()) {
      sendMin = (ot::TreeNode *)&(*(globalCoarse.begin()));
    } else {
      sendMin = &(rootNode);
    }

    par::Mpi_Allgather<ot::TreeNode>(sendMin, &(* vtkDist.begin()), 1, comm);

    minsAllBlocks.clear();
    for(int j = 0; j < npes; j++) {
      if(vtkDist[j] != rootNode) {
        minsAllBlocks.push_back(vtkDist[j]);
      }
    }//end for j

#else

    int *sendCnt = new int[npes];
    int *recvCnt = new int[npes];
    int *sendOffsets = new int[npes];
    int *recvOffsets = new int[npes];

    //1. Now compute the wts of these cells ...
    //    a. Get the min and max nodes at each processor.
    std::vector<TreeNode> _mins_maxs(2*npes);

    // communicate ...
    ot::TreeNode sendMinMax[2];

    assert(par::test::isSorted(nodes, comm));
    assert(par::test::isSorted(globalCoarse, comm));

    if (!nodes.empty()) {
      sendMinMax[0] =  nodes[0];
      sendMinMax[1] =  nodes[nodes.size()-1];
    } else {
      sendMinMax[0] = rootNode;
      sendMinMax[1] = rootNode;
    }

    //std::cout << rank << ": min= " << sendMinMax[0] << ", max= " << sendMinMax[1] << std::endl;
    par::Mpi_Allgather<ot::TreeNode>(sendMinMax, &(*_mins_maxs.begin()), 2, comm);

    std::vector<std::vector<TreeNode> > sendNodes(npes);
    std::vector<std::vector<unsigned int> > keymap(npes);





    for (int i=0; i<npes; i++) {
      sendCnt[i] = 0;
    }

    //    b. Now compute which cells go to which cells ...
    //       logic is that if the coarse cell is between the min and max at a
    //       processor or if it is an ancestor of min, then it is sent to that
    //       processor.
    //Naive Logic:
    for (unsigned int i=0; i<globalCoarse.size(); i++) {
      for (int p=0;p<npes;p++) {
#ifdef __DEBUG_OCT__
        assert(areComparable(globalCoarse[i], _mins_maxs[2*p]));
#endif
         //if ( (globalCoarse[i].isAncestor(_mins_maxs[2*p])) || ( (globalCoarse[i] >= _mins_maxs[2*p]) && (globalCoarse[i] <=_mins_maxs[(2*p)+1]) ) ) {
        if ( (globalCoarse[i].isAncestor(_mins_maxs[2*p])) || ( (globalCoarse[i] >= _mins_maxs[2*p]) && (globalCoarse[i] <=_mins_maxs[(2*p)+1]) ) || (globalCoarse[i].isAncestor(_mins_maxs[2*p+1])) ) {
           sendNodes[p].push_back(globalCoarse[i]);
          // save keymap so that we can assign weights back to globalCoarse.
           // std::cout<<YLW<<"Global Coarse added:"<<globalCoarse[i]<<std::endl;
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
    par::Mpi_Alltoall<int>( sendCnt, recvCnt, 1, comm);

    //    b. Concatenate all nodes into one single Carray ...
    unsigned int totalSend = 0;
    unsigned int totalRecv = 0;
#ifdef __MEASURE_BPART_COMM__
    unsigned int numProcsSendI = 0;
    unsigned int numProcsRecvI = 0;
#endif
    for (unsigned int i = 0; i < npes; i++) {
      totalSend+= sendCnt[i];
      totalRecv+= recvCnt[i];
#ifdef __MEASURE_BPART_COMM__
      if(sendCnt[i]) {
        numProcsSendI++;
      }
      if(recvCnt[i]) {
        numProcsRecvI++;
      }
#endif
    }

    // create the send and recv buffers ...
    std::vector<ot::TreeNode> sendK (totalSend);
    std::vector<ot::TreeNode> recvK (totalRecv);

    // Now create sendK
    sendOffsets[0] = 0;
    recvOffsets[0] = 0;

    // compute offsets ...
    for (int i=1; i<npes; i++) {
      sendOffsets[i] = sendOffsets[i-1] + sendCnt[i-1];
      recvOffsets[i] = recvOffsets[i-1] + recvCnt[i-1];
    }

    for (int i=0; i<npes; i++) {
      for (unsigned int j=0; j<sendCnt[i]; j++) {
        //std::cout << rank << " -> " <<  i << " send: " << sendNodes[i][j] << std::endl;
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

//    if (! (seq::test::isSorted(sendK)) ) {
//      for (auto x: sendK)
//        std::cout << rank << ": " << x << std::endl;
//    }

    par::Mpi_Alltoallv_sparse<ot::TreeNode>(sendKptr, sendCnt, sendOffsets,
        recvKptr, recvCnt, recvOffsets, comm);

    sendK.clear();

#ifdef __MEASURE_BPART_COMM__
    MPI_Barrier(comm);
    unsigned int* allTotalSendI = new unsigned int[npes];
    unsigned int* allTotalRecvI = new unsigned int[npes];
    unsigned int* allNumProcsSendI = new unsigned int[npes];
    unsigned int* allNumProcsRecvI = new unsigned int[npes]; 
    par::Mpi_Gather<unsigned int>(&totalSend, allTotalSendI, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&totalRecv, allTotalRecvI, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&numProcsSendI, allNumProcsSendI, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&numProcsRecvI, allNumProcsRecvI, 1, 0, comm);
    if(!rank) {
      for(int i = 0; i < npes; i++) {
        std::cout<<" allTotalSendI["<<i<<"] in Bpart: "<<allTotalSendI[i]<<std::endl;
        std::cout<<" allTotalRecvI["<<i<<"] in Bpart: "<<allTotalRecvI[i]<<std::endl;
        std::cout<<" allNumProcsSendI["<<i<<"] in Bpart: "<<allNumProcsSendI[i]<<std::endl;
        std::cout<<" allNumProcsRecvI["<<i<<"] in Bpart: "<<allNumProcsRecvI[i]<<std::endl;
      }
    }
    delete [] allTotalSendI;
    delete [] allTotalRecvI;
    delete [] allNumProcsSendI;
    delete [] allNumProcsRecvI;
    MPI_Barrier(comm);
#endif

    //4. Now compute the wts of the locally received nodes ...
    //    a. loop through nodes and update the wts of the local chunks ...
    unsigned int *wts = NULL;
    if(totalRecv) {
      wts = new unsigned int [totalRecv];
    }
    for (unsigned int i=0; i<totalRecv; i++) {
      wts[i] = 0;
    }

    //decendants and chunks are both sorted at this point.
    unsigned int nextPt = 0;
    unsigned int nextNode = 0;

    // assert(seq::test::isSorted(nodes));
    // assert(seq::test::isSorted(recvK));

    //Every element in nodes is inside some element in recvK.
      std::cout<<RED<<"Block Part 2 While loop 1 Started"<<NRM<<std::endl;
    while (nextPt < nodes.size()) {
      //The first pt. lies in some block.
#ifdef __DEBUG_OCT__
      assert(areComparable(recvK[nextNode], nodes[nextPt]));
#endif
      if ((recvK[nextNode].isAncestor(nodes[nextPt])) || (recvK[nextNode] == nodes[nextPt])) {
        wts[nextNode]++;
        nextPt++;
      } else {
        nextNode++;
        if (nextNode == totalRecv) {
          //If this fails then either recvK and nodes are not sorted or
          //Some pt in nodes is not in any of recvK
          assert(false);
        }
      }//end if-else
    }//end while
      std::cout<<RED<<"Block Part 2 While loop 1 Ended"<<NRM<<std::endl;
    recvK.clear();

    //5. Now communicate the wts back to the procs ...
    unsigned int *recvWts = NULL;
    if(totalSend) {
      recvWts = new unsigned int[totalSend];
    }

    par::Mpi_Alltoallv_sparse<unsigned int>( wts, recvCnt, recvOffsets, 
        recvWts, sendCnt, sendOffsets, comm);

    //6. Now map them back to the blocks ...
    for (int i=0;i<globalCoarse.size();i++) {
      globalCoarse[i].setWeight(0);
    }

    for (int i=0; i<npes; i++) {
      for (int j=0; j<sendCnt[i]; j++) {
        globalCoarse[keymap[i][j]].addWeight(recvWts[sendOffsets[i] + j]);
      }
    }

    for (unsigned int i=0; i<npes; i++) {
      keymap[i].clear();
      sendNodes[i].clear();
    }
    sendNodes.clear();
    keymap.clear();

    if(recvWts) {
      delete [] recvWts;
      recvWts = NULL;
    }

    if(wts) {
      delete [] wts;
      wts = NULL;
    }

    /*
    int q=0;
    for (auto x: globalCoarse) {
      std::cout << rank << ": C[" << q++ << "] " << x << ", wt = " << x.getWeight() << std::endl;
    }
    */
    std::cout << RED " Before Partition" NRM << std::endl;
    par::partitionW<ot::TreeNode>(globalCoarse,getNodeWeight,comm);
    std::cout << RED " After Partition" NRM << std::endl;
    /*
    std::cout << RED " After Partition" NRM << std::endl;
    q=0;
    for (auto x: globalCoarse) {
      std::cout << rank << ": C[" << q++ << "] " << x << ", wt = " << x.getWeight() << std::endl;
    }
    */
    //Reset weights
    for (unsigned int i=0;i<globalCoarse.size(); i++) {
      globalCoarse[i].setWeight(1);
    }

    // Now communicate the nodes ...

    //7. Determine locally which keys to send to which proc ...
    //Compute Dist on globalCoarse....

    std::vector<TreeNode> vtkDist(npes);

    // ot::TreeNode *sendMin = NULL;
    ot::TreeNode sendMin;
    if (!globalCoarse.empty()) {
      sendMin = globalCoarse.begin()->getCFD();
    } else {
      sendMin = rootNode;
    }

    par::Mpi_Allgather<ot::TreeNode>(&sendMin, &(* vtkDist.begin()), 1, comm);

    minsAllBlocks.clear();
    for(int j = 0; j < npes; j++) {
      if(vtkDist[j] != rootNode) {
        minsAllBlocks.push_back(vtkDist[j]);
      }
    }//end for j

    for (unsigned int j = 1; j < npes ; j++) {
      if (vtkDist[j] == rootNode) {
        vtkDist[j] = vtkDist[j-1];
      }
    }//end for j

    //! correct dist ...
    if (npes>1) {
      if (vtkDist[npes-1] == vtkDist[npes-2]) {
        vtkDist[npes-1] = rootNode;
      }//end if

      for (unsigned int i=npes-2; i>0; i--) {
        if (vtkDist[i] == vtkDist[i-1]) {
          vtkDist[i] = vtkDist[i+1];
        }//end if
      }//end for
    }//end if npes > 1

    unsigned int *part = NULL;
    if(!nodes.empty()) {
      part = new unsigned int[nodes.size()];
    }

    if (npes > 1) {
      unsigned int pCnt=0;
      for (unsigned int i=0; i< nodes.size(); i++) {

        if ( (nodes[i] >= vtkDist[pCnt]) && ( (pCnt == (npes-1))
              || ( nodes[i] < vtkDist[pCnt+1] ) || (vtkDist[pCnt+1] == rootNode) ) ) {
          part[i] = pCnt;
        } else {

          while ( (pCnt < (npes -1)) && (nodes[i] >= vtkDist[pCnt+1])
              && (vtkDist[pCnt+1] != rootNode)  ) {
            pCnt++;
//              std::cout<<YLW<<"pCnt:"<<pCnt<<"/"<<(npes-1)<<NRM<<std::endl;
//              std::cout<<YLW<<"nodes "<<(i+1)<<"\t"<<nodes[i+1]<<"\t vtkDist[pCnt+1]:"<<vtkDist[pCnt+1]<<NRM<<std::endl;
          }//end while

          part[i] = pCnt;
        }//end if-else
      }//end for i
    }//end if np>1

      std::cout<<BLU<<"Line 1017"<<NRM<<std::endl;
    vtkDist.clear();
    //_________________________________________________________________________
    // Now the partitions should be contiguous since the two lists are globally
    // sorted ... and it's simply a shift between the two.
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //compute the total number of nodes being sent to each proc ...
    for (int i=0; i<npes; i++) {
      sendCnt[i]=0;
      recvCnt[i]=0;
    }

    if (npes > 1) {
      for (unsigned int i=0; i<nodes.size(); i++) {
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

    par::Mpi_Alltoall<int>( sendCnt, recvCnt, 1, comm);

    totalRecv=0;
    for (unsigned int i=0; i<npes; i++) {
      totalRecv += recvCnt[i];
    }//end for i

    sendOffsets[0] = 0;
    recvOffsets[0] = 0;

    // compute offsets ...
    for (int i=1; i<npes; i++) {
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

#ifdef __MEASURE_BPART_COMM__
    unsigned int totalSendSize = nodes.size();
    unsigned int totalRecvSize = newNodes.size();
    unsigned int numProcsSendF = 0;
    unsigned int numProcsRecvF = 0;
    for(int i = 0; i < npes; i++) {
      if(sendCnt[i]) {
        numProcsSendF++;
      }
      if(recvCnt[i]) {
        numProcsRecvF++;
      }
    }
#endif

    par::Mpi_Alltoallv_sparse<ot::TreeNode>( nodesPtr, sendCnt, sendOffsets,      
        newNodesPtr, recvCnt, recvOffsets, comm);

    // reset the pointer ...
    nodes = newNodes;
    newNodes.clear();

    // clean up ...
    delete [] sendCnt;
    sendCnt = NULL;

    delete [] recvCnt;
    recvCnt = NULL;

    delete [] sendOffsets;
    sendOffsets = NULL;

    delete [] recvOffsets;
    recvOffsets = NULL;

#ifdef __MEASURE_BPART_COMM__
    MPI_Barrier(comm);
    unsigned int* allTotalSendF = new unsigned int[npes];
    unsigned int* allTotalRecvF = new unsigned int[npes];
    unsigned int* allNumProcsSendF = new unsigned int[npes];
    unsigned int* allNumProcsRecvF = new unsigned int[npes]; 
    par::Mpi_Gather<unsigned int>(&totalSendSize, allTotalSendF, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&totalRecvSize, allTotalRecvF, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&numProcsSendF, allNumProcsSendF, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&numProcsRecvF, allNumProcsRecvF, 1, 0, comm);
    if(!rank) {
      for(int i = 0; i < npes; i++) {
        std::cout<<" allTotalSendF["<<i<<"] in Bpart: "<<allTotalSendF[i]<<std::endl;
        std::cout<<" allTotalRecvF["<<i<<"] in Bpart: "<<allTotalRecvF[i]<<std::endl;
        std::cout<<" allNumProcsSendF["<<i<<"] in Bpart: "<<allNumProcsSendF[i]<<std::endl;
        std::cout<<" allNumProcsRecvF["<<i<<"] in Bpart: "<<allNumProcsRecvF[i]<<std::endl;
      }
    }
    delete [] allTotalSendF;
    delete [] allTotalRecvF;
    delete [] allNumProcsSendF;
    delete [] allNumProcsRecvF;
    MPI_Barrier(comm);
#endif

#endif

    PROF_BLKPART2_END
  } // end blockPart

}//end namespace


