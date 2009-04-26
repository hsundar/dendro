
/**
  @file oda.txx
  @author Hari Sundar, hsundar@gmail.com
  @author Rahul S. Sampath, rahul.sampath@gmail.com
 */

#include "TreeNode.h"
#include "indexHolder.h"
#include "Sort.h" /**< The non-standard sort algorithm. */
#include "dtypes.h"
#include "odaUtils.h"
#include "parUtils.h"
#include "cnumEtypes.h"
#include "dendro.h"

#ifdef __DEBUG__
#ifndef __DEBUG_DA__
#define __DEBUG_DA__
#endif
#endif

namespace ot {

  inline unsigned int DA::curr() {
#ifdef __DEBUG_DA__
    assert(m_bIamActive);
#endif
    return m_uiCurrent;
  }

  inline unsigned int DA::currWithInfo() {	
#ifdef __DEBUG_DA__
    assert(m_bIamActive);
#endif
    //Store current info...
    m_lcLoopInfo.currentOffset = m_ptCurrentOffset;
    m_lcLoopInfo.currentIndex = m_uiCurrent;
    m_lcLoopInfo.qCounter = m_uiQuotientCounter;
    m_lcLoopInfo.pgQcounter = m_uiPreGhostQuotientCnt;
    return m_uiCurrent;
  }

  inline unsigned char DA::getFlag(unsigned int i) {
#ifdef __DEBUG_DA__
    assert(m_bIamActive);
#endif
    return m_ucpOctLevels[i] ;
  }

  inline unsigned char DA::getLevel(unsigned int i) {
#ifdef __DEBUG_DA__
    assert(m_bIamActive);
#endif
    return(m_ucpOctLevels[i] & ot::TreeNode::MAX_LEVEL );
  }

  inline bool DA::isNode(unsigned int i) {
#ifdef __DEBUG_DA__
    assert(m_bIamActive);
#endif
    return( m_ucpOctLevels[i] & ot::TreeNode::NODE);
  }

  inline bool DA::isGhost(unsigned int i) {
#ifdef __DEBUG_DA__
    assert(m_bIamActive);
#endif
    return( ( i < m_uiElementBegin) || (i >= m_uiPostGhostBegin) );
  }

  inline bool DA::isHanging(unsigned int i) {
#ifdef __DEBUG_DA__
    assert(m_bIamActive);
#endif
    return( m_ucpLutMasksPtr[(i<<1) + 1] );
  } 

  inline bool DA::computedLocalToGlobalElems() {
    return m_bComputedLocalToGlobalElems;
  }

  inline bool DA::computedLocalToGlobal() {
    return m_bComputedLocalToGlobal;
  }

  inline unsigned int DA::getDimension() {
    return m_uiDimension;
  }

  inline unsigned int DA::getMaxDepth() {
    return m_uiMaxDepth;
  }

  inline std::vector<ot::TreeNode> DA::getBlocks() {
    return m_tnBlocks;
  }

  inline unsigned int DA::getNumBlocks() {
    return static_cast<unsigned int>(m_tnBlocks.size());
  }

  inline std::vector<ot::TreeNode> DA::getMinAllBlocks() {
    return m_tnMinAllBlocks;
  }

  inline bool DA::isLUTcompressed() {
    return m_bCompressLut;
  }

  inline unsigned int DA::getLocalBufferSize() { 
    return m_uiLocalBufferSize; 
  }

  inline unsigned int DA::getElementSize() {
    return m_uiElementSize;
  }

  inline unsigned int DA::getPreGhostElementSize() {
    return m_uiPreGhostElementSize;
  }

  inline unsigned int DA::getIndependentSize() {
    return m_uiIndependentElementSize;
  }

  inline unsigned int DA::getGhostedNodeSize() {
    return (m_uiNodeSize + m_uiBoundaryNodeSize +
        m_uiPreGhostNodeSize + m_uiPreGhostBoundaryNodeSize +
        m_uiPostGhostNodeSize);
  }

  inline unsigned int DA::getGhostedElementSize() {
    return m_uiElementSize + m_uiPreGhostElementSize;
  }

  inline unsigned int DA::getInternalNodeSize() {
    return m_uiNodeSize;
  }

  inline unsigned int DA::getNodeSize() {
    return (m_uiNodeSize + m_uiBoundaryNodeSize);
  }

  inline unsigned int DA::getBoundaryNodeSize() {
    return m_uiBoundaryNodeSize;
  }

  inline unsigned int DA::getInputSize() {
    return m_uiInputSize;
  }

  inline unsigned int DA::getIdxElementBegin() {
    return m_uiElementBegin;
  }

  inline unsigned int DA::getIdxElementEnd() {
    return m_uiElementEnd;
  }

  inline unsigned int DA::getIdxPostGhostBegin() {
    return m_uiPostGhostBegin;
  }

  inline bool DA::iAmActive() {
    return m_bIamActive;
  }

  inline MPI_Comm DA::getComm() {
    return m_mpiCommAll;
  }

  inline MPI_Comm DA::getCommActive() {
    return m_mpiCommActive;
  }

  inline int DA::getNpesAll() {
    return m_iNpesAll;
  }

  inline int DA::getNpesActive() {
    return m_iNpesActive;
  }

  inline int DA::getRankAll() {
    return m_iRankAll;
  }

  inline int DA::getRankActive() {
    return m_iRankActive;
  }

  inline Point DA::getCurrentOffset() {
#ifdef __DEBUG_DA__
    assert(m_bIamActive);
#endif
    return m_ptCurrentOffset;
  }

  inline Point DA::getGhostedOffset() {
#ifdef __DEBUG_DA__
    assert(m_bIamActive);
#endif
    return m_ptGhostedOffset;
  }

  inline Point DA::getOffset() {
#ifdef __DEBUG_DA__
    assert(m_bIamActive);
#endif
    return m_ptOffset;
  } 

  inline DendroIntL* DA::getLocalToGlobalMap() {
    return m_dilpLocalToGlobal;
  }

  inline DendroIntL* DA::getLocalToGlobalElemsMap() {
    return m_dilpLocalToGlobalElems;
  }


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // 
  // Implementation...

  //Next functions...
  template<>
    inline unsigned int DA::next<ot::DA_FLAGS::ALL>() {
#ifdef __DEBUG_DA__
      assert(m_bIamActive);
#endif
      if ( m_uiElementBegin && (m_uiCurrent < (m_uiElementBegin - 1)) ) {
        incrementPreGhostOffset();
      } else {
        incrementCurrentOffset();
      }
      m_uiCurrent++;

      while ( (m_uiCurrent < (m_uiPreGhostElementSize + m_uiElementSize)) &&
          ( m_ucpLutMasks[(2*m_uiCurrent) + 1]  == ot::DA_FLAGS::FOREIGN) ) {

        if(m_bCompressLut) {
          updateQuotientCounter();
        }
        if ( m_uiElementBegin && (m_uiCurrent < (m_uiElementBegin - 1)) ) {
          incrementPreGhostOffset();
        } else {
          incrementCurrentOffset();
        }
        m_uiCurrent++;
      }
      //std::cout << YLW"Index is not Foreign element " << m_uiCurrent << std::endl;
      return m_uiCurrent;
    }//end function

  template<>
    inline unsigned int DA::next<ot::DA_FLAGS::INDEPENDENT>() {
#ifdef __DEBUG_DA__
      assert(m_bIamActive);
#endif
      incrementCurrentOffset();
      m_uiCurrent++;
      while ( (m_uiCurrent < m_uiIndependentElementEnd) &&
          (m_ucpOctLevels[m_uiCurrent] & ot::DA_FLAGS::DEP_ELEM) ) {
        // std::cout << RED"Skipping Dependent Element"NRM <<std::endl;
        if(m_bCompressLut) {
          updateQuotientCounter();        
        }
        incrementCurrentOffset();
        m_uiCurrent++;
      }
      return m_uiCurrent;
    }//end function

  template<>
    inline unsigned int DA::next<ot::DA_FLAGS::DEPENDENT>() {
#ifdef __DEBUG_DA__
      assert(m_bIamActive);
#endif
      if ( m_uiElementBegin && (m_uiCurrent < (m_uiElementBegin - 1)) ) {
        incrementPreGhostOffset();
      } else {
        incrementCurrentOffset();
      }
      m_uiCurrent++;

      while ( (m_uiCurrent < (m_uiPreGhostElementSize + m_uiElementSize)) &&
          (  (!(m_ucpOctLevels[m_uiCurrent] & ot::DA_FLAGS::DEP_ELEM )) ||
             ( m_ucpLutMasks[(2*m_uiCurrent) + 1] == ot::DA_FLAGS::FOREIGN) ) ) {               
        if(m_bCompressLut) {
          updateQuotientCounter();  
        }
        if ( m_uiElementBegin && (m_uiCurrent < (m_uiElementBegin - 1)) ) {
          incrementPreGhostOffset();
        } else {          
          incrementCurrentOffset();
        }
        m_uiCurrent++;
      }
      return m_uiCurrent;
    }//end function

  template<>
    inline unsigned int DA::next<ot::DA_FLAGS::W_DEPENDENT>() {
#ifdef __DEBUG_DA__
      assert(m_bIamActive);
#endif
      incrementCurrentOffset();
      m_uiCurrent++;
      while ( (m_uiCurrent < (m_uiPreGhostElementSize + m_uiElementSize)) &&
          ( (!(m_ucpOctLevels[m_uiCurrent] & ot::DA_FLAGS::DEP_ELEM )) ||
            ( m_ucpLutMasks[(2*m_uiCurrent) + 1] == ot::DA_FLAGS::FOREIGN) ) ) { 

        if(m_bCompressLut) {
          updateQuotientCounter();  
        }
        incrementCurrentOffset();          
        m_uiCurrent++;
      } 
      return m_uiCurrent;
    }//end function

  template<>
    inline unsigned int DA::next<ot::DA_FLAGS::WRITABLE>() {
#ifdef __DEBUG_DA__
      assert(m_bIamActive);
#endif
      incrementCurrentOffset();
      m_uiCurrent++;
      return m_uiCurrent;
    }//end function


  //Init functions...
  template<>	
    inline void DA::init<ot::DA_FLAGS::ALL>() {
#ifdef __DEBUG_DA__
      assert(m_bIamActive);
#endif
      m_ptCurrentOffset = m_ptGhostedOffset;
      m_uiCurrent = 0;
      m_uiQuotientCounter = 0;
      m_uiPreGhostQuotientCnt= 0;
      if ( (m_uiCurrent < (m_uiPreGhostElementSize + m_uiElementSize))
          && ( m_ucpLutMasks[(2*m_uiCurrent) + 1]  == ot::DA_FLAGS::FOREIGN) ) {
        if(m_bCompressLut) {
          updateQuotientCounter();
        }
        next<ot::DA_FLAGS::ALL>();
      }
    }//end function

  template<>	
    inline void DA::init<ot::DA_FLAGS::INDEPENDENT>() {
#ifdef __DEBUG_DA__
      assert(m_bIamActive);
#endif
      m_ptCurrentOffset = m_ptIndependentOffset;
      m_uiCurrent = m_uiIndependentElementBegin;
      m_uiQuotientCounter = m_uiIndependentElementQuotient;
      m_uiPreGhostQuotientCnt = static_cast<unsigned int>(m_ptsPreGhostOffsets.size());      
    }//end function

  template<>	
    inline void DA::init<ot::DA_FLAGS::DEPENDENT>() {
#ifdef __DEBUG_DA__
      assert(m_bIamActive);
#endif
      m_ptCurrentOffset = m_ptGhostedOffset;
      m_uiCurrent = 0;
      m_uiQuotientCounter = 0;
      m_uiPreGhostQuotientCnt= 0;
      //Process for the case where m_idxCurrent is actually independent ...
      if (m_iNpesActive==1) {
        m_uiCurrent = m_uiPreGhostElementSize + m_uiElementSize;           
      }else {
        if ( (m_uiCurrent < (m_uiPreGhostElementSize + m_uiElementSize) ) &&
            ( ( m_ucpLutMasks[(2*m_uiCurrent) + 1]  == ot::DA_FLAGS::FOREIGN) ||
              (!(m_ucpOctLevels[m_uiCurrent] & ot::DA_FLAGS::DEP_ELEM)) ) ) {
          if(m_bCompressLut) {
            updateQuotientCounter();  
          }
          next<ot::DA_FLAGS::DEPENDENT>(); 
        }
      }
    }//end function

  template<>	
    inline void DA::init<ot::DA_FLAGS::WRITABLE>() {
#ifdef __DEBUG_DA__
      assert(m_bIamActive);
#endif
      m_ptCurrentOffset = m_ptOffset;
      m_uiCurrent = m_uiElementBegin;
      m_uiQuotientCounter = m_uiElementQuotient;
      m_uiPreGhostQuotientCnt = m_ptsPreGhostOffsets.size();
    }//end function

  template<>	
    inline void DA::init<ot::DA_FLAGS::FROM_STORED>() {
#ifdef __DEBUG_DA__
      assert(m_bIamActive);
#endif
      //Use Stored info...
      m_ptCurrentOffset = m_lcLoopInfo.currentOffset;
      m_uiCurrent = m_lcLoopInfo.currentIndex;
      m_uiQuotientCounter = m_lcLoopInfo.qCounter;
      m_uiPreGhostQuotientCnt = m_lcLoopInfo.pgQcounter;			
    }//end function

  template<>	
    inline void DA::init<ot::DA_FLAGS::W_DEPENDENT>() {
#ifdef __DEBUG_DA__
      assert(m_bIamActive);
#endif
      m_ptCurrentOffset = m_ptOffset;
      m_uiCurrent = m_uiElementBegin;
      m_uiQuotientCounter = m_uiElementQuotient;
      m_uiPreGhostQuotientCnt = m_ptsPreGhostOffsets.size();
      //Process for the case where m_idxCurrent is actually independent ...
      if (m_iNpesActive==1) {
        m_uiCurrent = m_uiPreGhostElementSize + m_uiElementSize;
      }else {
        if ( (m_uiCurrent < (m_uiPreGhostElementSize + m_uiElementSize) ) &&
            ( ( m_ucpLutMasks[(2*m_uiCurrent) + 1]  == ot::DA_FLAGS::FOREIGN) ||
              (!(m_ucpOctLevels[m_uiCurrent] & ot::DA_FLAGS::DEP_ELEM)) ) ) {
          if(m_bCompressLut) {
            updateQuotientCounter();  
          }
          next<ot::DA_FLAGS::W_DEPENDENT>(); 
        }
      }
    }//end function

  //End functions...
  template<>
    inline unsigned int DA::end<ot::DA_FLAGS::ALL>() {
#ifdef __DEBUG_DA__
      assert(m_bIamActive);
#endif
      return (m_uiPreGhostElementSize + m_uiElementSize);
    }

  template<>
    inline unsigned int DA::end<ot::DA_FLAGS::DEPENDENT>() {
#ifdef __DEBUG_DA__
      assert(m_bIamActive);
#endif
      return (m_uiPreGhostElementSize + m_uiElementSize);
    }

  template<>
    inline unsigned int DA::end<ot::DA_FLAGS::WRITABLE>() {
#ifdef __DEBUG_DA__
      assert(m_bIamActive);
#endif
      return (m_uiPreGhostElementSize + m_uiElementSize);
    }

  template<>
    inline unsigned int DA::end<ot::DA_FLAGS::W_DEPENDENT>() {
#ifdef __DEBUG_DA__
      assert(m_bIamActive);
#endif
      return (m_uiPreGhostElementSize + m_uiElementSize);
    }

  template<>
    inline unsigned int DA::end<ot::DA_FLAGS::INDEPENDENT>() {
#ifdef __DEBUG_DA__
      assert(m_bIamActive);
#endif
      return m_uiIndependentElementEnd;
    }

  //Functions for communicating ghost nodes...

  template <typename T>
    int DA::ReadFromGhostsBegin ( T* arr, unsigned int dof) {
      PROF_READ_GHOST_NODES_BEGIN_BEGIN

        // first need to create contiguous list of boundaries ...
        T* sendK = NULL;
      if(m_uipScatterMap.size()) {
        sendK = new T[dof*m_uipScatterMap.size()];
      }
      for (unsigned int i = 0; i < m_uipScatterMap.size(); i++ ) {
        // add dof loop ...
        for (unsigned int j = 0; j < dof; j++) {
          sendK[(dof*i) + j] = arr[(dof*m_uipScatterMap[i]) + j];
        }
      }

      // create a new context ...
      updateContext ctx;
      ctx.buffer = arr;
      ctx.keys = sendK;

      // Post Recv ...
      for (unsigned int i = 0; i < m_uipRecvProcs.size(); i++) {
        MPI_Request *req = new MPI_Request();
        par::Mpi_Irecv<T>(arr + (dof*m_uipRecvOffsets[i]), (dof*m_uipRecvCounts[i]), 
            m_uipRecvProcs[i], m_uiCommTag, m_mpiCommActive, req );
        ctx.requests.push_back(req);
      }

      //*********** Send ****************//
      for (unsigned int i = 0; i < m_uipSendProcs.size(); i++) {
        MPI_Request *req = new MPI_Request();
        par::Mpi_Isend<T>( sendK + (dof*m_uipSendOffsets[i]), (dof*m_uipSendCounts[i]),
            m_uipSendProcs[i], m_uiCommTag, m_mpiCommActive, req );
        ctx.requests.push_back(req);
      }

      // Increment tag ....
      m_uiCommTag++;

      m_mpiContexts.push_back(ctx);

      PROF_READ_GHOST_NODES_BEGIN_END
    }

  template <typename T>
    int DA::ReadFromGhostsEnd(T* arr) {
      PROF_READ_GHOST_NODES_END_BEGIN

        // find the context ...
        unsigned int ctx;
      for ( ctx = 0; ctx < m_mpiContexts.size(); ctx++) {
        if ( m_mpiContexts[ctx].buffer == arr) {
          break;
        }
      }

      MPI_Status status;
      // need to wait for the commns to finish ...
      for (unsigned int i = 0; i < m_mpiContexts[ctx].requests.size(); i++) {
        MPI_Wait(m_mpiContexts[ctx].requests[i], &status);
        delete m_mpiContexts[ctx].requests[i];
      }

      // delete the sendkeys ...
      T *sendK = static_cast<T *>(m_mpiContexts[ctx].keys);

      if(sendK) {
        delete [] sendK;
        sendK = NULL;
      }

      // clear the Requests ...
      assert(ctx < m_mpiContexts.size());
      m_mpiContexts[ctx].requests.clear();

      // remove the context ...
      m_mpiContexts.erase(m_mpiContexts.begin() + ctx);

      PROF_READ_GHOST_NODES_END_END
    }


  template <typename T>
    int DA::WriteToGhostsBegin ( T* arr, unsigned int dof) {
      PROF_WRITE_GHOST_NODES_BEGIN_BEGIN

        // first need to create contiguous list of boundaries ...
        T* recvK = NULL;
      if(m_uipScatterMap.size()) {
        recvK = new T[dof*(m_uipScatterMap.size())];
      }

      // create a new context ...
      updateContext ctx;
      ctx.buffer = arr;
      ctx.keys = recvK;

      // Post Recv ...
      for (unsigned int i = 0; i < m_uipSendProcs.size(); i++) {
        MPI_Request *req = new MPI_Request();
        par::Mpi_Irecv<T>( recvK + (dof*m_uipSendOffsets[i]), (dof*m_uipSendCounts[i]), 
            m_uipSendProcs[i], m_uiCommTag, m_mpiCommActive, req );
        ctx.requests.push_back(req);
      }

      //The communication here is just the opposite of the communication in readFromGhosts...
      //*********** Send ****************//
      for (unsigned int i = 0; i < m_uipRecvProcs.size(); i++) {
        MPI_Request *req = new MPI_Request();
        par::Mpi_Isend<T>( arr + (dof*m_uipRecvOffsets[i]), (dof*m_uipRecvCounts[i]), 
            m_uipRecvProcs[i], m_uiCommTag, m_mpiCommActive, req );
        ctx.requests.push_back(req);
      }

      // Increment tag ....
      m_uiCommTag++;

      m_mpiContexts.push_back(ctx);

      PROF_WRITE_GHOST_NODES_BEGIN_END
    }

  template <typename T>
    int DA::WriteToGhostsEnd(T* arr, unsigned int dof) {
      PROF_WRITE_GHOST_NODES_END_BEGIN

        // find the context ...
        unsigned int ctx;
      for ( ctx = 0; ctx < m_mpiContexts.size(); ctx++) {
        if ( m_mpiContexts[ctx].buffer == arr) {
          break;
        }
      }

      MPI_Status status;
      // need to wait for the commns to finish ...
      for (unsigned int i = 0; i < m_mpiContexts[ctx].requests.size(); i++) {
        MPI_Wait(m_mpiContexts[ctx].requests[i], &status);
        delete m_mpiContexts[ctx].requests[i];
      }

      //Add ghost values to the local vector.
      T *recvK = static_cast<T *>(m_mpiContexts[ctx].keys);
      for (unsigned int i=0; i<m_uipScatterMap.size(); i++ ) {
        for (unsigned int j=0; j<dof; j++) {
          arr[(dof*m_uipScatterMap[i]) + j] += recvK[(dof*i) + j];
        }
      }

      // delete the keys ...
      if(recvK) {
        delete [] recvK;
        recvK = NULL;
      }

      // clear the Requests ...
      assert(ctx < m_mpiContexts.size());
      m_mpiContexts[ctx].requests.clear();

      // remove the context ...
      m_mpiContexts.erase(m_mpiContexts.begin() + ctx);

      PROF_WRITE_GHOST_NODES_END_END
    }

  //Functions for communicating ghost elements...

  template <typename T>
    int DA::ReadFromGhostElemsBegin ( T* arr, unsigned int dof) {
      PROF_READ_GHOST_ELEMS_BEGIN_BEGIN

        // first need to create contiguous list of boundaries ...
        T* sendK = NULL;

      if(m_uipElemScatterMap.size()) {
        sendK = new T[dof*m_uipElemScatterMap.size()];
      }

      for (unsigned int i = 0; i < m_uipElemScatterMap.size(); i++ ) {
        // add dof loop ...
        for (unsigned int j = 0; j < dof; j++) {
          sendK[dof*i + j] = arr[dof*m_uipElemScatterMap[i]+j];
        }
      }

      // create a new context ...
      updateContext ctx;
      ctx.buffer = arr;
      ctx.keys = sendK;

      // Post Recv ...
      for (unsigned int i = 0; i < m_uipElemRecvProcs.size(); i++) {
        MPI_Request *req = new MPI_Request();
        par::Mpi_Irecv<T>( arr + dof*m_uipElemRecvOffsets[i], dof*m_uipElemRecvCounts[i],
            m_uipElemRecvProcs[i], m_uiCommTag, m_mpiCommActive, req );
        ctx.requests.push_back(req);
      }

      //*********** Send ****************//
      for (unsigned int i = 0; i < m_uipElemSendProcs.size(); i++) {
        MPI_Request *req = new MPI_Request();
        par::Mpi_Issend<T>( sendK + dof*m_uipElemSendOffsets[i], dof*m_uipElemSendCounts[i],
            m_uipElemSendProcs[i], m_uiCommTag, m_mpiCommActive, req );
        ctx.requests.push_back(req);
      }

      // Increment tag ....
      m_uiCommTag++;

      m_mpiContexts.push_back(ctx);

      PROF_READ_GHOST_ELEMS_BEGIN_END
    }

  template <typename T>
    int DA::ReadFromGhostElemsEnd(T* arr) {
      PROF_READ_GHOST_ELEMS_END_BEGIN

        // find the context ...
        unsigned int ctx;
      for ( ctx = 0; ctx < m_mpiContexts.size(); ctx++) {
        if ( m_mpiContexts[ctx].buffer == arr) {
          break;
        }
      }

      MPI_Status status;
      // need to wait for the commns to finish ...
      for (unsigned int i = 0; i < m_mpiContexts[ctx].requests.size(); i++) {
        MPI_Wait(m_mpiContexts[ctx].requests[i], &status);
        delete m_mpiContexts[ctx].requests[i];
      }

      // delete the sendkeys ...
      T *sendK = static_cast<T *>(m_mpiContexts[ctx].keys);

      if(sendK) {
        delete [] sendK;
        sendK = NULL;
      }

      // clear the Requests ...
      assert(ctx < m_mpiContexts.size());
      m_mpiContexts[ctx].requests.clear();

      // remove the context ...
      m_mpiContexts.erase(m_mpiContexts.begin() + ctx);

      PROF_READ_GHOST_ELEMS_END_END
    }


  template <typename T>
    int DA::WriteToGhostElemsBegin ( T* arr, unsigned int dof) {

      PROF_WRITE_GHOST_ELEMS_BEGIN_BEGIN

        // first need to create contiguous list of boundaries ...
        T* recvK = NULL;
      if(m_uipElemScatterMap.size()) {
        recvK = new T[dof*m_uipElemScatterMap.size()];
      }

      // create a new context ...
      updateContext ctx;
      ctx.buffer = arr;
      ctx.keys = recvK;

      // and Recv ...
      for (unsigned int i = 0; i < m_uipElemSendProcs.size(); i++) {
        MPI_Request *req = new MPI_Request();
        par::Mpi_Irecv<T>( recvK + dof*m_uipElemSendOffsets[i], dof*m_uipElemSendCounts[i],
            m_uipElemSendProcs[i], m_uiCommTag, m_mpiCommActive, req );
        ctx.requests.push_back(req);
      }

      //The communication here is just the opposite of the communication in readFromGhosts...
      //*********** Send ****************//
      for (unsigned int i = 0; i < m_uipElemRecvProcs.size(); i++) {
        MPI_Request *req = new MPI_Request();
        par::Mpi_Issend<T>( arr + dof*m_uipElemRecvOffsets[i], dof*m_uipElemRecvCounts[i],
            m_uipElemRecvProcs[i], m_uiCommTag, m_mpiCommActive, req );
        ctx.requests.push_back(req);
      }

      // Increment tag ....
      m_uiCommTag++;

      m_mpiContexts.push_back(ctx);

      PROF_WRITE_GHOST_ELEMS_BEGIN_END
    }

  template <typename T>
    int DA::WriteToGhostElemsEnd(T* arr, unsigned int dof) {

      PROF_WRITE_GHOST_ELEMS_END_BEGIN

        // find the context ...
        unsigned int ctx;
      for ( ctx = 0; ctx < m_mpiContexts.size(); ctx++) {
        if ( m_mpiContexts[ctx].buffer == arr) {
          break;
        }
      }

      MPI_Status status;
      // need to wait for the commns to finish ...
      for (unsigned int i = 0; i < m_mpiContexts[ctx].requests.size(); i++) {
        MPI_Wait(m_mpiContexts[ctx].requests[i], &status);
        delete m_mpiContexts[ctx].requests[i];
      }

      //Add ghost values to the local vector.
      T *recvK = static_cast<T *>(m_mpiContexts[ctx].keys);
      for (unsigned int i = 0; i < m_uipElemScatterMap.size(); i++ ) {
        for (unsigned int j = 0; j < dof; j++) {
          arr[dof*m_uipElemScatterMap[i]+j] += recvK[dof*i + j];
        }
      }

      // delete the keys ...
      if(recvK) {
        delete [] recvK;
        recvK = NULL;
      }

      // clear the Requests ...
      assert(ctx < m_mpiContexts.size());
      m_mpiContexts[ctx].requests.clear();

      // remove the context ...
      m_mpiContexts.erase(m_mpiContexts.begin() + ctx);

      PROF_WRITE_GHOST_ELEMS_END_END
    }


  template <typename T>
    int DA::createVector(std::vector<T> &arr, bool isElemental,
        bool isGhosted, unsigned int dof) {
      // first determine the length of the vector ...
      unsigned int sz = 0;
      if(m_bIamActive) {
        if (isElemental) {
          sz = m_uiElementSize;
          if (isGhosted) {
            sz += m_uiPreGhostElementSize;
          }
        } else {
          sz = m_uiNodeSize + m_uiBoundaryNodeSize;
          if (isGhosted) {
            sz += (m_uiPreGhostNodeSize + m_uiPreGhostBoundaryNodeSize + m_uiPostGhostNodeSize);
          }
        }
        // now for dof ...
        sz *= dof;
      }//end if active

      // now create the vector
      arr.resize(sz);

      return 0;
    }

  template < typename T >
    int DA::vecGetBuffer(std::vector<T> &in, T* &out, bool isElemental,
        bool isGhosted, bool isReadOnly, unsigned int dof) {

      // Some error checks ... make sure the size of Vec in matches those implied
      // by the other params ...
      unsigned int sz = 0;
      if(m_bIamActive) {
        if (isElemental) {
          sz = m_uiElementSize;
          if (isGhosted) {
            sz += m_uiPreGhostElementSize;
          }
        } else {
          sz = m_uiNodeSize + m_uiBoundaryNodeSize;
          if (isGhosted) {
            sz += (m_uiPreGhostNodeSize + m_uiPreGhostBoundaryNodeSize + m_uiPostGhostNodeSize);
          }
        }
        // now for dof ...
        sz *= dof;
      }//end if active

      unsigned int vecSz = static_cast<unsigned int>(in.size());

      if ( sz != vecSz) {
        std::cerr  << "In function " << __func__ << " sizes are unequal, sz is  " << sz
          << " and vecSz is " << vecSz << std::endl;
        assert(false);
        return -1;
      };

      if(!m_bIamActive) {
        assert(m_uiLocalBufferSize == 0);
        assert(m_uiElementBegin == 0);
        assert(m_uiElementEnd == 0);
        assert(m_uiPostGhostBegin == 0);
      }

      // get the local Arrray,
      T *array = NULL;
      if(!in.empty()) {
        array = &(*(in.begin()));
      }

      if(isGhosted && isElemental) {
        //simply copy the pointer
        //This is the only case where the buffer will not be the size of the
        //fullLocalBufferSize. 
        out = array;
      }else {
        // First let us allocate for the buffer ... the local buffer will be of full
        // length.
        sz = dof*m_uiLocalBufferSize;
        //The default constructor of datatype T is responsible of initializing
        //out with the appropriate zero entries. 
        if(sz) {
          out = new T[sz];
        } else {
          out = NULL;
        }
      }

      unsigned int vecCnt=0;
      // Now we can populate the out buffer ... and that needs a loop through the
      // elements ...
      if (isGhosted) {
        if (isElemental) {
          //nothing to be done here.
        } else {
          //nodal and ghosted
          // now copy ...
          for (unsigned int i=0; i<m_uiLocalBufferSize; i++) {
            // skip the ones that are not nodes ...
            if ( ! (m_ucpOctLevels[i] & ot::TreeNode::NODE ) ) {
              continue;
            }
            for (unsigned int j=0; j<dof; j++) {
              out[dof*i+j] = array[dof*vecCnt + j];
            }
            vecCnt++;
          }
        }
      } else {
        if (isElemental) {
          //elemental and non-ghosted
          // is a simple copy ... 
          for (unsigned int i = m_uiElementBegin; i < m_uiElementEnd; i++) {
            for (unsigned int j = 0; j < dof; j++) {
              out[dof*i+j] = array[dof*vecCnt + j];
            }
            vecCnt++;
          }
        } else {
          //nodal and non-ghosted
          for (unsigned int i = m_uiElementBegin; i < m_uiElementEnd; i++) {
            if ( ! (m_ucpOctLevels[i] & ot::TreeNode::NODE ) ) {
              continue;
            }
            for (unsigned int j = 0; j < dof; j++) {
              out[dof*i+j] = array[dof*vecCnt + j];
            }
            vecCnt++;
          }
          for (unsigned int i = m_uiElementEnd; i < m_uiPostGhostBegin; i++) {
            // add the remaining boundary nodes ...
            if ( ! ( (m_ucpOctLevels[i] & ot::TreeNode::NODE ) &&
                  (m_ucpOctLevels[i] & ot::TreeNode::BOUNDARY ) ) ) {
              continue;
            }
            for (unsigned int j = 0; j < dof; j++) {
              out[dof*i+j] = array[dof*vecCnt + j];
            }
            vecCnt++;
          }
        }
      }
      return 0;
    }

  template < typename T >
    int DA::vecRestoreBuffer(std::vector<T> &in, T* out, bool isElemental,
        bool isGhosted, bool isReadOnly, unsigned int dof) {
      // Some error checks ... make sure the size of Vec in matches those implied
      // by the other params ...
      unsigned int sz = 0;
      if(m_bIamActive) {
        if (isElemental) {
          sz = m_uiElementSize;
          if (isGhosted) {
            sz += m_uiPreGhostElementSize;
          }
        } else {
          sz = m_uiNodeSize + m_uiBoundaryNodeSize;
          if (isGhosted) {
            sz += (m_uiPreGhostNodeSize + m_uiPreGhostBoundaryNodeSize + m_uiPostGhostNodeSize);
          }
        }
        // now for dof ...
        sz *= dof;
      }//end if active

      unsigned int vecSz = static_cast<unsigned int>(in.size());

      if ( sz != vecSz) {
        std::cerr  << "In function STL::" << __func__ << " sizes are unequal, sz is  " <<
          sz << " and vecSz is " << vecSz << std::endl;
        assert(false);
        return -1;
      };

      if(!m_bIamActive) {
        assert(m_uiLocalBufferSize == 0);
        assert(m_uiElementBegin == 0);
        assert(m_uiElementEnd == 0);
        assert(m_uiPostGhostBegin == 0);
      }

      unsigned int vecCnt=0;

      // if is readonly, then simply deallocate and return ...
      if ( isGhosted && isElemental ) {
        out = NULL;
      } else if ( isReadOnly ) {
        // no need to write back ... simply clean up and return
        if(out) {
          delete [] out;
          out = NULL;
        }
      } else {
        // need to write back ...
        // get the local Arrray,
        T *array = NULL;
        if(!in.empty()) {
          array = &(*in.begin());
        }

        if ( isElemental ) {
          //elemental and non-ghosted
          for (unsigned int i = m_uiElementBegin; i < m_uiElementEnd; i++) {
            for (unsigned int j = 0; j < dof; j++) {
              array[dof*vecCnt + j] = out[dof*i+j];
            }
            vecCnt++;
          }
        } else if ( isGhosted ) {
          // nodal and ghosted ...
          for (unsigned int i = 0; i < m_uiLocalBufferSize; i++) {
            // skip the ones that are not nodes ...
            if ( ! (m_ucpOctLevels[i] & ot::TreeNode::NODE ) ) {
              continue;
            }
            for (unsigned int j = 0; j < dof; j++) {
              array[dof*vecCnt + j] = out[dof*i+j];
            }
            vecCnt++;
          }
        } else {
          // nodal non ghosted ...
          for (unsigned int i = m_uiElementBegin; i < m_uiElementEnd; i++) {
            if ( ! (m_ucpOctLevels[i] & ot::TreeNode::NODE ) ) {
              continue;
            }
            for (unsigned int j = 0; j < dof; j++) {
              array[dof*vecCnt + j] = out[dof*i+j];
            }
            vecCnt++;
          }
          for (unsigned int i = m_uiElementEnd; i < m_uiPostGhostBegin; i++) {
            // add the remaining boundary nodes ...
            if ( ! ( (m_ucpOctLevels[i] & ot::TreeNode::NODE ) &&
                  (m_ucpOctLevels[i] & ot::TreeNode::BOUNDARY ) ) ) {
              continue;
            }
            for (unsigned int j = 0; j < dof; j++) {
              array[dof*vecCnt + j] = out[dof*i+j];
            }
            vecCnt++;
          }
        }

        if(out) {
          delete [] out;
          out = NULL;
        }  
      }
      return 0;
    }

  inline unsigned char DA::getChildNumber() {
#ifdef __DEBUG_DA__
    assert(m_bIamActive);
#endif

    if( m_bCompressLut && m_ucpLutMasksPtr[(m_uiCurrent << 1) + 1] ) {
      return ( m_ucpSortOrdersPtr[m_uiCurrent] );
    }else {
      unsigned int lev = m_ucpOctLevels[m_uiCurrent] & ot::TreeNode::MAX_LEVEL;
      unsigned int len = (unsigned int)(1u << ( m_uiMaxDepth - lev ) );
      unsigned int len_par = (unsigned int)(1u << ( m_uiMaxDepth  - lev +1 ) );

      unsigned int i = m_ptCurrentOffset.xint() % len_par;
      unsigned int j = m_ptCurrentOffset.yint() % len_par;
      unsigned int k = m_ptCurrentOffset.zint() % len_par;

      i /= len;
      j /= len;
      k /= len;

      return static_cast<unsigned char>((k << 2) + (j << 1) + i); //4*k+2*j+i
    }
  }

  inline int DA::getNodeIndices(unsigned int* nodes) {
#ifdef __DEBUG_DA__
    assert(m_bIamActive);
#endif
#ifdef __DEBUG_DA__    
    if ( m_ucpLutMasksPtr[2*m_uiCurrent+1] == ot::DA_FLAGS::FOREIGN ) {
      std::cerr << RED"Passing foreign mask for index "NRM << m_uiCurrent << std::endl;
      assert(false);
    }
#endif    

    unsigned int ii = (m_uiCurrent << 3); //*8

    if(!m_bCompressLut) {
      for(unsigned int j = 0; j < 8; j++) {
        nodes[j] = m_uiNlistPtr[ii + j];
      }
    }else {
      // m_ucpLutRemainders, m_uspLutQuotients.
      unsigned char _mask = m_ucpLutMasksPtr[(m_uiCurrent << 1)];

      // get the index into the 8 nodes ...
      unsigned int nn[8];

      // first let us get the offsets ..

      nn[0] = m_ucpLutRemaindersPtr[ii];
      if ( _mask & 1) {
        nn[0] += (static_cast<unsigned int>(m_uspLutQuotientsPtr[m_uiQuotientCounter++]) << 8);
      }
      nn[0] = m_uiCurrent - nn[0];

      for (unsigned int j = 1; j < 8; j++) {
        nn[j] = m_ucpLutRemaindersPtr[ii + j];
        if ( _mask & (1 << j ) ) {
          nn[j] += (static_cast<unsigned int>(m_uspLutQuotientsPtr[m_uiQuotientCounter++]) << 8); 
        }
        nn[j] +=  nn[j-1];
      }

      // Unsorting

      unsigned char hnMask = m_ucpLutMasksPtr[(m_uiCurrent << 1) + 1];

      if ( hnMask ) {

        unsigned int d = (m_ucpOctLevels[m_uiCurrent] & ot::TreeNode::MAX_LEVEL);
        unsigned int sz = ( 1u << (m_uiMaxDepth - d) );

        unsigned int x = m_ptCurrentOffset.xint();
        unsigned int y = m_ptCurrentOffset.yint();
        unsigned int z = m_ptCurrentOffset.zint();

        unsigned int xp = x + sz;
        unsigned int yp = y + sz;
        unsigned int zp = z + sz;

        unsigned int childNum = m_ucpSortOrdersPtr[m_uiCurrent];

        ot::TreeNode vertices[8];
        //initialize vertices...
        {
          vertices[0] = TreeNode(x,y,z,m_uiMaxDepth,m_uiDimension,m_uiMaxDepth);

          vertices[1] = TreeNode(xp,y,z,m_uiMaxDepth,m_uiDimension,m_uiMaxDepth);

          vertices[2] = TreeNode(x,yp,z,m_uiMaxDepth,m_uiDimension,m_uiMaxDepth);

          vertices[3] = TreeNode(xp,yp,z,m_uiMaxDepth,m_uiDimension,m_uiMaxDepth);

          vertices[4] = TreeNode(x,y,zp,m_uiMaxDepth,m_uiDimension,m_uiMaxDepth);

          vertices[5] = TreeNode(xp,y,zp,m_uiMaxDepth,m_uiDimension,m_uiMaxDepth);

          vertices[6] = TreeNode(x,yp,zp,m_uiMaxDepth,m_uiDimension,m_uiMaxDepth);

          vertices[7] = TreeNode(xp,yp,zp,m_uiMaxDepth,m_uiDimension,m_uiMaxDepth);
        }//end init Vtx Block

        switch (childNum) {
          case 0:
            {     
              //0 and 7 can't be hanging
              //+x,+y+z
              switch( hnMask ) {
                case  cNumEtype<0>::ET_Y: 
                  {
                    //2y
                    vertices[2] = TreeNode(x ,yp + sz ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
                    break;
                  }
                case  cNumEtype<0>::ET_X:
                  {
                    //1x
                    vertices[1] = TreeNode(xp + sz, y ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
                    break;
                  }
                case  cNumEtype<0>::ET_XY:
                  {
                    //2y 1x
                    vertices[1] = TreeNode(xp + sz, y ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    vertices[2] = TreeNode(x ,yp + sz ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
                    break;
                  }
                case  cNumEtype<0>::ET_Z:
                  {
                    //4z
                    vertices[4] = TreeNode(x,
                        y,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<0>::ET_ZY:
                  {
                    //4z 2y
                    vertices[4] = TreeNode(x,
                        y,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[2] = TreeNode(x ,yp + sz ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<0>::ET_ZX:
                  {
                    //4z 1x
                    vertices[4] = TreeNode(x,
                        y,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[1] = TreeNode(xp + sz, y ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);


                    break;
                  }
                case  cNumEtype<0>::ET_ZXY:
                  {
                    //4z 2y 1x
                    vertices[4] = TreeNode(x,
                        y,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[1] = TreeNode(xp + sz, y ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    vertices[2] = TreeNode(x ,yp + sz ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<0>::ET_XY_XY:
                  { 
                    //3xy 2y 1x
                    vertices[3] = TreeNode(xp + sz,
                        yp + sz,z, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[1] = TreeNode(xp + sz, y ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    vertices[2] = TreeNode(x ,yp + sz ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<0>::ET_XY_ZXY:
                  { 
                    //3xy 4z 2y 1x
                    vertices[3] = TreeNode(xp + sz,
                        yp + sz,z, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[4] = TreeNode(x,
                        y,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[1] = TreeNode(xp + sz, y ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    vertices[2] = TreeNode(x ,yp + sz ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<0>::ET_YZ_ZY: 
                  {
                    //6yz 4z 2y
                    vertices[6] = TreeNode(x,
                        yp + sz,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[4] = TreeNode(x,
                        y,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[2] = TreeNode(x ,yp + sz ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<0>::ET_YZ_ZXY: 
                  {
                    //6yz 4z 2y 1x
                    vertices[6] = TreeNode(x,
                        yp + sz,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[4] = TreeNode(x,
                        y,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[1] = TreeNode(xp + sz, y ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    vertices[2] = TreeNode(x ,yp + sz ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<0>::ET_YZ_XY_ZXY: 
                  {
                    //6yz 3xy 4z 2y 1x
                    vertices[6] = TreeNode(x,
                        yp + sz,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[3] = TreeNode(xp + sz,
                        yp + sz,z, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[4] = TreeNode(x,
                        y,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[1] = TreeNode(xp + sz, y ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    vertices[2] = TreeNode(x ,yp + sz ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<0>::ET_ZX_ZX: 
                  {
                    //5zx 4z 1x
                    vertices[5] = TreeNode(xp + sz,
                        y,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[4] = TreeNode(x,
                        y,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[1] = TreeNode(xp + sz, y ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);


                    break;
                  }
                case  cNumEtype<0>::ET_ZX_ZXY: 
                  {
                    //5zx 4z 2y 1x
                    vertices[5] = TreeNode(xp + sz,
                        y,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[4] = TreeNode(x,
                        y,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[1] = TreeNode(xp + sz, y ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    vertices[2] = TreeNode(x ,yp + sz ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<0>::ET_ZX_XY_ZXY: 
                  {
                    //5zx 3xy 4z 2y 1x
                    vertices[5] = TreeNode(xp + sz,
                        y,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[3] = TreeNode(xp + sz,
                        yp + sz,z, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[4] = TreeNode(x,
                        y,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[1] = TreeNode(xp + sz, y ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    vertices[2] = TreeNode(x ,yp + sz ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<0>::ET_ZX_YZ_ZXY: 
                  {
                    //5zx 6yz 4z 2y 1x
                    vertices[5] = TreeNode(xp + sz,
                        y,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[6] = TreeNode(x,
                        yp + sz,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[4] = TreeNode(x,
                        y,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[1] = TreeNode(xp + sz, y ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    vertices[2] = TreeNode(x ,yp + sz ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<0>::ET_ZX_YZ_XY_ZXY: 
                  {
                    //5zx 6yz 3xy 
                    //4z 2y 1x
                    vertices[5] = TreeNode(xp + sz,
                        y,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[6] = TreeNode(x,
                        yp + sz,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[3] = TreeNode(xp + sz,
                        yp + sz,z, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[4] = TreeNode(x,
                        y,zp + sz, m_uiMaxDepth,
                        m_uiDimension, m_uiMaxDepth);

                    vertices[1] = TreeNode(xp + sz, y ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    vertices[2] = TreeNode(x ,yp + sz ,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                default:
                  assert(false);

              }//end 17 types
              break;
            }
          case 1: 
            {
              //1 and 6 can't be hanging
              //-x,+y,+z
              switch( hnMask ) {
                case  cNumEtype<1>::ET_Y:  
                  {
                    //0x
                    vertices[0] = TreeNode(x - sz,y,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<1>::ET_X: 
                  {
                    //3y
                    vertices[3] = TreeNode(xp,yp + sz,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<1>::ET_XY: 
                  {
                    //0x 3y
                    vertices[0] = TreeNode(x - sz,y,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    vertices[3] = TreeNode(xp,yp + sz,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<1>::ET_Z: 
                  {
                    //5z
                    vertices[5] = TreeNode(xp,y,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<1>::ET_ZY: 
                  {
                    //5z 0x
                    vertices[0] = TreeNode(x - sz,y,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<1>::ET_ZX: 
                  {
                    //5z 3y
                    vertices[5] = TreeNode(xp,y,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[3] = TreeNode(xp,yp + sz,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<1>::ET_ZXY: 
                  {
                    //5z 3y 0x
                    vertices[5] = TreeNode(xp,y,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x - sz,y,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    vertices[3] = TreeNode(xp,yp + sz,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<1>::ET_XY_XY:
                  { 
                    //2xy 3y 0x
                    vertices[2] = TreeNode(x - sz,yp + sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x - sz,y,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    vertices[3] = TreeNode(xp,yp + sz,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<1>::ET_XY_ZXY:
                  { 
                    //2xy 5z 3y 0x
                    vertices[2] = TreeNode(x - sz,yp + sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x - sz,y,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    vertices[3] = TreeNode(xp,yp + sz,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<1>::ET_YZ_ZY: 
                  {
                    //4xz 5z 0x
                    vertices[4] = TreeNode(x - sz,y,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x - sz,y,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<1>::ET_YZ_ZXY: 
                  {
                    //4xz 5z 3y 0x
                    vertices[4] = TreeNode(x - sz,y,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x - sz,y,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    vertices[3] = TreeNode(xp,yp + sz,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<1>::ET_YZ_XY_ZXY: 
                  {
                    //4xz 2xy 5z 3y 0x
                    vertices[4] = TreeNode(x - sz,y,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x - sz,yp + sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x - sz,y,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    vertices[3] = TreeNode(xp,yp + sz,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<1>::ET_ZX_ZX: 
                  {
                    //7yz 5z 3y
                    vertices[7] = TreeNode(xp,yp + sz,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[3] = TreeNode(xp,yp + sz,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<1>::ET_ZX_ZXY: 
                  {
                    //7yz 5z 3y 0x
                    vertices[7] = TreeNode(xp,yp + sz,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x - sz,y,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    vertices[3] = TreeNode(xp,yp + sz,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<1>::ET_ZX_XY_ZXY: 
                  {
                    //7yz 2xy 5z 3y 0x
                    vertices[7] = TreeNode(xp,yp + sz,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x - sz,yp + sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x - sz,y,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    vertices[3] = TreeNode(xp,yp + sz,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<1>::ET_ZX_YZ_ZXY: 
                  {
                    //7yz 4xz 5z 3y 0x
                    vertices[7] = TreeNode(xp,yp + sz,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[4] = TreeNode(x - sz,y,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x - sz,y,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    vertices[3] = TreeNode(xp,yp + sz,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<1>::ET_ZX_YZ_XY_ZXY: 
                  {
                    //7yz 4xz 2xy 5z 3y 0x
                    vertices[7] = TreeNode(xp,yp + sz,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[4] = TreeNode(x - sz,y,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x - sz,yp + sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x - sz,y,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    vertices[3] = TreeNode(xp,yp + sz,
                        z , m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);

                    break;
                  }
                default:
                  assert(false);

              }//end 17 types

              break;
            }
          case 2: 
            {
              //2 and 5 can't be hanging
              //+x,-y,+z
              switch( hnMask ) {
                case  cNumEtype<2>::ET_Y:  
                  {
                    //3x
                    vertices[3] = TreeNode(xp + sz,yp,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<2>::ET_X: 
                  {
                    //0y
                    vertices[0] = TreeNode(x,y - sz,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<2>::ET_XY: 
                  {
                    //0y 3x
                    vertices[0] = TreeNode(x,y - sz,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[3] = TreeNode(xp + sz,yp,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<2>::ET_Z: 
                  {
                    //6z
                    vertices[6] = TreeNode(x,yp,
                        zp + sz , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<2>::ET_ZY: 
                  {
                    //6z 3x
                    vertices[6] = TreeNode(x,yp,
                        zp + sz , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[3] = TreeNode(xp + sz,yp,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<2>::ET_ZX: 
                  {
                    //6z 0y
                    vertices[6] = TreeNode(x,yp,
                        zp + sz , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y - sz,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<2>::ET_ZXY: 
                  {
                    //6z 0y 3x
                    vertices[6] = TreeNode(x,yp,
                        zp + sz , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y - sz,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[3] = TreeNode(xp + sz,yp,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<2>::ET_XY_XY:
                  { 
                    //1xy 0y 3x
                    vertices[1] = TreeNode(xp + sz,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y - sz,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[3] = TreeNode(xp + sz,yp,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<2>::ET_XY_ZXY:
                  { 
                    //1xy 6z 0y 3x
                    vertices[1] = TreeNode(xp + sz,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x,yp,
                        zp + sz , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y - sz,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[3] = TreeNode(xp + sz,yp,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<2>::ET_YZ_ZY: 
                  {
                    //7xz 6z 3x
                    vertices[7] = TreeNode(xp + sz,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x,yp,
                        zp + sz , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[3] = TreeNode(xp + sz,yp,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<2>::ET_YZ_ZXY: 
                  {
                    //7xz 6z 0y 3x
                    vertices[7] = TreeNode(xp + sz,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x,yp,
                        zp + sz , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y - sz,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[3] = TreeNode(xp + sz,yp,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<2>::ET_YZ_XY_ZXY: 
                  {
                    //7xz 1xy 6z 0y 3x
                    vertices[7] = TreeNode(xp + sz,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp + sz,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x,yp,
                        zp + sz , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y - sz,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[3] = TreeNode(xp + sz,yp,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<2>::ET_ZX_ZX: 
                  {
                    //4zy 6z 0y
                    vertices[4] = TreeNode(x,y - sz,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x,yp,
                        zp + sz , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y - sz,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<2>::ET_ZX_ZXY: 
                  {
                    //4zy 6z 0y 3x
                    vertices[4] = TreeNode(x,y - sz,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x,yp,
                        zp + sz , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y - sz,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[3] = TreeNode(xp + sz,yp,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<2>::ET_ZX_XY_ZXY: 
                  {
                    //4zy 1xy 6z 0y 3x
                    vertices[4] = TreeNode(x,y - sz,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp + sz,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x,yp,
                        zp + sz , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y - sz,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[3] = TreeNode(xp + sz,yp,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<2>::ET_ZX_YZ_ZXY: 
                  {
                    //4zy 7xz 6z 0y 3x
                    vertices[4] = TreeNode(x,y - sz,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp + sz,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x,yp,
                        zp + sz , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y - sz,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[3] = TreeNode(xp + sz,yp,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<2>::ET_ZX_YZ_XY_ZXY: 
                  {
                    //4zy 7xz 1xy 6z 0y 3x
                    vertices[4] = TreeNode(x,y - sz,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp + sz,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp + sz,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x,yp,
                        zp + sz , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y - sz,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[3] = TreeNode(xp + sz,yp,
                        z , m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                default:
                  assert(false);

              }//end 17 types


              break;
            }
          case 3: 
            {
              //3 and 4 can't be hanging
              //-x,-y,+z
              switch( hnMask ) {
                case  cNumEtype<3>::ET_Y:  
                  {
                    //1y
                    vertices[1] = TreeNode(xp,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<3>::ET_X: 
                  {
                    //2x
                    vertices[2] = TreeNode(x - sz,yp,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<3>::ET_XY: 
                  {
                    //2x 1y
                    vertices[2] = TreeNode(x - sz,yp,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<3>::ET_Z: 
                  {
                    //7z
                    vertices[7] = TreeNode(xp,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<3>::ET_ZY: 
                  {
                    //7z 1y
                    vertices[7] = TreeNode(xp,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<3>::ET_ZX: 
                  {
                    //7z 2x
                    vertices[7] = TreeNode(xp,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x - sz,yp,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<3>::ET_ZXY: 
                  {
                    //7z 2x 1y
                    vertices[7] = TreeNode(xp,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x - sz,yp,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<3>::ET_XY_XY:
                  { 
                    //0xy 2x 1y
                    vertices[0] = TreeNode(x - sz,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x - sz,yp,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<3>::ET_XY_ZXY:
                  { 
                    //0xy 7z 2x 1y
                    vertices[0] = TreeNode(x - sz,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x - sz,yp,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<3>::ET_YZ_ZY: 
                  {
                    //5yz 7z 1y
                    vertices[5] = TreeNode(xp,y - sz,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<3>::ET_YZ_ZXY: 
                  {
                    //5yz 7z 2x 1y
                    vertices[5] = TreeNode(xp,y - sz,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x - sz,yp,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<3>::ET_YZ_XY_ZXY: 
                  {
                    //5yz 0xy 7z 2x 1y
                    vertices[5] = TreeNode(xp,y - sz,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x - sz,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x - sz,yp,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<3>::ET_ZX_ZX: 
                  {
                    //6xz 7z 2x
                    vertices[6] = TreeNode(x - sz,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x - sz,yp,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<3>::ET_ZX_ZXY: 
                  {
                    //6xz 7z 2x 1y
                    vertices[6] = TreeNode(x - sz,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x - sz,yp,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<3>::ET_ZX_XY_ZXY: 
                  {
                    //6xz 0xy 7z 2x 1y
                    vertices[6] = TreeNode(x - sz,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x - sz,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x - sz,yp,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<3>::ET_ZX_YZ_ZXY: 
                  {
                    //6xz 5yz 7z 2x 1y
                    vertices[6] = TreeNode(x - sz,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y - sz,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x - sz,yp,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<3>::ET_ZX_YZ_XY_ZXY: 
                  {
                    //6xz 5yz 0xy 7z 2x 1y
                    vertices[6] = TreeNode(x - sz,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y - sz,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x - sz,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp,yp,
                        zp + sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x - sz,yp,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y - sz,
                        z, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                default:
                  assert(false);

              }//end 17 types


              break;
            }
          case 4: 
            {
              //4 and 3 can't be hanging
              //+x,+y,-z
              switch( hnMask ) {
                case  cNumEtype<4>::ET_Y:  
                  {
                    //0z
                    vertices[0] = TreeNode(x,y,
                        z - sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<4>::ET_X: 
                  {
                    //5x
                    vertices[5] = TreeNode(xp+sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<4>::ET_XY: 
                  {
                    //5x 0z
                    vertices[5] = TreeNode(xp+sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y,
                        z - sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<4>::ET_Z: 
                  {
                    //6y
                    vertices[6] = TreeNode(x,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<4>::ET_ZY: 
                  {
                    //6y 0z
                    vertices[6] = TreeNode(x,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y,
                        z - sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<4>::ET_ZX: 
                  {
                    //6y 5x
                    vertices[6] = TreeNode(x,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp+sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<4>::ET_ZXY: 
                  {
                    //6y 5x 0z
                    vertices[6] = TreeNode(x,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp+sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y,
                        z - sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<4>::ET_XY_XY:
                  { 
                    //1xz 5x 0z
                    vertices[1] = TreeNode(xp+sz,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp+sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<4>::ET_XY_ZXY:
                  { 
                    //1xz 6y 5x 0z
                    vertices[1] = TreeNode(xp+sz,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp+sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y,
                        z - sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<4>::ET_YZ_ZY: 
                  {
                    //2yz 6y 0z
                    vertices[2] = TreeNode(x,yp+sz,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y,
                        z - sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<4>::ET_YZ_ZXY: 
                  {
                    //2yz 6y 5x 0z
                    vertices[2] = TreeNode(x,yp+sz,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp+sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y,
                        z - sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<4>::ET_YZ_XY_ZXY: 
                  {
                    //2yz 1xz 6y 5x 0z
                    vertices[2] = TreeNode(x,yp+sz,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp+sz,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp+sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y,
                        z - sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<4>::ET_ZX_ZX: 
                  {
                    //7xy 6y 5x
                    vertices[7] = TreeNode(xp+sz,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp+sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<4>::ET_ZX_ZXY: 
                  {
                    //7xy 6y 5x 0z
                    vertices[7] = TreeNode(xp+sz,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp+sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y,
                        z - sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<4>::ET_ZX_XY_ZXY: 
                  {
                    //7xy 1xz 6y 5x 0z
                    vertices[7] = TreeNode(xp+sz,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp+sz,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp+sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y,
                        z - sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<4>::ET_ZX_YZ_ZXY: 
                  {
                    //7xy 2yz 6y 5x 0z
                    vertices[7] = TreeNode(xp+sz,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x,yp+sz,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp+sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y,
                        z - sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<4>::ET_ZX_YZ_XY_ZXY: 
                  {
                    //7xy 2yz 1xz 6y 5x 0z
                    vertices[7] = TreeNode(xp+sz,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x,yp+sz,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp+sz,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp+sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y,
                        z - sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                default:
                  assert(false);

              }//end 17 types


              break;
            }
          case 5: 
            {
              //5 and 2 can't be hanging
              //-x,+y,-z
              switch( hnMask ) {
                case  cNumEtype<5>::ET_Y:  
                  {
                    //1z
                    vertices[1] = TreeNode(xp,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<5>::ET_X: 
                  {
                    //7y
                    vertices[7] = TreeNode(xp,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<5>::ET_XY: 
                  {
                    //7y 1z
                    vertices[7] = TreeNode(xp,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<5>::ET_Z: 
                  {
                    //4x
                    vertices[4] = TreeNode(x-sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<5>::ET_ZY: 
                  {
                    //4x 1z
                    vertices[4] = TreeNode(x-sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<5>::ET_ZX: 
                  {
                    //4x 7y
                    vertices[4] = TreeNode(x-sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<5>::ET_ZXY: 
                  {
                    //4x 7y 1z
                    vertices[4] = TreeNode(x-sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<5>::ET_XY_XY:
                  { 
                    //3yz 7y 1z
                    vertices[3] = TreeNode(xp,yp+sz,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<5>::ET_XY_ZXY:
                  { 
                    //3yz 4x 7y 1z
                    vertices[3] = TreeNode(xp,yp+sz,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[4] = TreeNode(x-sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<5>::ET_YZ_ZY: 
                  {
                    //0xz 4x 1z
                    vertices[4] = TreeNode(x-sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x-sz,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<5>::ET_YZ_ZXY: 
                  {
                    //0xz 4x 7y 1z
                    vertices[4] = TreeNode(x-sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x-sz,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<5>::ET_YZ_XY_ZXY: 
                  {
                    //0xz 3yz 4x 7y 1z
                    vertices[3] = TreeNode(xp,yp+sz,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[4] = TreeNode(x-sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x-sz,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<5>::ET_ZX_ZX: 
                  {
                    //6xy 4x 7y
                    vertices[4] = TreeNode(x-sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x-sz,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<5>::ET_ZX_ZXY: 
                  {
                    //6xy 4x 7y 1z
                    vertices[4] = TreeNode(x-sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x-sz,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<5>::ET_ZX_XY_ZXY: 
                  {
                    //6xy 3yz 4x 7y 1z
                    vertices[3] = TreeNode(xp,yp+sz,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[4] = TreeNode(x-sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x-sz,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<5>::ET_ZX_YZ_ZXY: 
                  {
                    //6xy 0xz 4x 7y 1z
                    vertices[4] = TreeNode(x-sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x-sz,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x-sz,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<5>::ET_ZX_YZ_XY_ZXY: 
                  {
                    //6xy 0xz 3yz 4x 7y 1z
                    vertices[3] = TreeNode(xp,yp+sz,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[4] = TreeNode(x-sz,y,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x-sz,yp+sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x-sz,y,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                default:
                  assert(false);

              }//end 17 types


              break;
            }
          case 6: 
            {
              //6 and 1 can't be hanging
              //+x,-y,-z
              switch( hnMask ) {
                case  cNumEtype<6>::ET_Y:  
                  {
                    //2z
                    vertices[2] = TreeNode(x,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<6>::ET_X: 
                  {
                    //4y
                    vertices[4] = TreeNode(x,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<6>::ET_XY: 
                  {
                    //2z 4y
                    vertices[4] = TreeNode(x,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<6>::ET_Z: 
                  {
                    //7x
                    vertices[7] = TreeNode(xp+sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<6>::ET_ZY: 
                  {
                    //7x 2z
                    vertices[7] = TreeNode(xp+sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<6>::ET_ZX: 
                  {
                    //7x 4y
                    vertices[4] = TreeNode(x,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp+sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<6>::ET_ZXY: 
                  {
                    //7x 4y 2z
                    vertices[4] = TreeNode(x,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp+sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<6>::ET_XY_XY:
                  {
                    //0yz 4y 2z
                    vertices[4] = TreeNode(x,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y-sz,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<6>::ET_XY_ZXY:
                  { 
                    //0yz 7x 4y 2z
                    vertices[4] = TreeNode(x,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp+sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y-sz,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<6>::ET_YZ_ZY: 
                  {
                    //3xz 7x 2z
                    vertices[3] = TreeNode(xp+sz,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp+sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<6>::ET_YZ_ZXY: 
                  {
                    //3xz 7x 4y 2z
                    vertices[3] = TreeNode(xp+sz,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[4] = TreeNode(x,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp+sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<6>::ET_YZ_XY_ZXY: 
                  {
                    //3xz 0yz 7x 4y 2z
                    vertices[3] = TreeNode(xp+sz,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[4] = TreeNode(x,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp+sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y-sz,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<6>::ET_ZX_ZX: 
                  {
                    //5xy 7x 4y
                    vertices[4] = TreeNode(x,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp+sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp+sz,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<6>::ET_ZX_ZXY: 
                  {
                    //5xy 7x 4y 2z
                    vertices[4] = TreeNode(x,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp+sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp+sz,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<6>::ET_ZX_XY_ZXY: 
                  {
                    //5xy 0yz 7x 4y 2z
                    vertices[4] = TreeNode(x,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp+sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp+sz,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y-sz,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<6>::ET_ZX_YZ_ZXY: 
                  {
                    //5xy 3xz 7x 4y 2z
                    vertices[3] = TreeNode(xp+sz,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[4] = TreeNode(x,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp+sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp+sz,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<6>::ET_ZX_YZ_XY_ZXY: 
                  {
                    //5xy 3xz 0yz 7x 4y 2z
                    vertices[3] = TreeNode(xp+sz,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[4] = TreeNode(x,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[7] = TreeNode(xp+sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp+sz,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[0] = TreeNode(x,y-sz,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                default:
                  assert(false);

              }//end 17 types

              break;
            }
          case 7: 
            {
              //0 and 7 can't be hanging			
              //-x,-y,-z
              switch( hnMask ) {
                case  cNumEtype<7>::ET_Y:  
                  {
                    //3z
                    vertices[3] = TreeNode(xp,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<7>::ET_X: 
                  {
                    //6x
                    vertices[6] = TreeNode(x-sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<7>::ET_XY: 
                  {
                    //6x 3z
                    vertices[3] = TreeNode(xp,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x-sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<7>::ET_Z: 
                  {
                    //5y
                    vertices[5] = TreeNode(xp,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<7>::ET_ZY: 
                  {
                    //5y 3z
                    vertices[3] = TreeNode(xp,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<7>::ET_ZX: 
                  {
                    //5y 6x
                    vertices[6] = TreeNode(x-sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<7>::ET_ZXY: 
                  {
                    //5y 6x 3z
                    vertices[3] = TreeNode(xp,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x-sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<7>::ET_XY_XY:
                  { 
                    //2xz 6x 3z
                    vertices[3] = TreeNode(xp,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x-sz,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x-sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<7>::ET_XY_ZXY:
                  { 
                    //2xz 5y 6x 3z
                    vertices[3] = TreeNode(xp,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x-sz,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x-sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<7>::ET_YZ_ZY: 
                  {
                    //1yz 5y 3z
                    vertices[3] = TreeNode(xp,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y-sz,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<7>::ET_YZ_ZXY: 
                  {
                    //1yz 5y 6x 3z
                    vertices[3] = TreeNode(xp,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y-sz,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x-sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<7>::ET_YZ_XY_ZXY: 
                  {
                    //1yz 2xz 5y 6x 3z
                    vertices[3] = TreeNode(xp,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x-sz,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y-sz,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x-sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<7>::ET_ZX_ZX: 
                  {
                    //4xy 5y 6x
                    vertices[4] = TreeNode(x-sz,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x-sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<7>::ET_ZX_ZXY: 
                  {
                    //4xy 5y 6x 3z
                    vertices[3] = TreeNode(xp,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[4] = TreeNode(x-sz,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x-sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<7>::ET_ZX_XY_ZXY: 
                  {
                    //4xy 2xz 5y 6x 3z
                    vertices[3] = TreeNode(xp,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[4] = TreeNode(x-sz,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x-sz,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x-sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<7>::ET_ZX_YZ_ZXY: 
                  {
                    //4xy 1yz 5y 6x 3z
                    vertices[3] = TreeNode(xp,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[4] = TreeNode(x-sz,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y-sz,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x-sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                case  cNumEtype<7>::ET_ZX_YZ_XY_ZXY: 
                  {
                    //4xy 1yz 2xz 5y 6x 3z
                    vertices[3] = TreeNode(xp,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[4] = TreeNode(x-sz,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[2] = TreeNode(x-sz,yp,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[1] = TreeNode(xp,y-sz,
                        z-sz, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[6] = TreeNode(x-sz,yp,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    vertices[5] = TreeNode(xp,y-sz,
                        zp, m_uiMaxDepth, m_uiDimension,
                        m_uiMaxDepth);

                    break;
                  }
                default:
                  assert(false);

              }//end 17 types

              break; 
            }
          default:
            assert(false);
        }//end case childnum

        //init idxHolders
        std::vector<seq::IndexHolder<ot::TreeNode> > idxHolder(8);
        {
          idxHolder[0].value = &vertices[0];
          idxHolder[0].index = 0;

          idxHolder[1].value = &vertices[1];
          idxHolder[1].index = 1;

          idxHolder[2].value = &vertices[2];
          idxHolder[2].index = 2;

          idxHolder[3].value = &vertices[3];
          idxHolder[3].index = 3;

          idxHolder[4].value = &vertices[4];
          idxHolder[4].index = 4;

          idxHolder[5].value = &vertices[5];
          idxHolder[5].index = 5;

          idxHolder[6].value = &vertices[6];
          idxHolder[6].index = 6;

          idxHolder[7].value = &vertices[7];
          idxHolder[7].index = 7;
        }//end idxHolders Block

        sort < seq::IndexHolder<ot::TreeNode> *> (&(*idxHolder.begin()), 8); 

        for (int k = 0; k < 8; k++) {
          nodes[idxHolder[k].index] = nn[k];
        }

      } else {
        //No hanging nodes...
        nodes[0] = nn[0];
        switch(m_ucpSortOrdersPtr[m_uiCurrent]) {
          case ot::DA_FLAGS::ZYX: {
                                    nodes[4] = nn[1];
                                    nodes[2] = nn[2];
                                    nodes[6] = nn[3];
                                    nodes[1] = nn[4];
                                    nodes[5] = nn[5];
                                    nodes[3] = nn[6];
                                    break;
                                  }
          case ot::DA_FLAGS::YZX: {
                                    nodes[2] = nn[1];
                                    nodes[4] = nn[2];
                                    nodes[6] = nn[3];
                                    nodes[1] = nn[4];
                                    nodes[3] = nn[5];
                                    nodes[5] = nn[6];
                                    break;
                                  }
          case ot::DA_FLAGS::YXZ: {
                                    nodes[2] = nn[1];
                                    nodes[1] = nn[2];
                                    nodes[3] = nn[3];
                                    nodes[4] = nn[4];
                                    nodes[6] = nn[5];
                                    nodes[5] = nn[6];
                                    break;
                                  }
          case ot::DA_FLAGS::ZXY: {
                                    nodes[4] = nn[1];
                                    nodes[1] = nn[2];
                                    nodes[5] = nn[3];
                                    nodes[2] = nn[4];
                                    nodes[6] = nn[5];
                                    nodes[3] = nn[6];
                                    break;
                                  }
          case ot::DA_FLAGS::XZY: {
                                    nodes[1] = nn[1];
                                    nodes[4] = nn[2];
                                    nodes[5] = nn[3];
                                    nodes[2] = nn[4];
                                    nodes[3] = nn[5];
                                    nodes[6] = nn[6];
                                    break;
                                  }
          case ot::DA_FLAGS::XYZ: {
                                    nodes[1] = nn[1];
                                    nodes[2] = nn[2];
                                    nodes[3] = nn[3];
                                    nodes[4] = nn[4];
                                    nodes[5] = nn[5];
                                    nodes[6] = nn[6];
                                    break;
                                  }
        }//end switch-case order type
        nodes[7] = nn[7];
      }//end if hanging
#ifdef __DEBUG_DA__
      for ( unsigned int i=0; i<8; i++) {
        if (nodes[i] >= m_uiLocalBufferSize) {
          std::cout << m_iRankActive << RED" Node index is out of bounds( "<<nodes[i]<<" ) for element "NRM<< m_uiCurrent << std::endl;
          std::cout << m_iRankActive << " preElem: "<<m_uiPreGhostElementSize<<" elemBeg: "<<m_uiElementBegin<<" elemEnd: "<<m_uiElementEnd<<" postBeg: "<<m_uiPostGhostBegin<<" bufSz: "<<m_uiLocalBufferSize << std::endl;
          std::cout<<"CNum: "<<this->getChildNumber()<<std::endl;
          std::cout<<"nlist after decompression..."<<std::endl;
          for(unsigned int j=0; j< 8; j++) {
            std::cout<<nodes[j]<<", ";
          }
          std::cout<<std::endl;
          assert(false);
        }
      }
#endif
    }//end if compress
    return(0);
  }//end function

  template <typename T>
    void injectNodalVector(ot::DA* dac, ot::DA* daf, unsigned int dof,
        std::vector<T>& fVec, std::vector<T>& cVec, void (*setZero)(T&)) {

      dac->createVector<T>(cVec, false, false, dof);

      T* cArr = NULL;
      T* fArr = NULL;
      dac->vecGetBuffer<T>(cVec, cArr, false, false, false, dof);
      daf->vecGetBuffer<T>(fVec, fArr, false, false, true, dof);

      for(int i = 0; i < (dof*(dac->getLocalBufferSize())); i++) {
        (*setZero)(cArr[i]);
      }

      daf->ReadFromGhostsBegin<T>(fArr, dof);
      daf->ReadFromGhostsEnd<T>(fArr);

      if(dac->iAmActive()) {
        for(dac->init<ot::DA_FLAGS::WRITABLE>(), daf->init<ot::DA_FLAGS::WRITABLE>();
            dac->curr() < dac->end<ot::DA_FLAGS::WRITABLE>();
            dac->next<ot::DA_FLAGS::WRITABLE>()) {
          unsigned int idxC = dac->curr();
          Point cPt = dac->getCurrentOffset();
          Point fPt = daf->getCurrentOffset();
          assert(cPt == fPt);
          unsigned char currentFlags;
          unsigned char hnMask = dac->getHangingNodeIndex(idxC);
          bool isBoundary = dac->isBoundaryOctant(&currentFlags);
          unsigned int cIndices[8];
          if(currentFlags > ot::TreeNode::NEG_POS_DEMARCATION) {
            //has atleast one positive boundary
            dac->getNodeIndices(cIndices);
          } else {
            if(dac->isLUTcompressed()) {
              dac->updateQuotientCounter();
            }
          }
          int xPosBdy = (currentFlags & ot::TreeNode::X_POS_BDY);
          int yPosBdy = (currentFlags & ot::TreeNode::Y_POS_BDY);
          int zPosBdy = (currentFlags & ot::TreeNode::Z_POS_BDY);
          if(daf->getLevel(daf->curr()) == dac->getLevel(idxC)) {
            unsigned int idxF = daf->curr();
            unsigned int fIndices[8];
            if(currentFlags > ot::TreeNode::NEG_POS_DEMARCATION) {
              daf->getNodeIndices(fIndices);
            } else {
              if(daf->isLUTcompressed()) {
                daf->updateQuotientCounter();
              }
            }
            if(xPosBdy) {
              if(!(hnMask & (1<<1))) {
                for(int d = 0; d < dof; d++) {
                  cArr[(dof*cIndices[1]) + d] = fArr[(dof*fIndices[1]) + d];
                }
              }
            }
            if(yPosBdy) {
              if(!(hnMask & (1<<2))) {
                for(int d = 0; d < dof; d++) {
                  cArr[(dof*cIndices[2]) + d] = fArr[(dof*fIndices[2]) + d];
                }
              }
            }
            if(zPosBdy) {
              if(!(hnMask & (1<<4))) {
                for(int d = 0; d < dof; d++) {
                  cArr[(dof*cIndices[4]) + d] = fArr[(dof*fIndices[4]) + d];
                }
              }
            }
            if(xPosBdy && yPosBdy) {
              if(!(hnMask & (1<<3))) {
                for(int d = 0; d < dof; d++) {
                  cArr[(dof*cIndices[3]) + d] = fArr[(dof*fIndices[3]) + d];
                }
              }
            }
            if(xPosBdy && zPosBdy) {
              if(!(hnMask & (1<<5))) {
                for(int d = 0; d < dof; d++) {
                  cArr[(dof*cIndices[5]) + d] = fArr[(dof*fIndices[5]) + d];
                }
              }
            }
            if(yPosBdy && zPosBdy) {
              if(!(hnMask & (1<<6))) {
                for(int d = 0; d < dof; d++) {
                  cArr[(dof*cIndices[6]) + d] = fArr[(dof*fIndices[6]) + d];
                }
              }
            }
            if(xPosBdy && yPosBdy && zPosBdy) {
              if(!(hnMask & (1<<7))) {
                for(int d = 0; d < dof; d++) {
                  cArr[(dof*cIndices[7]) + d] = fArr[(dof*fIndices[7]) + d];
                }
              }
            }
            if(!(hnMask & 1)) {
              //Anchor is not hanging
              for(int d = 0; d < dof; d++) {
                cArr[(dof*idxC) + d] = fArr[(dof*idxF) + d];
              }
            }
            daf->next<ot::DA_FLAGS::WRITABLE>();
          } else {
            for(unsigned int cNumFine = 0; cNumFine < 8; cNumFine++) {
              unsigned int idxF = daf->curr();
              unsigned int fIndices[8];
              if(currentFlags > ot::TreeNode::NEG_POS_DEMARCATION) {
                daf->getNodeIndices(fIndices);
              } else {
                if(daf->isLUTcompressed()) {
                  daf->updateQuotientCounter();
                }
              }
              if(xPosBdy) {
                if(cNumFine == 1) {
                  if(!(hnMask & (1<<1))) {
                    for(int d = 0; d < dof; d++) {
                      cArr[(dof*cIndices[1]) + d] = fArr[(dof*fIndices[1]) + d];
                    }
                  }
                }
              }
              if(yPosBdy) {
                if(cNumFine == 2) {
                  if(!(hnMask & (1<<2))) {
                    for(int d = 0; d < dof; d++) {
                      cArr[(dof*cIndices[2]) + d] = fArr[(dof*fIndices[2]) + d];
                    }
                  }
                }
              }
              if(zPosBdy) {
                if(cNumFine == 4) {
                  if(!(hnMask & (1<<4))) {
                    for(int d = 0; d < dof; d++) {
                      cArr[(dof*cIndices[4]) + d] = fArr[(dof*fIndices[4]) + d];
                    }
                  }
                }
              }
              if(xPosBdy && yPosBdy) {
                if(cNumFine == 3) {
                  if(!(hnMask & (1<<3))) {
                    for(int d = 0; d < dof; d++) {
                      cArr[(dof*cIndices[3]) + d] = fArr[(dof*fIndices[3]) + d];
                    }
                  }
                }
              }
              if(xPosBdy && zPosBdy) {
                if(cNumFine == 5) {
                  if(!(hnMask & (1<<5))) {
                    for(int d = 0; d < dof; d++) {
                      cArr[(dof*cIndices[5]) + d] = fArr[(dof*fIndices[5]) + d];
                    }
                  }
                }
              }
              if(yPosBdy && zPosBdy) {
                if(cNumFine == 6) {
                  if(!(hnMask & (1<<6))) {
                    for(int d = 0; d < dof; d++) {
                      cArr[(dof*cIndices[6]) + d] = fArr[(dof*fIndices[6]) + d];
                    }
                  }
                }
              }
              if(xPosBdy && yPosBdy && zPosBdy) {
                if(cNumFine == 7) {
                  if(!(hnMask & (1<<7))) {
                    for(int d = 0; d < dof; d++) {
                      cArr[(dof*cIndices[7]) + d] = fArr[(dof*fIndices[7]) + d];
                    }
                  }
                }
              }
              if(cNumFine == 0) {
                if(!(hnMask & 1)) {
                  //Anchor is not hanging
                  for(int d = 0; d < dof; d++) {
                    cArr[(dof*idxC) + d] = fArr[(dof*idxF) + d];
                  }
                }
              }
              daf->next<ot::DA_FLAGS::WRITABLE>();
            }
          }
        }
      }

      dac->WriteToGhostsBegin<T>(cArr, dof);
      dac->WriteToGhostsEnd<T>(cArr, dof);

      dac->vecRestoreBuffer<T>(cVec, cArr, false, false, false, dof);
      daf->vecRestoreBuffer<T>(fVec, fArr, false, false, true, dof);
    }

} // end namespace

