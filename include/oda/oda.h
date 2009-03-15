
/**
  @file oda.h
  @brief 		The class that manages the octree mesh that supports
  trilinear shape functions.
  @author		Hari Sundar, hsundar@gmail.com
  @author		Rahul S. Sampath, rahul.sampath@gmail.com
  */ 

#ifndef __OCT_DA_H__
#define __OCT_DA_H__

#include "mpi.h"
#include "matRecord.h"
#include "loopCounters.h"
#include "updateCtx.h"
#include "odaUtils.h"
#include "cnumEtypes.h"
#include "Point.h"
#include <vector>
#include "petsc.h"
#include "petscvec.h"
#include "petscmat.h"
#include "dendro.h"

#ifndef iC
#define iC(fun) {CHKERRQ(fun);}
#endif

#ifdef __DEBUG__
#ifndef __DEBUG_DA__
#define __DEBUG_DA__
#endif
#endif

#ifdef PETSC_USE_LOG

#include "petscsys.h"

namespace ot {
  extern int readFromGhostNodesBeginEvent;
  extern int readFromGhostNodesEndEvent;
  extern int readFromGhostElemsBeginEvent;
  extern int readFromGhostElemsEndEvent;
  extern int writeToGhostNodesBeginEvent;
  extern int writeToGhostNodesEndEvent;
  extern int writeToGhostElemsBeginEvent;
  extern int writeToGhostElemsEndEvent;
  extern int buildDaEvent;
  extern int setMatValuesEvent;
  extern int buildNlistEvent;
  extern int buildNlistCommEvent;

  extern int buildDAstage1Event;
  extern int buildDAstage2Event;
  extern int buildDAstage3Event;
  extern int buildDAstage4Event;
  extern int buildDAstage5Event;
  extern int buildDAstage6Event;
  extern int buildDAstage7Event;
  extern int buildDAstage8Event;
  extern int buildDAstage9Event;
}

#define PROF_WRITE_GHOST_ELEMS_BEGIN_BEGIN PetscLogEventBegin(writeToGhostElemsBeginEvent,0,0,0,0);
#define PROF_WRITE_GHOST_ELEMS_BEGIN_END \
  PetscLogEventEnd(writeToGhostElemsBeginEvent,0,0,0,0); \
return 0;

#define PROF_WRITE_GHOST_ELEMS_END_BEGIN PetscLogEventBegin(writeToGhostElemsEndEvent,0,0,0,0);
#define PROF_WRITE_GHOST_ELEMS_END_END \
  PetscLogEventEnd(writeToGhostElemsEndEvent,0,0,0,0); \
return 0;

#define PROF_WRITE_GHOST_NODES_BEGIN_BEGIN PetscLogEventBegin(writeToGhostNodesBeginEvent,0,0,0,0);
#define PROF_WRITE_GHOST_NODES_BEGIN_END \
  PetscLogEventEnd(writeToGhostNodesBeginEvent,0,0,0,0); \
return 0;

#define PROF_WRITE_GHOST_NODES_END_BEGIN PetscLogEventBegin(writeToGhostNodesEndEvent,0,0,0,0);
#define PROF_WRITE_GHOST_NODES_END_END \
  PetscLogEventEnd(writeToGhostNodesEndEvent,0,0,0,0); \
return 0;

#define PROF_READ_GHOST_ELEMS_BEGIN_BEGIN PetscLogEventBegin(readFromGhostElemsBeginEvent,0,0,0,0);
#define PROF_READ_GHOST_ELEMS_BEGIN_END \
  PetscLogEventEnd(readFromGhostElemsBeginEvent,0,0,0,0); \
return 0;

#define PROF_READ_GHOST_ELEMS_END_BEGIN PetscLogEventBegin(readFromGhostElemsEndEvent,0,0,0,0);
#define PROF_READ_GHOST_ELEMS_END_END \
  PetscLogEventEnd(readFromGhostElemsEndEvent,0,0,0,0); \
return 0;

#define PROF_READ_GHOST_NODES_BEGIN_BEGIN PetscLogEventBegin(readFromGhostNodesBeginEvent,0,0,0,0);
#define PROF_READ_GHOST_NODES_BEGIN_END \
  PetscLogEventEnd(readFromGhostNodesBeginEvent,0,0,0,0); \
return 0;

#define PROF_READ_GHOST_NODES_END_BEGIN PetscLogEventBegin(readFromGhostNodesEndEvent,0,0,0,0);
#define PROF_READ_GHOST_NODES_END_END \
  PetscLogEventEnd(readFromGhostNodesEndEvent,0,0,0,0); \
return 0;

#define PROF_SET_MAT_VALUES_BEGIN PetscLogEventBegin(setMatValuesEvent,0,0,0,0);
#define PROF_SET_MAT_VALUES_END \
  PetscLogEventEnd(setMatValuesEvent,0,0,0,0); \
return 0;

#define PROF_BUILD_DA_BEGIN PetscLogEventBegin(buildDaEvent,0,0,0,0);
#define PROF_BUILD_DA_END \
  PetscLogEventEnd(buildDaEvent,0,0,0,0); \
return;

#define PROF_BUILD_NLIST_BEGIN PetscLogEventBegin(buildNlistEvent,0,0,0,0);
#define PROF_BUILD_NLIST_END \
  PetscLogEventEnd(buildNlistEvent,0,0,0,0); \
return ;

#define PROF_BUILD_NLIST_COMM_BEGIN PetscLogEventBegin(buildNlistCommEvent,0,0,0,0);
#define PROF_BUILD_NLIST_COMM_END PetscLogEventEnd(buildNlistCommEvent,0,0,0,0); 

#define PROF_BUILD_DA_STAGE1_BEGIN PetscLogEventBegin(buildDAstage1Event,0,0,0,0);
#define PROF_BUILD_DA_STAGE1_END PetscLogEventEnd(buildDAstage1Event,0,0,0,0); 

#define PROF_BUILD_DA_STAGE2_BEGIN PetscLogEventBegin(buildDAstage2Event,0,0,0,0);
#define PROF_BUILD_DA_STAGE2_END PetscLogEventEnd(buildDAstage2Event,0,0,0,0); 

#define PROF_BUILD_DA_STAGE3_BEGIN PetscLogEventBegin(buildDAstage3Event,0,0,0,0);
#define PROF_BUILD_DA_STAGE3_END PetscLogEventEnd(buildDAstage3Event,0,0,0,0); 

#define PROF_BUILD_DA_STAGE4_BEGIN PetscLogEventBegin(buildDAstage4Event,0,0,0,0);
#define PROF_BUILD_DA_STAGE4_END PetscLogEventEnd(buildDAstage4Event,0,0,0,0); 

#define PROF_BUILD_DA_STAGE5_BEGIN PetscLogEventBegin(buildDAstage5Event,0,0,0,0);
#define PROF_BUILD_DA_STAGE5_END PetscLogEventEnd(buildDAstage5Event,0,0,0,0); 

#define PROF_BUILD_DA_STAGE6_BEGIN PetscLogEventBegin(buildDAstage6Event,0,0,0,0);
#define PROF_BUILD_DA_STAGE6_END PetscLogEventEnd(buildDAstage6Event,0,0,0,0); 

#define PROF_BUILD_DA_STAGE7_BEGIN PetscLogEventBegin(buildDAstage7Event,0,0,0,0);
#define PROF_BUILD_DA_STAGE7_END PetscLogEventEnd(buildDAstage7Event,0,0,0,0); 

#define PROF_BUILD_DA_STAGE8_BEGIN PetscLogEventBegin(buildDAstage8Event,0,0,0,0);
#define PROF_BUILD_DA_STAGE8_END PetscLogEventEnd(buildDAstage8Event,0,0,0,0); 

#define PROF_BUILD_DA_STAGE9_BEGIN PetscLogEventBegin(buildDAstage9Event,0,0,0,0);
#define PROF_BUILD_DA_STAGE9_END PetscLogEventEnd(buildDAstage9Event,0,0,0,0); 

#else

#define PROF_READ_GHOST_NODES_BEGIN_BEGIN 
#define PROF_READ_GHOST_NODES_BEGIN_END return 0;

#define PROF_READ_GHOST_NODES_END_BEGIN 
#define PROF_READ_GHOST_NODES_END_END return 0;

#define PROF_READ_GHOST_ELEMS_BEGIN_BEGIN 
#define PROF_READ_GHOST_ELEMS_BEGIN_END return 0;

#define PROF_READ_GHOST_ELEMS_END_BEGIN 
#define PROF_READ_GHOST_ELEMS_END_END return 0;

#define PROF_WRITE_GHOST_NODES_BEGIN_BEGIN 
#define PROF_WRITE_GHOST_NODES_BEGIN_END return 0;

#define PROF_WRITE_GHOST_NODES_END_BEGIN 
#define PROF_WRITE_GHOST_NODES_END_END return 0;

#define PROF_WRITE_GHOST_ELEMS_BEGIN_BEGIN 
#define PROF_WRITE_GHOST_ELEMS_BEGIN_END return 0;

#define PROF_WRITE_GHOST_ELEMS_END_BEGIN 
#define PROF_WRITE_GHOST_ELEMS_END_END return 0;

#define PROF_SET_MAT_VALUES_BEGIN 
#define PROF_SET_MAT_VALUES_END return 0;

#define PROF_BUILD_DA_BEGIN 
#define PROF_BUILD_DA_END return;

#define PROF_BUILD_NLIST_BEGIN 
#define PROF_BUILD_NLIST_END return ;

#define PROF_BUILD_NLIST_COMM_BEGIN 
#define PROF_BUILD_NLIST_COMM_END 

#define PROF_BUILD_DA_STAGE1_BEGIN 
#define PROF_BUILD_DA_STAGE1_END 

#define PROF_BUILD_DA_STAGE2_BEGIN 
#define PROF_BUILD_DA_STAGE2_END 

#define PROF_BUILD_DA_STAGE3_BEGIN 
#define PROF_BUILD_DA_STAGE3_END 

#define PROF_BUILD_DA_STAGE4_BEGIN 
#define PROF_BUILD_DA_STAGE4_END 

#define PROF_BUILD_DA_STAGE5_BEGIN 
#define PROF_BUILD_DA_STAGE5_END 

#define PROF_BUILD_DA_STAGE6_BEGIN 
#define PROF_BUILD_DA_STAGE6_END 

#define PROF_BUILD_DA_STAGE7_BEGIN 
#define PROF_BUILD_DA_STAGE7_END 

#define PROF_BUILD_DA_STAGE8_BEGIN 
#define PROF_BUILD_DA_STAGE8_END 

#define PROF_BUILD_DA_STAGE9_BEGIN 
#define PROF_BUILD_DA_STAGE9_END 

#endif

namespace ot {

  //Forward Declaration
  class TreeNode;

  /** 
   * @brief 		Class that manages the octree mesh.
   * @author		Hari Sundar, hsundar@seas.upenn.edu 
   * @author		Rahul S. Sampath, rahulss@seas.upenn.edu 
   * @date 		2007-01-17
   *
   * Since, the input octree represents elements sorted based on the Morton
   * ordering, we need to created a sorted list of nodes. All nodes in the
   * balanced octree, except for hanging nodes and positive boundary nodes, are
   * associated with exactly one element (octant) in octree. Psuedo-elements are
   * created for all the additional nodes, which are marked. This can be done
   * with a single linear pass (\f$\mathcal{O}( n / p )\f$) followed by a
   * parallel sort (\f$\mathcal{O}( n / p \log n / p )\f$)). We get two major
   * benefits from this. Firstly, we store only one ordering for both the
   * elements and the nodes, and the size of the combined array is \f$\Theta ( 2
   * n )\f$ for an input balanced octree of \f$n\f$ elements. Secondly, we are
   * now in a position to compress this data structure, from one utilizing 4
   * <tt>int</tt> to only one <tt>char</tt> per node. On systems with 4
   * <tt>byte</tt> words, this translated to a 16x compression. The algorithm
   * does not add any overheads for in order iterations over either the nodes or
   * the elements.
   *
   * Only one <tt>char</tt> (8 bits) is stored per element/node. 
   * This stores a number of things, including the level of the 
   * element (5 bits), whether it is a Node (N, 1-bit), whether it 
   * is a Boundary node (BN, 1-bit), and whether it is a dependent
   * element(DE, 1-bit). A dependent element is one which
   * has at least 1 node, which is owned by another processor.
   * Whether the element/node is a ghost or not can be determined
   * purely by it's index.
   *
   * The bit structure for this <tt>char</tt> is
   *
   *  N BN DE LEVEL
   *
   * where the highest bit determines whether it is an element or not. 
   *
   * Random access into the elemental/nodal array is however not currently
   * possible, and should not be required. The DA can be extended if required to
   * support random access, by having additional markers, for every \f$k\f$
   * nodes, within the nodal array, which will enable searching in
   * \f$\mathcal{O}( \log n / k + k )\f$ and increase storage by \f$3 n / k\f$
   * <tt>int</tt>. For the search to be efficient, we will need to choose \f$k\f$
   * such that \f$k = \log n\f$.  This will allow us to search within the DA in
   * \f$\mathcal{O}( \log n )\f$, but will increase storage from \f$\Theta ( 2 n
   * )\f$ to \f$\Theta ( 2 n + 12 n / \log n )\f$.
   *
   * In addition we need the neighbourhood data to be prepared for each
   * element/node. In case of and elemental iteration, we need to store the
   * information pertaining to the 7 additional nodes (in 3D) defining this
   * hexahedral element. In the case of a node-based iteration, we will need to
   * store the information about 26 neighbouring nodes, connected to this node
   * via adjacent elements.
   * 
   * The major issue with this is the storage of the indices, 
   * which is 8 times the number of nodes in the octree. In order 
   * to reduce storage costs only the offsets from the current 
   * index are stored.  In addition, these offsets are compressed 
   * by using a variant of the Golomb-Rice encoding. 
   *  
   * In order to overlap computation and communication, it is 
   * preferable that the Finite Element calculations be done in 
   * two stages. The communication is initialized first, and 
   * simultaneously we process the independent elements. We then 
   * wait for the communication to finish before processing the 
   * dependent element. A sample loop will look like, 
   *  
   * @code 
   * // Start communication 
   * readFromGhostsBegin(array, true); 
   * // Process independent elements 
  * for ( init<ot::DA_FLAGS::INDEPENDENT> ;
      curr() < end<ot::DA_FLAGS::INDEPENDENT>(); 
      next<ot::DA_FLAGS::INDEPENDENT>() ) { 
    *     // process current element ...
      * } 
  *  
    * // wait for communication to end.
    * readFromGhostsEnd(); 
  *  
    * // Now process the dependent elements 
    * for ( init<ot::DA_FLAGS::DEPENDENT>; 
        curr() < end<ot::DA_FLAGS::DEPENDENT>(); 
        next<ot::DA_FLAGS::DEPENDENT>() ) { 
      *      // process current element ...
        * } 
  * @endcode 
    */
    class DA {

      protected:

        LoopCounters m_lcLoopInfo;

        //The compressed Octree. Only the levels of the octants are stored.
        unsigned char*          m_ucpOctLevels;
        bool                    m_bCompressLut;
        bool                    m_bComputedLocalToGlobal;
        bool                    m_bComputedLocalToGlobalElems;
        bool                    m_bIamActive;
        unsigned int            m_uiInputSize;
        std::vector<ot::TreeNode>	   m_tnBlocks;
        std::vector<ot::TreeNode>	   m_tnMinAllBlocks;
        std::vector<unsigned int>          m_uiNlist;  
        unsigned int*                      m_uiNlistPtr;
        DendroIntL*                      m_dilpLocalToGlobal;
        DendroIntL*                      m_dilpLocalToGlobalElems;

        // The remainders from the Golomb-Rice encoding of the neighbour
        // look-up-table. 
        std::vector<unsigned char>          m_ucpLutRemainders;
        unsigned char*                      m_ucpLutRemaindersPtr;
        std::vector<unsigned char>          m_ucpSortOrders;
        unsigned char*                      m_ucpSortOrdersPtr;

        // The quotient from the Golomb-Rice encoding of the neighbour
        // look-up-table. 
        std::vector<unsigned short>         m_uspLutQuotients;
        unsigned short*                     m_uspLutQuotientsPtr;

        // The masks for the look-up table ...
        // the size of the mask will be 2*nelem, the first character storing the 
        // compression mask and the second the hanging node mask.
        std::vector<unsigned char>      m_ucpLutMasks;
        unsigned char*                  m_ucpLutMasksPtr;

        std::vector<unsigned char>      m_ucpPreGhostConnectivity;
        unsigned char*                  m_ucpPreGhostConnectivityPtr;

        std::vector<Point>              m_ptsPreGhostOffsets;
        Point*                          m_ptsPreGhostOffsetsPtr;

        std::vector<Point>              PreGhostAnchors;

        // The counter that maintains the current
        //quotient index while looping through the DA.

        unsigned int                    m_uiQuotientCounter;
        unsigned int                    m_uiPreGhostQuotientCnt;

        // The quotients for element begin and independent begin.
        unsigned int                    m_uiElementQuotient;
        unsigned int                    m_uiIndependentElementQuotient;

        // The offset (x,y,z) for the region spanned by the
        //octants on the current processor, including ghost octants.
        Point                           m_ptGhostedOffset;

        // The offset (x,y,z) for the region spanned by the octants on the current
        // processor.
        Point                           m_ptOffset;

        Point                           m_ptIndependentOffset;
        Point                           m_ptCurrentOffset;

        //------------------------------------------------------------------------

        // The number of nodes owned by the current processor.
        unsigned int                    m_uiNodeSize;
        unsigned int                    m_uiBoundaryNodeSize;
        // The number of elements owned by the current processor.
        unsigned int                    m_uiElementSize;

        unsigned int                    m_uiIndependentElementSize;
        // The number of pre-ghosts elements, i.e., ghost elements with morton ids less than
        //  those belonging to this processor.
        unsigned int                    m_uiPreGhostElementSize;

        // The number of pre-ghosts nodes, i.e., ghost nodes with morton ids less than
        //  those belonging to this processor.
        unsigned int                    m_uiPreGhostNodeSize;
        unsigned int                    m_uiPreGhostBoundaryNodeSize;
        // The number of post-ghosts nodes, i.e., ghost nodes with morton ids greater
        //  than those belonging to this processor.
        unsigned int                    m_uiPostGhostNodeSize;    

        unsigned int                    m_uiLocalBufferSize;

        // indices ...
        unsigned int                           m_uiElementBegin;
        unsigned int                           m_uiElementEnd;
        unsigned int                           m_uiPostGhostBegin;

        unsigned int                           m_uiIndependentElementBegin;
        unsigned int                           m_uiIndependentElementEnd;

        unsigned int                           m_uiCurrent;
        //------------------------------------------------------------------------

        // The dimensionality of the underlying space. Only tested in 3D.
        unsigned int                    m_uiDimension;
        // The maximum depth of the octree discretization. Is usually obtained from
        // the octree generation and balancing procedure.
        unsigned int                    m_uiMaxDepth;

        // The MPI Communication context. 
        MPI_Comm                        m_mpiCommAll;
        int                             m_iRankAll;
        int                             m_iNpesAll;
        MPI_Comm                        m_mpiCommActive;
        int                             m_iRankActive;
        int                             m_iNpesActive;
        // Stuff for communication and synchronization ...
        //-------------------------------------------------------------------------
        std::vector<unsigned int>               m_uipScatterMap;
        std::vector<unsigned int>               m_uipSendOffsets;
        std::vector<unsigned int>               m_uipSendProcs;
        std::vector<unsigned int>               m_uipSendCounts;

        std::vector<unsigned int>               m_uipElemScatterMap;
        std::vector<unsigned int>               m_uipElemSendOffsets;
        std::vector<unsigned int>               m_uipElemSendProcs;
        std::vector<unsigned int>               m_uipElemSendCounts;

        std::vector<unsigned int>               m_uipRecvOffsets;
        std::vector<unsigned int>               m_uipRecvProcs;
        std::vector<unsigned int>               m_uipRecvCounts;

        std::vector<unsigned int>               m_uipElemRecvOffsets;
        std::vector<unsigned int>               m_uipElemRecvProcs;
        std::vector<unsigned int>               m_uipElemRecvCounts;

        // tags and contexts ...
        std::vector<updateContext>              m_mpiContexts;
        unsigned int                            m_uiCommTag;

      public:

        /**
          @return the list of blocks owned by the calling processor
          @author Rahul Sampath
          */
        std::vector<ot::TreeNode> getBlocks();

        /**
          @return the number of blocks owned by the calling processor
          @author Rahul Sampath
          */
        unsigned int getNumBlocks();

        /**
          @brief A function to get information about the partitioning of the mesh.
          @return The smallest (in Morton ordering not size) block on each processor
          @author Rahul Sampath
          */
        std::vector<ot::TreeNode> getMinAllBlocks();

        /**
          @return The size of the input octree while calling the DA constructor 
          @author Rahul Sampath 
          */
        unsigned int getInputSize();

        /**
          @name Constructors and destructors
          */
        //@{

        /**
          @author Hari Sundar
          @author Rahul Sampath
          @brief  The constructor that builds the DA from a Sorted, Linear,
          Complete, 2:1 Balanced, Octree.
          @param	comm The MPI communication context which is to be used for setting up the DA.
          @param iAmActive If an address is provided then the status of the calling processor is returned in that address.
          'true' if the mesh is distributed on the calling processor and 'false' otherwise.
          @param blocksPtr If this parameter is NOT NULL, then the input is assumed to be a valid partition of the octree that 
          can be used for the mesh as it is. If this parameter is NULL, then a fresh partition is computed inside the
          constructor itself. Typical users should pass NULL for this parameter.
          @param compressLut Pass 'true' to compress the element-to-node mappings using Goloumb-Rice encoding and 'false' otherwise.
          Note that uncompressing these mappings during the Matvecs will incur a small overhead.
          @param	in   The sorted, complete, linear, 2:1 balanced
          octree that is to be used for constructing the DA.
          @param inputActiveComm An MPI_Comm consisting of only those processors which call the constructor with an non-empty 'in' vector. This is a subset of comm, which can be got by calling the function par::splitComm2way.

          The constructor builds the DA based on the input octree. It constructs all the required neighbour 
          information and also (optionally) compresses the structure to reduce the memory footprint.
          */
        DA(std::vector<ot::TreeNode> &in, MPI_Comm comm, MPI_Comm inputActiveComm,  
            bool compressLut = false, const std::vector<ot::TreeNode> *blocksPtr = NULL,
            bool *iAmActive = NULL);

        /**
          @author Rahul Sampath
          @brief Private constructor used from within DAMG.
          This assumes that the input already included positive boundaries,
          the original octree has been embedded into the larger octree and the
          maxdepth has been set apropriately and the octants without
          hanging anchors have been tagged.
          */
        DA(unsigned int dummy, std::vector<ot::TreeNode> &in, MPI_Comm comm, 
            MPI_Comm inputActiveComm,  bool compressLut = false,
            const std::vector<ot::TreeNode> *blocksPtr = NULL, bool *iAmActive = NULL);

        /**
          @author Hari Sundar
          @brief The destructor for the DA object
          */   
        ~DA();
        //@}

        /**
          @name Information about the DA / domain
          */
        //@{

        /** 
          @brief Get the MPI communicator from the DA.
          @return MPI_Comm
          @author Rahul Sampath
          */
        MPI_Comm getComm();

        /** 
          @brief Get the communicator only containing the active processors.
          @return MPI_Comm
          @author Rahul Sampath
          */
        MPI_Comm getCommActive();

        /**
          @return the total number of processors (including inactive processors)
          @author Rahul Sampath
          */
        int getNpesAll();

        /**
          @return the number of active processors
          @author Rahul Sampath
          */
        int getNpesActive();

        /**
          @return the rank of the calling processor, when both
          active and inactive processors are involved.  
          @author Rahul Sampath
          */
        int getRankAll();

        /**
          @return the rank of the calling processors, when only the active processors are involved
          @author Rahul Sampath
          */
        int getRankActive();

        /** 
          @author Hari Sundar
          @brief Get the offset for the current index. 
          @return The offset. 
          Must not be called by inactive processors

          @see getOffset() 
          @see getGhostedOffset() 
          @see iAmActive()
          */
        Point getCurrentOffset();

        /**
          @author Hari Sundar
          @brief Get the offset for the smallest element on this processor, including ghost elements.
          @return The offset. 

          Must not be called by inactive processors
          @see getOffset() 
          @see getCurrentOffset() 
          @see iAmActive()
          */
        Point getGhostedOffset();

        /**
          @author Hari Sundar
          @brief Get the offset for the smallest element on this processor, not including ghost elements.
          @return The offset. 

          Must not be called by inactive processors
          @see getGhostedOffset() 
          @see getCurrentOffset() 
          @see iAmActive()
          */
        Point getOffset();

        /**
          @author Hari Sundar
          @brief Given an octant specified by a point (its anchor) and
          its level it returns the anchor of the octant
          that immediately follows the given point in the Morton ordering.
          This function assumes that the octree is linear.        
          @param p The anchor of the current octant
          @param d The level of the current octant
          @return the anchor of the next octant 
          */
        Point getNextOffset(Point p, unsigned char d);

        /**
          @author Hari Sundar
          @brief Points to the next anchor. This function is required because we only 
          store the levels of the octants and not their anchors. So the anchors are
          computed on the fly within the loops.
          */
        void incrementCurrentOffset();

        /**
          @author Rahul Sampath
          @brief Points to the anchor of the next pre-ghost octant.
          This function is required because we only 
          store the levels of the octants and not their anchors. So the anchors are
          computed on the fly within the loops.
          */
        void incrementPreGhostOffset();

        /**
          @author Rahul Sampath
          @brief Call this function to check if curr() points to an octant
          touching the domain boundaries from the inside. This function is for real octants only, 
          pseudo-octants can not be tested using this function.
          @param flags The type of boundary is returned in this variable.
          The type is one of the enumerations in BoundaryType2
          @return true if curr() points to an internal boundary octant      
          @see curr()
          @see TreeNode::BoundaryType2
          */ 
        bool isBoundaryOctant(unsigned char *flags=NULL);

        /**
          @author Rahul Sampath
          @return The local to global map computed using the function computeLocalToGlobalMappings()      
          @see computeLocalToGlobalMappings
          @see computedLocalToGlobal
          */
        DendroIntL* getLocalToGlobalMap();

        /**
          @author Rahul Sampath
          @return The local to global map computed using the function computeLocalToGlobalElemMappings()      
          @see computeLocalToGlobalElemMappings
          @see computedLocalToGlobalElems
          */
        DendroIntL* getLocalToGlobalElemsMap();

        /**
          @author Hari Sundar
          @brief Returns the total number of Nodes belonging to this processor. 
          This does not include the ghost nodes. This will include the 
          positive boundary nodes if any belong to this processor.
          @return The number of local nodes.
          */
        unsigned int getNodeSize();

        /**
          @author Rahul Sampath
          @brief Returns the total number of positive Boundary Nodes belonging to this processor. 
          This does not include the ghost nodes. 
          @return The number of local (positive) boundary nodes.
          */
        unsigned int getBoundaryNodeSize();

        /**
          @author Rahul Sampath
          @brief Returns the total number of internal Nodes belonging to this processor. 
          This does not include the ghost nodes and positive boundaries . 
          @return The number of local internal nodes.
          */
        unsigned int getInternalNodeSize();

        /**
          @author Hari Sundar
          @brief Returns the total number of elements belonging to this processor. 
          This does not include the ghost elements.  
          @return The number of local elements.
          */
        unsigned int getElementSize();

        /**
          @author Rahul Sampath
          @brief Returns the total number of pre-ghost elements. 
          */
        unsigned int getPreGhostElementSize();

        /**
          @author Rahul Sampath
          @brief Returns the number of INDEPENDENT elements belonging to this processor. 
          @return The number of local elements.
          */
        unsigned int getIndependentSize();

        /**
          @author Hari Sundar
          @brief Returns the total number of Nodes on this processor, 
          including the ghost nodes. This will include the 
          boundary nodes if any belong to this processor.
          @return The number of nodes.
          */
        unsigned int getGhostedNodeSize();

        /**
          @author Hari Sundar
          @brief Returns the total number of elements on this processor, 
          including the ghost elements. 
          @return The number of nodes.
          */
        unsigned int getGhostedElementSize();

        /**
          @author Hari Sundar
          @return the index of the first local element
          */
        unsigned int getIdxElementBegin();

        /**
          @author Rahul Sampath
          @return the index of the last local element
          */
        unsigned int getIdxElementEnd();

        /**
          @author Hari Sundar
          @return the index of the first post ghost element
          */
        unsigned int getIdxPostGhostBegin();

        /**
          @author Hari Sundar
          @brief Returns the maximum depth (level) of the octree from which this DA was created.

          @return The maximum depth of the octree.
          The return value is the maximum depth in the modified octree that includes 'pseudo-octants' for boundary nodes. This octree has
          a maximum depth equal to 1 more than that of the input octree used to construct the finite element mesh. Hence, the value
          returned by this function will be 1 more than the true maximum depth of the input octree.
          */ 
        unsigned int getMaxDepth();

        /**
          @author Hari Sundar
          @brief Returns the dimension of the octree from which this DA was created.
          @return The dimension of the octree.
          */
        unsigned int getDimension();

        //@}

        /** 
          @name Communication functions 
          */
        //@{

        /**
          @author Hari Sundar
          @author Rahul Sampath	 
          @brief Updates the ghost values by obtaining values from the processors which own them.

          @param arr		the local buffer which needs to be updated. This must be obtained with a call to
          vecGetBuffer().
          @param isElemental	specifies whether the current buffer is elemental (true) or nodal (false).
          @param dof		The degrees of freedom for the current vector, default is 1.
          @see ReadFromGhostsEnd()
          Updates the ghost values by obtaining values from the processors which own them.
          ReadFromGhostsEnd()
          must be called before the ghosted values can be used. 
          */
        //Communicating Ghost Nodes
        template <typename T>
          int ReadFromGhostsBegin ( T* arr, unsigned int dof=1);

        /**
          @author Hari Sundar
          @author Rahul Sampath
         * @brief Waits for updates of the ghost values to finish.
        **/ 
        template <typename T>
          int ReadFromGhostsEnd(T* arr);


        /**
          @author Rahul Sampath
         * @brief Send the ghost values to the processors that own them so that these
         values can be added.  
        **/ 
        template <typename T>
          int WriteToGhostsBegin ( T* arr, unsigned int dof=1);

        /**
          @author Rahul Sampath
         * @brief Waits for updates of the ghost values to finish.
        **/ 
        template <typename T>
          int WriteToGhostsEnd(T* arr, unsigned int dof=1);


        /**
          @author Rahul Sampath
          Counterpart of ReadFromGhostsBegin for elemental arrays
          @see ReadFromGhostsBegin()
          */
        template <typename T>
          int ReadFromGhostElemsBegin ( T* arr, unsigned int dof=1);

        /**
          @author Rahul Sampath
          Counterpart of ReadFromGhostsEnd() for elemental arrays
          @see ReadFromGhostsEnd()
          */
        template <typename T>
          int ReadFromGhostElemsEnd(T* arr);

        /**
          @author Rahul Sampath
          Counterpart of WriteToGhostsBegin() for elemental arrays
          @see WriteToGhostsBegin()
          */
        template <typename T>
          int WriteToGhostElemsBegin ( T* arr, unsigned int dof=1);

        /**
          @author Rahul Sampath
          Counterpart of WriteToGhostsEnd for elemental arrays
          @see WriteToGhostsEnd()
          */
        template <typename T>
          int WriteToGhostElemsEnd(T* arr, unsigned int dof=1);

        //@}

        /**
          @name Array access functions 
          */
        //@{

        /**
          @author Hari Sundar
          @brief Returns a PETSc vector of appropriate size of the requested type.
          @param local     the local vector, a PETSc vector that may be used with the PETSc routines.
          @param isElemental true if an elemental vector is desired, 
          false for a nodal vector.
          @param isGhosted true if memory is to be allocated for ghost values.
          @param dof       the degrees of freedom for the vector. The default is 1.
          @return PETSc error code.
          */
        int createVector(Vec &local, bool isElemental, bool isGhosted, unsigned int dof=1);

        /**
          @author Rahul Sampath
          @brief Similar to createVector(), except the vector is only distributed on the active processors.
          @see createVector()
          */
        int createActiveVector(Vec &local, bool isElemental, bool isGhosted, unsigned int dof=1);

        /**
          @author Rahul Sampath
          @brief Returns a PETSc Matrix of appropriate size of the requested type.
          @param M the matrix
          @param mtype the type of matrix
          @param dof the number of degrees of freedom per node.
          */
        int createMatrix(Mat &M, MatType mtype, unsigned int dof=1);

        /**
          @author Rahul Sampath
          @brief Similar to createMatrix, except the matrix is only distributed on the active processors.
          @see createMatrix()
          */
        int createActiveMatrix(Mat &M, MatType mtype, unsigned int dof=1);

        /**
          @author Rahul Sampath
          @brief Computes mappings between the local and global numberings for nodal buffers.
          @see setValuesInMatrix()
          Call this function only if you need to create Matrices using this mesh. This function must be called
          once before calling setValuesInMatrix(). This function should not
          be called more than once for a given mesh.
          */
        int computeLocalToGlobalMappings();

        /**
          @author Rahul Sampath
          @brief Computes mappings between the local and global numberings for elemental buffers.
          This function is probably required only for developers. Typical users will not need this.
          This function should not be called more than once for a given mesh.
          */
        int computeLocalToGlobalElemMappings();

        /**
          @author Rahul Sampath
          @return 'true' if the function computeLocalToGlobalMappings() was called for this mesh.
          @see computeLocalToGlobalMappings()
          */
        bool computedLocalToGlobal();

        /**
          @author Rahul Sampath
          @return 'true' if the function computeLocalToGlobalElemMappings() was called for this mesh.
          @see computeLocalToGlobalElemMappings()
          */
        bool computedLocalToGlobalElems();

        /**
          @author Rahul Sampath
          @brief a wrapper for setting values into the Matrix.
          This internally calls PETSc's MatSetValues() function.
          @param mat The matrix
          @param records The values and their indices
          @param dof the number of degrees of freedom per node
          @param mode Either INSERT_VALUES or ADD_VALUES
          @return an error flag
          Call PETSc's MatAssembly routines to assemble the matrix after setting the values.
          'records' will be cleared inside the function. It would be more efficient to set values in chunks by
          calling this function multiple times with different sets of
          values instead of a single call at the end of the loop. One
          can use the size of 'records' to determine the number of
          such chunks. Calls to this function with the INSERT_VALUES and ADD_VALUES
          options cannot be mixed without intervening calls to PETSc's MatAssembly routines.
          */
        int setValuesInMatrix(Mat mat, std::vector<ot::MatRecord>& records,
            unsigned int dof, InsertMode mode);

        /**
          @author Hari Sundar
          @brief Returns a std. vector of appropriate size of the requested type. 

          @param local the local vector.
          @param isElemental true if an elemental vector is desired, 
          false for a nodal vector.
          @param isGhosted true if memory is to be allocated for ghost 
          values.
          @param dof the degrees of freedom for the vector. The default is 1.
          @return PETSc error code.
          */
        template <typename T>
          int  createVector(std::vector<T> &local, bool isElemental,
              bool isGhosted, unsigned int dof=1);

        /**
          @author Hari Sundar
          @brief Returns a C-array of type PetscScalar from a PETSc Vec for quick local access. 
          @param in The PETSc Vec which needs to be accessed localy. 
          @param out The local C-array which is used to access data 
          localy.
          @param isElemental true if in is an elemental vector, false 
          if it is a nodal vector.
          @param isGhosted true if in contains ghost values.
          @param isReadOnly true if the buffer is required only for 
          reading, should be set to false if writes
          will be performed.
          @param dof the degrees of freedom for the vector. The default is 1.
          @see vecRestoreBuffer()
          Returns a C-array of type PetscScalar from a PETSc Vec for
          quick local access. In addition, this operation is 
          required to use the oda based indexing. vecRestoreBuffer() must be
          called when the buffer is no longer needed.
          If isReadOnly is true, this involves a simple copy of local values from in.
          The ghosts will have junk value in this case. If isReadOnly is false,
          the buffer will be zeroed out first and then the local values from in
          will be copied. The ghosts will have 0 values in this case.
          */
        int vecGetBuffer(Vec in, PetscScalar* &out, bool isElemental,
            bool isGhosted, bool isReadOnly, unsigned int dof=1); 

        /**
          @author Hari Sundar
          @brief Returns a C-array of type T from a distributed std vector for quick local access. 
          @param in The std::vector which needs to be accessed localy. 
          @param out The local C-array which is used to access data 
          localy. 
          @param isElemental true if in is an elemental vector, false 
          if it is a nodal vector.
          @param isGhosted true if in contains ghost values.
          @param isReadOnly true if the buffer is required only for 
          reading, should be set to false if writes
          will be performed.
          @param dof the degrees of freedom for the vector. The default is 1.
          @see vecRestoreBuffer()
          Returns a C-array of type T from a distributed std vector 
          for quick local access. In addition, this operation is 
          required to use the oda based indexing. vecRestoreBuffer() must be 
          called when the buffer is no longer needed.
          T must have a default constructor, which should zero out the object.
          */
        template < typename T >
          int vecGetBuffer(std::vector<T> &in, T* &out, bool isElemental,
              bool isGhosted, bool isReadOnly, unsigned int dof=1); 

        /**
          @author Hari Sundar
          @brief Restores the C-array of type PetscScalar to a PETSc Vec after quick local access. 
          @param in The PETSc Vec which was accessed localy. 
          @param out The local C-array which is used to access data 
          localy. 
          @param isElemental true if in is an elemental vector, false 
          if it is a nodal vector.
          @param isGhosted true if in contains ghost values.
          @param isReadOnly true if the buffer was used only for 
          reading, should be set to false if writes
          will be performed.
          @param dof the degrees of freedom for the vector. The default is 1.
          @see vecGetBuffer()
          Restores the C-array of type PetscScalar to a PETSc Vec after quick local access. 
          */
        int vecRestoreBuffer(Vec in, PetscScalar* out, bool isElemental, 
            bool isGhosted, bool isReadOnly, unsigned int dof=1); 

        /**
          @author Hari Sundar
          @brief Restores the C-array of type T to a distributed std vector after quick local access. 
          @param in The std::vector which was accessed localy. 
          @param out The local C-array which is used to access data 
          localy. 
          @param isElemental true if in is an elemental vector, false 
          if it is a nodal vector.
          @param isGhosted true if in contains ghost values.
          @param isReadOnly true if the buffer was used only for 
          reading, should be set to false if writes
          will be performed.
          @param dof the degrees of freedom for the vector. The default is 1.
          @see vecGetBuffer()
          Restores the C-array of type T to a distributed std vector after quick local access. 
          */
        template < typename T >
          int vecRestoreBuffer(std::vector<T> &in, T* out, bool isElemental,
              bool isGhosted, bool isReadOnly, unsigned int dof=1); 
        //@}

        //----------------------------------

        /**
          @name Element/Node access and iterators
          */
        //@{

        /**
          @author Rahul Sampath
          @author Hari Sundar
          @brief Initializes the internal counters for a new loop. Remember that the DA
          currently only supports elemental loops.
          @param LoopType valid types are All, Local, Independent,
          Dependent, and Ghosted.

          Sample loop through the elements:

          @code
          for ( init<loopType>; curr() < end<loopType>(); next<loopType>() ) { 
        // Do whatever is required ...
        } 
        @endcode

        @see next()
        @see curr()
        @see end()
        */
        template<ot::DA_FLAGS::loopType type>
          void init();

        /**

          @author Hari Sundar
          @author Rahul Sampath
          @brief Returns an index to the begining of the current loop.
          The loop needs to be initialized using a call to
          init().
          @return the index to the begining of the loop.

          Sample loop through the elements:

          @code 
          for ( init<loopType>; curr() < end<loopType>(); next<loopType>() ) { 
        // Do whatever is required ...
        } 
        @endcode 

        @see end()
        @see next()
        @see init()
        */		
        unsigned int curr();

        /**

          @author Rahul Sampath
          @brief Returns an index to the begining of the current loop.
          The loop needs to be initialized using a call to
          init(). This also stores the current position within an
          internal structure. This can be used to re-start a loop using the FROM_STORED loopType.
          @return the index to the begining of the loop.
          @see loopType
          @see curr()            
          */
        unsigned int currWithInfo();

        /**
          @author Rahul Sampath
          @author Hari Sundar
          @brief Returns an index to the end of the current loop.
          The loop needs to be initialized using a call to
          init().

          @return the index to the end of the loop. 

          Sample loop through the elements: 

          @code 
          for ( init<loopType>; curr() < end<loopType>(); next<loopType>() ) { 
        // Do whatever is required ...
        } 
        @endcode 

        @see curr()
        @see init()
        @see next()
        */
        template<ot::DA_FLAGS::loopType type>
          unsigned int end();

        /**
          @author Rahul Sampath
          @author Hari Sundar
          @brief Returns an index to the next element of the current 
          loop. The loop needs to be initialized using a call to
          initializeLoopCounter().

          @return the index to the end of the loop.

          Sample loop through the elements:

          @code
          for ( init<loopType>; curr() < end<loopType>(); next<loopType>() ) { 
        // Do whatever is required ...
        } 
        @endcode 

        @see init()
        @see curr()
        @see end()
        */
        template<ot::DA_FLAGS::loopType type>
          unsigned int next();

        /**
          @author Hari Sundar
          @brief Returns the child number of the current element.
          @return  the child number of the current element
          */
        unsigned char getChildNumber();

        /**
          @author Hari Sundar
          @brief Returns the compressed representation of the octant at index i.
          @param i the index of the octant.
          @return the compressed representation of the octant at index i.
          */
        unsigned char getFlag(unsigned int i);

        /** 
          @author Hari Sundar
          @brief Returns true if the element specified by the index contains a hanging node.
          @param i the index to the element.
          */
        bool isHanging(unsigned int i);

        /**
          @author Hari Sundar
          @brief Returns true if the element/node specified by the index is a Ghost.
          @param i the index to the element/node.
          */
        bool isGhost(unsigned int i);

        /**
          @author Hari Sundar
          @brief Returns true if the element specified by the index corresponds to a node.
          @param i the index to the element.
          */
        bool isNode(unsigned int i);

        /** 
          @author Hari Sundar
          @brief Returns information pertaining to which of the elements 8 nodes are hanging.
          @param i the index to the element.
          @return the bitmask specifying which of the nodes are hanging.

          Returns information pertaining to which of the elements 8 nodes are hanging.	
          */
        unsigned char getHangingNodeIndex(unsigned int i);

        /** 
          @author Hari Sundar
          @brief Returns the type mask for the given element.
          @param i the index to the element.

          Returns the type mask for the given element. The type mask is
          used to identify what kind of an element the element in 
          question is. Information can be obtained from the bits of the 
          mask. The information is stored as NHDL, where N is 1 bit to 
          identify a node, H is 1 bit to identify a hanging element, D 
          is one bit to identify a dependent element, and the remaining 
          5 bits are used to detect the level of the octant. 
          */
        unsigned char getTypeMask(unsigned int i);

        /** 
          @author Hari Sundar
          @brief Returns the level of the octant specified by the index.
          @param i the index to the element/node.
          @return the level of the octant.
          The return value is the level of the octant in the modified octree that includes 'pseudo-octants' for boundary nodes. This octree has
          a maximum depth equal to 1 more than that of the input octree used to construct the finite element mesh. Hence, the value
          returned by this function will be 1 more than the true level of the octant in the input octree.
          */
        unsigned char getLevel(unsigned int i);

        /**
          @author Hari Sundar
          @author Rahul Sampath
          @brief Returns the indices to the nodes of the current 
          element.
          @param nodes   Indices into the nodes of the given element. Should be
          allocated by the user prior to calling.
          @return Error code.
          */ 
        int getNodeIndices(unsigned int* nodes);
        //@}

        /**
          @author Rahul Sampath
          @return the total number of octants (local, ghosts and FOREIGN) stored on the calling processor.
          */
        unsigned int getLocalBufferSize(); 


        /**
          @author Hari Sundar
          @brief Call this function, if a call to getNodeIndices() 
          is skipped within the loop and if the element-to-node mappings are compressed.
          @see getNodeIndices()
          */
        void updateQuotientCounter();

        /**
          @author Rahul Sampath
          @return true if the element-to-node mappings were compressed using Goloumb-Rice encoding
          */
        bool isLUTcompressed();

        /**
          @author Rahul Sampath
          @return true if the calling processor is active
          */
        bool iAmActive();

      protected:

        /**
          @author Rahul Sampath
          @author Hari Sundar
          @param in The octree, including pesudo-octants for 
          positive boundaries and ghosts (both pre and post).
          @brief Builds the element-to-node mappings and (optionally)
          compresses it using Goloumb-Rice encoding.
          */
        void buildNodeList(std::vector<ot::TreeNode> &in);

        /**
          @author Rahul Sampath
          @brief This function is called from within the constructor. 
          It includes all components of the constructor. 
          */
        void DA_FactoryPart0(std::vector<ot::TreeNode>& in, MPI_Comm comm,
            MPI_Comm activeInputComm, bool compressLut, bool* iAmActive);

        void DA_FactoryPart1(std::vector<ot::TreeNode>& in); 

        void DA_FactoryPart2(std::vector<ot::TreeNode>& in); 

        void DA_FactoryPart3(std::vector<ot::TreeNode>& in, MPI_Comm comm, bool compressLut, 
            const std::vector<ot::TreeNode>* blocksPtr, bool* iAmActive);

    };//end DA class definition

}// end namespace

#include "oda.txx"

#endif


