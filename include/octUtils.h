
/**
 * @file octUtils.h
 * @author		Hari Sundar, hsundar@gmail.com
 * @author		Rahul S. Sampath, rahul.sampath@gmail.com 
 * @brief A list of non-member functions for the TreeNode class.
**/ 

#ifndef _OCTUTILS_H_
#define _OCTUTILS_H_

#include "mpi.h"
#include <vector>

#ifdef PETSC_USE_LOG

#include "petscsys.h"

namespace ot {
  extern int markHangingEvent;
  extern int addBdyEvent;
  extern int FLNstage1Event;
  extern int FLNstage2Event;
  extern int FLNstage3Event;
  extern int FLNstage4Event;
  extern int FLNstage5Event;
  extern int FLNstage6Event;
  extern int FLNstage7Event;
  extern int FLNstage8Event;
  extern int FLNstage9Event;
  extern int FLNstage10Event;
  extern int FLNstage11Event;

  extern int mergeRecvKeysBalEvent; 
  extern int prepBalWlistEvent;
  extern int prepBalComm1MssgEvent;
  extern int prepBalComm2MssgEvent;
  extern int finalBalMergeEvent;
  extern int mergeComboBalEvent;
  extern int pickNhBlocksEvent;
  extern int balanceEvent;
  extern int balSubtreeEvent;
  extern int completeSubtreeEvent;
  extern int coarsenEvent;
  extern int coarsenSeqEvent;
  extern int simpleCoarsenEvent;
  extern int mergeOctreesEvent;
  extern int rg2oEvent;
  extern int p2oEvent;
  extern int p2oSeqEvent;
  extern int p2oLocalEvent;
  extern int n2oEvent;
  extern int n2oSeqEvent;
  extern int completeRegionEvent;
  extern int blockPart1Event;
  extern int blockPart2Event;
  extern int blockPart3Event;
  extern int conBalEvent;
  extern int rippleBalEvent;
  extern int ptrRippleBalEvent;
  extern int parRippleType3Event;
  extern int parRippleType2Event;
  extern int parRippleType1Event;
  extern int comboRippleEvent;
  extern int pickBndEvent;
  extern int balCommEvent;
  extern int balScatterEvent;
  extern int balSplitCommEvent;
  extern int balBpart1Event;
  extern int balBpart2Event;
}

#define PROF_FLN_STAGE1_BEGIN PetscLogEventBegin(FLNstage1Event,0,0,0,0);
#define PROF_FLN_STAGE1_END PetscLogEventEnd(FLNstage1Event,0,0,0,0); 

#define PROF_FLN_STAGE2_BEGIN PetscLogEventBegin(FLNstage2Event,0,0,0,0);
#define PROF_FLN_STAGE2_END PetscLogEventEnd(FLNstage2Event,0,0,0,0); 

#define PROF_FLN_STAGE3_BEGIN PetscLogEventBegin(FLNstage3Event,0,0,0,0);
#define PROF_FLN_STAGE3_END PetscLogEventEnd(FLNstage3Event,0,0,0,0); 

#define PROF_FLN_STAGE4_BEGIN PetscLogEventBegin(FLNstage4Event,0,0,0,0);
#define PROF_FLN_STAGE4_END PetscLogEventEnd(FLNstage4Event,0,0,0,0); 

#define PROF_FLN_STAGE5_BEGIN PetscLogEventBegin(FLNstage5Event,0,0,0,0);
#define PROF_FLN_STAGE5_END PetscLogEventEnd(FLNstage5Event,0,0,0,0); 

#define PROF_FLN_STAGE6_BEGIN PetscLogEventBegin(FLNstage6Event,0,0,0,0);
#define PROF_FLN_STAGE6_END PetscLogEventEnd(FLNstage6Event,0,0,0,0); 

#define PROF_FLN_STAGE7_BEGIN PetscLogEventBegin(FLNstage7Event,0,0,0,0);
#define PROF_FLN_STAGE7_END PetscLogEventEnd(FLNstage7Event,0,0,0,0); 

#define PROF_FLN_STAGE8_BEGIN PetscLogEventBegin(FLNstage8Event,0,0,0,0);
#define PROF_FLN_STAGE8_END PetscLogEventEnd(FLNstage8Event,0,0,0,0); 

#define PROF_FLN_STAGE9_BEGIN PetscLogEventBegin(FLNstage9Event,0,0,0,0);
#define PROF_FLN_STAGE9_END PetscLogEventEnd(FLNstage9Event,0,0,0,0); 

#define PROF_FLN_STAGE10_BEGIN PetscLogEventBegin(FLNstage10Event,0,0,0,0);
#define PROF_FLN_STAGE10_END PetscLogEventEnd(FLNstage10Event,0,0,0,0); 

#define PROF_FLN_STAGE11_BEGIN PetscLogEventBegin(FLNstage11Event,0,0,0,0);
#define PROF_FLN_STAGE11_END PetscLogEventEnd(FLNstage11Event,0,0,0,0); 

#define PROF_MARK_HANGING_BEGIN PetscLogEventBegin(markHangingEvent,0,0,0,0);
#define PROF_MARK_HANGING_END \
  PetscLogEventEnd(markHangingEvent,0,0,0,0); \
return ;

#define PROF_ADD_BDY_BEGIN PetscLogEventBegin(addBdyEvent,0,0,0,0);
#define PROF_ADD_BDY_END \
  PetscLogEventEnd(addBdyEvent,0,0,0,0); \
return ;

#define PROF_BAL_COMM_BEGIN PetscLogEventBegin(balCommEvent,0,0,0,0);
#define PROF_BAL_COMM_END PetscLogEventEnd(balCommEvent,0,0,0,0);

#define PROF_BAL_SCATTER_BEGIN PetscLogEventBegin(balScatterEvent,0,0,0,0);
#define PROF_BAL_SCATTER_END PetscLogEventEnd(balScatterEvent,0,0,0,0);

#define PROF_BAL_SPLIT_COMM_BEGIN PetscLogEventBegin(balSplitCommEvent,0,0,0,0);
#define PROF_BAL_SPLIT_COMM_END PetscLogEventEnd(balSplitCommEvent,0,0,0,0);

#define PROF_BAL_BPART1_BEGIN PetscLogEventBegin(balBpart1Event,0,0,0,0);
#define PROF_BAL_BPART1_END PetscLogEventEnd(balBpart1Event,0,0,0,0);

#define PROF_BAL_BPART2_BEGIN PetscLogEventBegin(balBpart2Event,0,0,0,0);
#define PROF_BAL_BPART2_END PetscLogEventEnd(balBpart2Event,0,0,0,0);

#define PROF_MERGE_RECV_KEYS_BAL_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(mergeRecvKeysBalEvent,0,0,0,0);
#define PROF_MERGE_RECV_KEYS_BAL_END \
  PetscLogEventEnd(mergeRecvKeysBalEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PREP_BAL_WLIST_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(prepBalWlistEvent,0,0,0,0);
#define PROF_PREP_BAL_WLIST_END \
  PetscLogEventEnd(prepBalWlistEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PREP_BAL_COMM1_MSSG_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(prepBalComm1MssgEvent,0,0,0,0);
#define PROF_PREP_BAL_COMM1_MSSG_END \
  PetscLogEventEnd(prepBalComm1MssgEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PREP_BAL_COMM2_MSSG_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(prepBalComm2MssgEvent,0,0,0,0);
#define PROF_PREP_BAL_COMM2_MSSG_END \
  PetscLogEventEnd(prepBalComm2MssgEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_FINAL_MERGE_IN_BAL_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(finalBalMergeEvent,0,0,0,0);
#define PROF_FINAL_MERGE_IN_BAL_END \
  PetscLogEventEnd(finalBalMergeEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_MERGE_COMBO_BAL_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(mergeComboBalEvent,0,0,0,0);
#define PROF_MERGE_COMBO_BAL_END \
  PetscLogEventEnd(mergeComboBalEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PICK_NH_BLOCKS_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(pickNhBlocksEvent,0,0,0,0);
#define PROF_PICK_NH_BLOCKS_END \
  PetscLogEventEnd(pickNhBlocksEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_COMBO_RIPPLE_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(comboRippleEvent,0,0,0,0);
#define PROF_COMBO_RIPPLE_END \
  PetscLogEventEnd(comboRippleEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_SIMPLE_COARSE_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(simpleCoarsenEvent,0,0,0,0);
#define PROF_SIMPLE_COARSE_END         \
  PetscLogEventEnd(simpleCoarsenEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_COARSE_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(coarsenEvent,0,0,0,0);
#define PROF_COARSE_END         \
  PetscLogEventEnd(coarsenEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_COARSE_SEQ_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(coarsenSeqEvent,0,0,0,0);
#define PROF_COARSE_SEQ_END         \
  PetscLogEventEnd(coarsenSeqEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_CON_BAL_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(conBalEvent,0,0,0,0);
#define PROF_CON_BAL_END         \
  PetscLogEventEnd(conBalEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_RIPPLE_BAL_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(rippleBalEvent,0,0,0,0);
#define PROF_RIPPLE_BAL_END         \
  PetscLogEventEnd(rippleBalEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PTR_RIPPLE_BAL_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(ptrRippleBalEvent,0,0,0,0);
#define PROF_PTR_RIPPLE_BAL_END         \
  PetscLogEventEnd(ptrRippleBalEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PICK_BND_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(pickBndEvent,0,0,0,0);
#define PROF_PICK_BND_END         \
  PetscLogEventEnd(pickBndEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_BAL_SUBTREE_BEGIN 	\
  PetscFunctionBegin; \
PetscLogEventBegin(balSubtreeEvent,0,0,0,0);
#define PROF_BAL_SUBTREE_END         \
  PetscLogEventEnd(balSubtreeEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_COMPLETE_SUBTREE_BEGIN 	\
  PetscFunctionBegin; \
PetscLogEventBegin(completeSubtreeEvent,0,0,0,0);
#define PROF_COMPLETE_SUBTREE_END         \
  PetscLogEventEnd(completeSubtreeEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_BAL_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(balanceEvent,0,0,0,0);
#define PROF_BAL_END \
  PetscLogEventEnd(balanceEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PAR_RIPPLE_TYPE3_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(parRippleType3Event,0,0,0,0);
#define PROF_PAR_RIPPLE_TYPE3_END \
  PetscLogEventEnd(parRippleType3Event,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PAR_RIPPLE_TYPE2_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(parRippleType2Event,0,0,0,0);
#define PROF_PAR_RIPPLE_TYPE2_END \
  PetscLogEventEnd(parRippleType2Event,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PAR_RIPPLE_TYPE1_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(parRippleType1Event,0,0,0,0);
#define PROF_PAR_RIPPLE_TYPE1_END \
  PetscLogEventEnd(parRippleType1Event,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_MERGE_OCTREES_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(mergeOctreesEvent,0,0,0,0);
#define PROF_MERGE_OCTREES_END \
  PetscLogEventEnd(mergeOctreesEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_RG2O_BEGIN	\
  PetscFunctionBegin; \
PetscLogEventBegin(rg2oEvent,0,0,0,0);
#define PROF_RG2O_END \
  PetscLogEventEnd(rg2oEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_P2O_BEGIN 	\
  PetscFunctionBegin; \
PetscLogEventBegin(p2oEvent,0,0,0,0);
#define PROF_P2O_END         \
  PetscLogEventEnd(p2oEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_P2O_SEQ_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(p2oSeqEvent,0,0,0,0);
#define PROF_P2O_SEQ_END \
  PetscLogEventEnd(p2oSeqEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_P2O_LOCAL_BEGIN 	\
  PetscFunctionBegin; \
PetscLogEventBegin(p2oLocalEvent,0,0,0,0);
#define PROF_P2O_LOCAL_END         \
  PetscLogEventEnd(p2oLocalEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_N2O_BEGIN 	\
  PetscFunctionBegin; \
PetscLogEventBegin(n2oEvent,0,0,0,0);
#define PROF_N2O_END         \
  PetscLogEventEnd(n2oEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_N2O_SEQ_BEGIN 	\
  PetscFunctionBegin; \
PetscLogEventBegin(n2oSeqEvent,0,0,0,0);
#define PROF_N2O_SEQ_END         \
  PetscLogEventEnd(n2oSeqEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_COMPLETE_REGION_BEGIN 	\
  PetscFunctionBegin; \
PetscLogEventBegin(completeRegionEvent,0,0,0,0);
#define PROF_COMPLETE_REGION_END         \
  PetscLogEventEnd(completeRegionEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_BLKPART1_BEGIN 	\
  PetscFunctionBegin; \
PetscLogEventBegin(blockPart1Event,0,0,0,0);
#define PROF_BLKPART1_END         \
  PetscLogEventEnd(blockPart1Event,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_BLKPART2_BEGIN 	\
  PetscFunctionBegin; \
PetscLogEventBegin(blockPart2Event,0,0,0,0);
#define PROF_BLKPART2_END         \
  PetscLogEventEnd(blockPart2Event,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_BLKPART3_BEGIN 	\
  PetscFunctionBegin; \
PetscLogEventBegin(blockPart3Event,0,0,0,0);
#define PROF_BLKPART3_END         \
  PetscLogEventEnd(blockPart3Event,0,0,0,0); \
PetscFunctionReturn(0);

#else

#define PROF_BAL_COMM_BEGIN 
#define PROF_BAL_SCATTER_BEGIN 
#define PROF_BAL_SPLIT_COMM_BEGIN 
#define PROF_BAL_BPART1_BEGIN 
#define PROF_BAL_BPART2_BEGIN 
#define PROF_FINAL_MERGE_IN_BAL_BEGIN
#define PROF_MERGE_COMBO_BAL_BEGIN
#define PROF_PICK_NH_BLOCKS_BEGIN
#define PROF_COARSE_BEGIN
#define PROF_COARSE_SEQ_BEGIN
#define PROF_SIMPLE_COARSE_BEGIN
#define PROF_BAL_BEGIN
#define PROF_BAL_SUBTREE_BEGIN
#define PROF_COMPLETE_SUBTREE_BEGIN
#define PROF_MERGE_OCTREES_BEGIN
#define PROF_RG2O_BEGIN
#define PROF_P2O_BEGIN
#define PROF_P2O_SEQ_BEGIN
#define PROF_P2O_LOCAL_BEGIN
#define PROF_N2O_BEGIN
#define PROF_N2O_SEQ_BEGIN
#define PROF_COMPLETE_REGION_BEGIN
#define PROF_BLKPART1_BEGIN
#define PROF_BLKPART2_BEGIN
#define PROF_BLKPART3_BEGIN
#define PROF_COMBO_RIPPLE_BEGIN
#define PROF_CON_BAL_BEGIN
#define PROF_RIPPLE_BAL_BEGIN
#define PROF_PTR_RIPPLE_BAL_BEGIN
#define PROF_PAR_RIPPLE_TYPE3_BEGIN
#define PROF_PAR_RIPPLE_TYPE2_BEGIN
#define PROF_PAR_RIPPLE_TYPE1_BEGIN
#define PROF_PICK_BND_BEGIN
#define PROF_PREP_BAL_COMM1_MSSG_BEGIN 
#define PROF_PREP_BAL_WLIST_BEGIN 
#define PROF_PREP_BAL_COMM2_MSSG_BEGIN 
#define PROF_MERGE_RECV_KEYS_BAL_BEGIN 

#define PROF_BAL_BPART2_END 
#define PROF_BAL_COMM_END 
#define PROF_BAL_SCATTER_END 
#define PROF_BAL_SPLIT_COMM_END 
#define PROF_BAL_BPART1_END 
#define PROF_MERGE_RECV_KEYS_BAL_END return 1;
#define PROF_PREP_BAL_COMM2_MSSG_END return 1;
#define PROF_PREP_BAL_WLIST_END return 1;
#define PROF_PREP_BAL_COMM1_MSSG_END return 1;
#define PROF_FINAL_MERGE_IN_BAL_END return 1;
#define PROF_MERGE_COMBO_BAL_END return 1;
#define PROF_PICK_NH_BLOCKS_END return 1;
#define PROF_COARSE_END return 1; 
#define PROF_COARSE_SEQ_END return 1; 
#define PROF_SIMPLE_COARSE_END return 1; 
#define PROF_BAL_END return 1; 
#define PROF_BAL_SUBTREE_END return 1; 
#define PROF_COMPLETE_SUBTREE_END return 1; 
#define PROF_MERGE_OCTREES_END return 1;
#define PROF_RG2O_END return 1; 
#define PROF_P2O_END return 1; 
#define PROF_P2O_SEQ_END return 1; 
#define PROF_P2O_LOCAL_END return 1; 
#define PROF_N2O_END return 1;
#define PROF_N2O_SEQ_END return 1;
#define PROF_COMPLETE_REGION_END return 1;
#define PROF_BLKPART1_END return 1;
#define PROF_BLKPART2_END return 1;
#define PROF_BLKPART3_END return 1;
#define PROF_COMBO_RIPPLE_END return 1;
#define PROF_CON_BAL_END return 1; 
#define PROF_RIPPLE_BAL_END return 1; 
#define PROF_PTR_RIPPLE_BAL_END return 1; 
#define PROF_PAR_RIPPLE_TYPE3_END return 1; 
#define PROF_PAR_RIPPLE_TYPE2_END return 1; 
#define PROF_PAR_RIPPLE_TYPE1_END return 1; 
#define PROF_PICK_BND_END return 1;

#define PROF_FLN_STAGE1_BEGIN
#define PROF_FLN_STAGE2_BEGIN
#define PROF_FLN_STAGE3_BEGIN
#define PROF_FLN_STAGE4_BEGIN
#define PROF_FLN_STAGE5_BEGIN
#define PROF_FLN_STAGE6_BEGIN
#define PROF_FLN_STAGE7_BEGIN
#define PROF_FLN_STAGE8_BEGIN 
#define PROF_FLN_STAGE9_BEGIN 
#define PROF_FLN_STAGE10_BEGIN 
#define PROF_FLN_STAGE11_BEGIN 

#define PROF_FLN_STAGE1_END
#define PROF_FLN_STAGE2_END
#define PROF_FLN_STAGE3_END
#define PROF_FLN_STAGE4_END
#define PROF_FLN_STAGE5_END
#define PROF_FLN_STAGE6_END
#define PROF_FLN_STAGE7_END
#define PROF_FLN_STAGE8_END 
#define PROF_FLN_STAGE9_END 
#define PROF_FLN_STAGE10_END 
#define PROF_FLN_STAGE11_END 

#define PROF_MARK_HANGING_BEGIN 
#define PROF_MARK_HANGING_END return ;

#define PROF_ADD_BDY_BEGIN
#define PROF_ADD_BDY_END return;

#endif

/**
  @brief The intra-processor balancing is done in two stages: a search-free intra-block balancing, followed by the ripple algorithm for inter-block balancing. This combined algorithm is actually implemented recursively and this factor controls the level of recursion.
  */
#define _COMBO_RIPPLE_FACTOR_ 1000000

/**
  @namespace ot
  @author Rahul Sampath
  @author Hari Sundar
  @brief A collection of Octree specific functions. 
  */
namespace ot {

  class TreeNode;

  struct TreeNodePointer;

  void appendOctantsAtLevel(const ot::TreeNodePointer & ptrOctree, 
      std::vector<ot::TreeNode> & wList, unsigned int lev);

  void findOctantOrFinestAncestor(ot::TreeNodePointer & octree,
      const ot::TreeNode & key, ot::TreeNodePointer* & result);

  void addOctantToTreeNodePointer(ot::TreeNodePointer & ptrOct, const ot::TreeNode & octant);

  /**
    @author Rahul Sampath
    @param linOct linear octree
    @param ptrOct linked list octree
    @brief Converts an octree from its linear representation to a linked-list representation
    */
  void convertLinearToPointer(const std::vector<ot::TreeNode> & linOct, ot::TreeNodePointer & ptrOct);

  /**
    @author Rahul Sampath
    @param linOct linear octree
    @param ptrOct linked list octree
    @brief Converts an octree from its linked-list representation to a linear representation
    */
  void convertPointerToLinear(std::vector<ot::TreeNode> & linOct, const ot::TreeNodePointer & ptrOct);

  void deleteTreeNodePointer(ot::TreeNodePointer & ptrOct);

  /**
    @author Hari Sundar
    @author Rahul Sampath
    @param in The input octree. This will be modified so
    that the original octree is embedded into a larger octree.
    @param bdy A list to store the pseudo-octants for the positive boundary 
    @brief generates psuedonodes for the positive boundary
    and inserts them into bdy
    */
  void addBoundaryNodesType1( std::vector<ot::TreeNode> &in,
      std::vector<ot::TreeNode> &bdy,
      unsigned int dim, unsigned int maxDepth);

  /**
    @author Rahul Sampath
    @param in The octree including pseudo-octants for
    the positive boundaries
    @brief Identifies hanging nodes. Uses two all2allv communications.
    No overlap of comm and comp.
    */
  void flagNodesType1(std::vector<ot::TreeNode> & in, MPI_Comm comm);

  /**
    @author Rahul Sampath
    @param in the Octree including pseudo-octants for positive boundaries
    @brief Identifies hanging nodes. It 
    just uses 1 communication by sending apriori results. Overlaps comm and comp.
    */
  void flagNodesType2(std::vector<ot::TreeNode> & in, MPI_Comm comm);

  /**
    @author Rahul Sampath
    @param in the Octree including pseudo-octants for positive boundaries
    @brief Identifies hanging nodes. 2 step communication, overlapping comm and comp. 
    */
  void flagNodesType3(std::vector<ot::TreeNode> & in, MPI_Comm comm);


  void addBoundaryNodesType2( std::vector<ot::TreeNode>& in,
      std::vector<ot::TreeNode> &bdy, unsigned int dim, unsigned int maxDepth);

  void discardExtraBoundaryOctants(std::vector<ot::TreeNode>& in,
      unsigned int dim, unsigned int maxDepth);

  void markBoundaryNodesAtAllLevels(std::vector<ot::TreeNode>& finestOctree, unsigned int nlevels,
      std::vector<ot::TreeNode>* coarserOctrees, unsigned int maxDepth);

  void markHangingNodesAtAllLevels(std::vector<ot::TreeNode>& finestOctree, unsigned int nlevels, 
      std::vector<ot::TreeNode>* coarserOctrees, MPI_Comm* activeComms,
      unsigned int dim, unsigned int maxDepth);

  /**
    @brief Makes the input linear. Removes duplicates and ancestors
    @param list the input vector (must be sorted)
    @param skipLast Pass 'true' if you do not wish to include the last element in the output and 'false' otherwise
    @author Rahul Sampath
    */
  int lineariseList(std::vector<ot::TreeNode> & list, bool skipLast = false);

  /**
    @brief Makes the input linear. Removes duplicates and ancestors
    @param list the input vector (must be sorted)  
    @author Rahul Sampath
    */
  int lineariseList(std::vector<ot::TreeNode> & list, MPI_Comm comm);

  /**
    @brief A comparator that uses the weights of the octants instead of the Morton ordering.	
    @return a < b based on their respective weigths.
    */
  bool lessThanUsingWts ( TreeNode  const & a,  TreeNode  const & b);

  /**
    @author Rahul Sampath
    @brief The region between min(first,second) and max(first,second) is appended to the output vector.
    Both ends could be inclusive, dependng on the options. The new elements are sorted, unique and linear.
    @param includeMin include min(first,second)
    @param includeMax include max(first,second)
    */
  int appendCompleteRegion(TreeNode first, TreeNode second, std::vector<ot::TreeNode>& out,
      bool includeMin, bool includeMax);

  /**
    @brief checks if the dim and maxdepths are the same.
    @return true if first and second are comparable
    */
  bool areComparable(TreeNode first, TreeNode second);

  /**
    @brief criteria for picking blocks inside blockPart
    @return true if a < b according to the criteria
    */
  bool bPartComparator(TreeNode a, TreeNode b) ;

  /**
    @author Rahul Sampath
    @return the nearest common ancestor of first and second.
    This works even if both the inputs are equal or if one 
    of them is the ancestor of the other.
*/
  TreeNode getNCA(TreeNode first, TreeNode second);

  /**
    @author Rahul Sampath
    @author Hari Sundar
    @brief A function to construct a complete, sorted, linear octree from a set of points.
    @param points the input points, (x,y,z) must be >=0.0 and < (gLens[0], gLens[1], gLens[2]) respectively
    @param gLens the global dimensions
    @param nodes the output octree
    @param dim the dimension of the tree (1 for binary trees, 2 for quadtrees and 3 for octrees)
    @param maxDepth the maximum depth of the octree must be <= _MAX_LEVEL_
    @param maxNumPtsPerOctant the maximum number of points per octant.
    @return an error flag
    @see _MAX_LEVEL_
    Points are cleared inside the function.
    */
  int points2Octree(std::vector<double>& points, double * gLens, std::vector<TreeNode> & nodes,
      unsigned int dim, unsigned int maxDepth, unsigned int maxNumPtsPerOctant, MPI_Comm comm ) ;

  /**
    @author Rahul Sampath
    @brief A function to construct a complete, sorted, linear octree from a
    regular grid by coarsening the regular grid elements based on some threshold.
    @param elementValues a value at each element ordered in the X,Y,Z order, i.e. X grows first
    length of elementValues must be equal to (nx*ny*nz)
    @param nx,ny,nz number of elements in each direction on the calling processor
    @param N number of elements in each direction in the entire regular grid. N must be a power of 2.
    The total number of elements in the regular grid will be N^3.
    @param xs,ys,zs the global index of the first coordinate on this processor 
    @param maxDepth the maximum depth of the octree must be <= _MAX_LEVEL_
    @see _MAX_LEVEL_
    */
  int regularGrid2Octree(const std::vector<double>& elementValues,
      unsigned int N, unsigned int nx, unsigned int ny, unsigned int nz,
      unsigned int xs, unsigned int ys, unsigned int zs, std::vector<TreeNode>& linOct,
      unsigned int dim, unsigned int maxDepth, double threshold, MPI_Comm comm);

  /**
    @author Rahul Sampath
    @brief Merges the 2 input octrees and linearizes the result
    @param inOct1, inOct2 input octrees
    @param outOct output octree
    */
  int mergeOctrees(std::vector<TreeNode>& inOct1, std::vector<TreeNode>& inOct2,
      std::vector<TreeNode>& outOct, MPI_Comm comm);

  /**
    @author Rahul Sampath
    @brief Sequential version of points2Octree
    @see points2Octree
    */
  int points2OctreeSeq(std::vector<double>& pts, double * gLens, std::vector<TreeNode> & nodes,
      unsigned int dim, unsigned int maxDepth, unsigned int maxNumPts);

  /**
    @author Rahul Sampath
    @author Hari Sundar
    @brief Sequential top-down loop inside points2Octree
    @see points2Octree
    */
  int p2oLocal(std::vector<TreeNode> & nodes, std::vector<TreeNode>& leaves,
      unsigned int maxNumPts, unsigned int dim, unsigned int maxDepth);

  /**
    @author Rahul Sampath
    @author Hari Sundar
    @brief Parallel 2:1 hybrid balance function. It uses a combination of search free and search based approaches to balancing. It does not use parallel searches. Instead, it uses a-priori communication and the concept of insulation layers to de-couple the problem of balancing.
    @param in a sorted, complete, linear octree
    @param out a sorted, complete, linear, balanced octree
    @param dim the dimension of the tree (1 for binary trees, 2 for quadtrees and 3 for octrees)
    @param maxDepth the maximum depth of the octree must be <= _MAX_LEVEL_   
    @param incCorner 'true' to balance across corners as well and 'false' otherwise.
    in is cleared inside the function
    @see _MAX_LEVEL_
    */
  int balanceOctree(std::vector<TreeNode > &in, std::vector<TreeNode > &out,
      unsigned int dim, unsigned int maxDepth, bool incCorner, 
      MPI_Comm comm, MPI_Comm* newCommPtr = NULL, bool* iAmActive = NULL);

  /**
    @author Rahul Sampath
    @brief Merge the input of ConBal and intra-processor ripple bal and 
    pick the inter-processor boundaries from the result
    */
  int mergeComboBalAndPickBoundary(std::vector<ot::TreeNode>& out, 
      std::vector<ot::TreeNode>& allBoundaryLeaves,
      const ot::TreeNode& firstBlock, const ot::TreeNode& lastBlock);

  /**
    @author Rahul Sampath
    */
  int finalMergeInBal(std::vector<ot::TreeNode>& out, std::vector<ot::TreeNode>& allBoundaryLeaves);

  /**
    @author Rahul Sampath
    */
  int prepareBalComm1MessagesType2(const std::vector<ot::TreeNode>& allBoundaryLeaves, 
      const std::vector<ot::TreeNode>& minsAllBlocks, int rank, unsigned int dim,
      unsigned int maxDepth, std::vector<TreeNode>* sendNodes,
      std::vector<unsigned int>* sentToPid, int* sendCnt);

  /**
    @author Rahul Sampath
    */
  int prepareBalComm1MessagesType1(const std::vector<ot::TreeNode>& allBoundaryLeaves, 
      const std::vector<ot::TreeNode>& myNhBlocks, int npes, unsigned int maxDepth, 
      std::vector<TreeNode>* sendNodes, std::vector<unsigned int>* sentToPid, int* sendCnt);

  /**
    @author Rahul Sampath
    */
  int prepareBalComm2Messages(const std::vector<ot::TreeNode>& allBoundaryLeaves,
      const std::vector<ot::TreeNode>& wList,
      const std::vector<std::vector<unsigned int> >& wListRanks,
      std::vector<TreeNode>* sendNodes, std::vector<unsigned int>* sentToPid, int* sendCnt);

  /**
    @author Rahul Sampath
    */
  int prepareWlistInBal(const std::vector<ot::TreeNode>& recvK1, 
      const int* recvCnt, int npes, const ot::TreeNode& myFirstBlock,
      const ot::TreeNode& myLastBlock, std::vector<TreeNode>& wList,
      std::vector<std::vector<unsigned int> >& wListRanks);

  /**
    @author Rahul Sampath
    */
  int mergeRecvKeysInBal(const std::vector<ot::TreeNode>& recvK1, const int* recvOffsets1,
      const std::vector<ot::TreeNode>& recvK2, const int* recvOffsets2, 
      int npes, std::vector<ot::TreeNode>& recvK);

  /**
    @author Rahul Sampath
    */
  int selectNeighboringBlocks(const std::vector<TreeNode>& allBlocks,
      const std::vector<TreeNode>& blocks, const std::vector<unsigned int>& maxBlockBndVec,
      int myRank, std::vector<TreeNode>& myNhBlocks);

  /**
    @author Rahul Sampath
    @brief An implementation of 2:1 balancing using parallel prioritized ripple propagation algorithm
    @see balanceOctree()
    */
  int parallelRippleType1(std::vector<TreeNode> & nodes,
      bool incCorners, bool checkBailOut, bool rePart,
      unsigned int dim, unsigned int maxDepth,
      MPI_Comm comm);

  /**
    @author Rahul Sampath
    @brief An implementation of 2:1 balancing using parallel prioritized ripple propagation algorithm
    @see balanceOctree()
    */
  int parallelRippleType2(std::vector<TreeNode> & nodes,
      bool incCorners, bool checkBailOut, bool rePart,
      unsigned int dim, unsigned int maxDepth,
      MPI_Comm comm);

  /**
    @author Rahul Sampath
    @brief An implementation of 2:1 balancing using parallel prioritized ripple propagation algorithm
    @see balanceOctree()
    */
  int parallelRippleType3(std::vector<TreeNode> & nodes,
      bool incCorners, bool checkBailOut, bool rePart,
      unsigned int dim, unsigned int maxDepth,
      MPI_Comm comm);

  /**
    @author Rahul Sampath
    @brief Constructs the complete, linear, sorted subtree of the block containing the octants in inp
    @param block the root of the subtree
    @param inp a set of decendants of block
    @param out the result
    @param isSorted 'true' if the input is sorted
    @param isUnique 'true' if the input is free from duplicates
    @param dim the dimension of the tree (1 for binary trees, 2 for quadtrees and 3 for octrees)
    @param maxDepth the maximum depth of the octree must be <= _MAX_LEVEL_   
    */
  int completeSubtree(TreeNode block, const std::vector<TreeNode > & inp,
      std::vector<TreeNode > & out, unsigned int dim,
      unsigned int maxDepth, bool isUnique, bool isSorted);

  /**
    @author Rahul Sampath
    @brief Constructs the complete, linear, sorted octree containing the octants in inp  
    @param inp a set of octants
    @param out the result
    @param isSorted 'true' if the input is sorted
    @param isUnique 'true' if the input is free from duplicates
    @param assertNoEmptyProcs 'true' if the input will be non-empty on every
    processor even after removing duplicates across processors. 
    @param dim the dimension of the tree (1 for binary trees, 2 for quadtrees and 3 for octrees)
    @param maxDepth the maximum depth of the octree must be <= _MAX_LEVEL_   
    @see completeSubtree()
    */
  int completeOctree (const std::vector<TreeNode > &in, std::vector<TreeNode > &out,
      unsigned int dim, unsigned int maxDepth, bool isUnique, bool isSorted,
      bool assertNoEmptyProcs, MPI_Comm comm);

  /**
    @author Rahul Sampath
    @brief Coarsens a given octree. Replaces every set of eight siblings with their parent. 
    @param in the input vector. This must be sorted, complete and linear. Need not be 2:1 balanced.
    @param out the output vector. This will be sorted, complete and linear.
    @param dim the dimension of the tree (1 for binary trees, 2 for quadtrees and 3 for octrees)
    @param maxDepth the maximum depth of the octree must be <= _MAX_LEVEL_   
    */
  int coarsenOctree(const std::vector<TreeNode > &in, std::vector<TreeNode> &out,
      unsigned int dim, unsigned int maxDepth, MPI_Comm comm,
      bool skipPartition = false, MPI_Comm* newCommPtr = NULL, bool* iAmActive = NULL);

  /**
    @author Rahul Sampath
    @brief A naive implementation of coarsenOctree
    @see coarsenOctree
    */
  int simpleCoarsen(const std::vector<TreeNode > &in, std::vector<TreeNode> &out, MPI_Comm comm);

  /**
    @author Rahul Sampath
    @brief A sequential version of coarsenOctree()
    @see coarsenOctree()
    */
  int coarsenOctree(const std::vector<TreeNode > &in, std::vector<TreeNode> &out); 

  /**
    @author Rahul Sampath
    @brief Replaces every octant with its eight children.
    @param in the input vector
    @param out the output vector
    */
  int refineOctree(const std::vector<TreeNode > &in, std::vector<TreeNode> &out); 

  /**
    @author Rahul Sampath
    @brief Replaces every octant with its eight children and partitions the result uniformly across the processors.
    @param in the input vector
    @param out the output vector
    @see refineOctre
    */
  int refineAndPartitionOctree(const std::vector<TreeNode > &in,
      std::vector<TreeNode> &out,MPI_Comm comm); 

  /**
    @author Rahul Sampath
    @brief Creates a regular grid octree.
    @param out the result
    @param lev the level of each octant
    @param dim the dimension of the tree (1 for binary trees, 2 for quadtrees and 3 for octrees)
    @param maxDepth the maximum depth of the octree must be <= _MAX_LEVEL_   
    */
  int createRegularOctree(std::vector<ot::TreeNode>& out, unsigned int lev, 
      unsigned int dim, unsigned int maxDepth, MPI_Comm comm);

  //nodes are partitioned (Used only for p2o) Basically uses a different
  //criteria for selecting blocks than the version used for balancing/meshing.
  //This is because, for balancing and meshing the input is complete. For p2o,
  //the input is typically not complete.
  int blockPartStage1_p2o(std::vector<TreeNode> &nodes, std::vector<TreeNode> &blocks,
      unsigned int dim, unsigned int maxDepth, MPI_Comm comm);

  int blockPartStage2_p2o(std::vector<TreeNode> &nodes, std::vector<TreeNode> &blocks,
      std::vector<ot::TreeNode>& minsAllBlocks, 
      unsigned int dim, unsigned int maxDepth, MPI_Comm comm);

  //nodes are partitioned (Used for balancing and meshing)
  int blockPartStage1(std::vector<TreeNode> &nodes, std::vector<TreeNode> &blocks,
      unsigned int dim, unsigned int maxDepth, MPI_Comm comm);

  //nodes and blocks are partitioned.
  int blockPartStage2(std::vector<TreeNode> &nodes, std::vector<TreeNode> &blocks,
      std::vector<ot::TreeNode>& minsAllBlocks, 
      unsigned int dim, unsigned int maxDepth, MPI_Comm comm);

  int balanceBlocks (const std::vector<TreeNode> &inp,
      const std::vector<TreeNode> &blocks,
      std::vector<TreeNode> &nodes, std::vector<TreeNode> &allBoundaryLeaves,
      bool incCorners, std::vector<unsigned int> *maxBlockBndVec = NULL);

  int comboRipple(std::vector<TreeNode> & in, bool incCorner,
      const unsigned int maxNum = _COMBO_RIPPLE_FACTOR_);

  /**
    @author Rahul Sampath
    @brief Sequential ripple propagation algorithm on linear octrees. The octree need not be complete. It can have holes.
    @param nodes a sorted, linear octree. The octree need not be complete. It can have holes.
    @param incCorners 'true' if you want to balance across corners as well.
    */
  int ripple(std::vector<TreeNode> & nodes, bool incCorners) ;

  /**
    @author Rahul Sampath
    @brief Sequential ripple propagation algorithm using pointer based octree representation (internally). The input (linear) octree need not be complete. It can have holes.
    @param nodes a sorted, linear octree. The octree need not be complete. It can have holes.
    @param incCorners 'true' if you want to balance across corners as well.
    */
  int pointerBasedRipple(std::vector<ot::TreeNode> & nodes, bool incCorners);

  int pickInterProcessorBoundaryNodes(const std::vector<ot::TreeNode> & nodes,
      std::vector<ot::TreeNode>& bndLeaves,
      const ot::TreeNode& firstBlock, const ot::TreeNode& lastBlock);

  int pickInterProcessorBoundaryNodes(const std::vector<ot::TreeNode> & nodes, 
      std::vector<unsigned int >& res,
      const ot::TreeNode& firstBlock, const ot::TreeNode& lastBlock);

  /*
  //old version
  int pickInterProcessorBoundaryNodes(const std::vector<ot::TreeNode> & blocks,
  const std::vector<ot::TreeNode> &nodes,
  std::vector<ot::TreeNode>& bndLeaves,
  unsigned int dim, unsigned int maxDepth);

  //old version
  int pickInterProcessorBoundaryNodes(const std::vector<ot::TreeNode> & blocks,
  const std::vector<ot::TreeNode> &nodes, 
  std::vector<unsigned int >& res,
  unsigned int dim, unsigned int maxDepth);
  */

  char int2char(int dec);
  int int2str(int dec,char* numStr) ;

  /**
    @author Rahul Sampath
    @brief Reads a list of points from a file
    @param filename the file name
    @param pts the points
    */
  int readPtsFromFile(char* filename, std::vector<double>& pts);
    int readNdsFromFile(char* filename, std::vector<double>& pts);
  /**
    @author Ilya Lashuk
    @brief Reads a list of points and corresponding values from a file
    @param filename the file name
    @param pts the points
    @param data the values
    */
  int readDataPtsFromFile(char* filename, std::vector<double>& pts, std::vector<double>& ptVals);

  /**
    @author Rahul Sampath
    @brief Writes a list of points to a file
    @param filename the file name
    @param pts the points
    */
  int writePtsToFile(char* filename,  std::vector<double>& pts);

  /**
    @author Ilya Lashuk
    @brief Writes a list of points and corresponding values from a file
    @param filename the file name
    @param pts the points
    @param data the values
    */
  int writeDataPtsToFile(char* filename, std::vector<double>& pts, std::vector<double>& data);


  /**
    @author Rahul Sampath
    @brief Writes a list of octants to a file
    @param filename the file name
    @param nodes the octants
    */
  int writeNodesToFile (char* filename, const std::vector<TreeNode> & nodes);
  int writeNdsToFile(char* filename, std::vector<double>& pts);


  /**
    @author Rahul Sampath
    @brief Reads a list of octants from a file
    @param filename the file name
    @param nodes the octants
    */
  int readNodesFromFile (char* filename,std::vector<TreeNode > & nodes );

  unsigned int getNodeWeight(const TreeNode * t);

}//end namespace

#endif /*OCTUTILS_H_*/



