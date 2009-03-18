
/**
  @file externVars.h
  @brief Global variables used for profiling the various functions in the library and for storing the stencils.
  @author Rahul S. Sampath
  This file must be included in the driver file (the one containing the main function) only. 
  */

#ifndef __EXTERN_VARS_H__
#define __EXTERN_VARS_H__

#include "petscmat.h"

#ifdef PETSC_USE_LOG

namespace par {

  /** @name Variables to profile the functions in the par module */
  //@{
  int searchEvent;
  int concatEvent;
  int parScatterEvent;
  int gatherEvent;
  int a2avWaitEvent;
  int all2AllvSparseEvent;
  int all2AllvDenseEvent;
  int allGatherEvent;
  int reduceEvent;
  int sendRecvEvent;
  int allReduceEvent;
  int all2AllEvent;
  int allGathervEvent;
  int bcastEvent;
  int scanEvent;
  int partwEvent;
  int sortEvent;
  int remdupEvent;
  int splitComm2wayEvent;
  int splitCommEvent;
  //@}
}

namespace ot {

  /** @name Variables to profile the functions in the oct module */
  //@{
  int mergeRecvKeysBalEvent; 
  int prepBalWlistEvent;
  int prepBalComm1MssgEvent;
  int prepBalComm2MssgEvent;
  int finalBalMergeEvent;
  int mergeComboBalEvent;
  int pickNhBlocksEvent;
  int comboRippleEvent;
  int balSubtreeEvent;
  int completeSubtreeEvent;
  int coarsenEvent;
  int coarsenSeqEvent;
  int simpleCoarsenEvent;
  int balanceEvent;
  int mergeOctreesEvent;
  int rg2oEvent;
  int p2oEvent;
  int p2oSeqEvent;
  int p2oLocalEvent;
  int n2oEvent;
  int n2oSeqEvent;
  int completeRegionEvent;
  int blockPart1Event;
  int blockPart2Event;
  int blockPart3Event;
  int conBalEvent;
  int rippleBalEvent;
  int ptrRippleBalEvent;
  int parRippleType3Event;
  int parRippleType2Event;
  int parRippleType1Event;
  int pickBndEvent;
  int balCommEvent;
  int balScatterEvent;
  int balSplitCommEvent;
  int balBpart1Event;
  int balBpart2Event;
  int markHangingEvent;
  int addBdyEvent;
  //@}

  /** @name Variables to profile the functions in the oda module  */
  //@{
  int DAbPart1Event;
  int DAbPart2Event;
  int DAbPart3Event;
  int DAaprioriCommEvent;
  int buildDaEvent;
  int buildNlistEvent;
  int buildNlistCommEvent;
  int pickGhostsEvent;
  int addBdySiblingsEvent;
  int setMatValuesEvent;
  int readFromGhostNodesBeginEvent;
  int readFromGhostNodesEndEvent;
  int readFromGhostElemsBeginEvent;
  int readFromGhostElemsEndEvent;
  int writeToGhostNodesBeginEvent;
  int writeToGhostNodesEndEvent;
  int writeToGhostElemsBeginEvent;
  int writeToGhostElemsEndEvent;
  int daInitEvent;
  int daFinalEvent;
  //@}

  /** @name Variables to profile the functions in the omg module */
  //@{
  int pcKspShellSetupEvent;
  int pcKspShellDestroyEvent;
  int pcKspShellApplyEvent;
  int setDaEvent;
  int setUpEvent;
  int createRp1Event;
  int createRp2Event;
  int setKspEvent;
  int restrictEvent;
  int dummyRestrictEvent;
  int prolongEvent;
  int scatterEvent;
  int damgInitEvent;
  int damgFinalEvent;
  //@}

  int buildDAstage1Event;
  int buildDAstage2Event;
  int buildDAstage3Event;
  int buildDAstage4Event;
  int buildDAstage5Event;
  int buildDAstage6Event;
  int buildDAstage7Event;
  int buildDAstage8Event;
  int buildDAstage9Event;

  int setDAstage1Event;
  int setDAstage2Event;
  int setDAstage3Event;
  int setDAstage4Event;
  int setDAstage5Event;
  int setDAstage6Event;

  int FLNstage1Event;
  int FLNstage2Event;
  int FLNstage3Event;
  int FLNstage4Event;
  int FLNstage5Event;
  int FLNstage6Event;
  int FLNstage7Event;
  int FLNstage8Event;
  int FLNstage9Event;
  int FLNstage10Event;
  int FLNstage11Event;

}//end namespace

#endif

namespace ot {

  /** @name Variables for storing the various stencils used in the oda and omg module */
  //@{
  double**** RmatType2Stencil = NULL;
  double***** RmatType1Stencil = NULL;
  unsigned short**** VtxMap1 = NULL; 
  unsigned short***** VtxMap2 = NULL; 
  unsigned short***** VtxMap3 = NULL; 
  unsigned short****** VtxMap4 = NULL; 
  double**** ShapeFnCoeffs = NULL; 
  //@}

  /** @name Global Function Handle used in PC_BlockDiag */
  //@{
  void (*getDofAndNodeSizeForPC_BlockDiag)(Mat pcMat,
      unsigned int & dof, unsigned int & nodeSize) = NULL;

  void (*computeInvBlockDiagEntriesForPC_BlockDiag)(Mat pcMat,
      double **invBlockDiagEntries) = NULL;
  //@}

  /**
    @name Global Function Handle for using KSP_Shell
    This will be used at the coarsest grid if not all processors are active on the coarsest grid
    */
  //@{
  void (*getPrivateMatricesForKSP_Shell)(Mat mat,
      Mat *AmatPrivate, Mat *PmatPrivate, MatStructure* pFlag) = NULL;
  //@}
}//end namespace

#endif


