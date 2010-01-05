
/**
  @file sys.C
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  */

#include "petsc.h"
#include "blockDiag.h"
#include "octUtils.h"
#include "parUtils.h"
#include "oda.h"
#include "odaUtils.h"
#include "omg.h"

namespace ot {

#undef __FUNCT__
#define __FUNCT__ "RegisterEvents"
  PetscErrorCode RegisterEvents()
  {
    PetscFunctionBegin;

    PCRegister("blockDiag","${DENDRO_DIR}/lib/libPC","ot::PCCreate_BlockDiag",ot::PCCreate_BlockDiag);

#ifdef PETSC_USE_LOG
    PetscLogEventRegister("splitComm-2way",PETSC_VIEWER_COOKIE,&par::splitComm2wayEvent);
    PetscLogEventRegister("splitComm",PETSC_VIEWER_COOKIE,&par::splitCommEvent);
    PetscLogEventRegister("partW",PETSC_VIEWER_COOKIE,&par::partwEvent);
    PetscLogEventRegister("sort",PETSC_VIEWER_COOKIE,&par::sortEvent);
    PetscLogEventRegister("remdup",PETSC_VIEWER_COOKIE,&par::remdupEvent);
    PetscLogEventRegister("ParallelSearch", PETSC_VIEWER_COOKIE,&par::searchEvent);
    PetscLogEventRegister("Par-Scatter", PETSC_VIEWER_COOKIE,&par::parScatterEvent);
    PetscLogEventRegister("Par-sendRecv", PETSC_VIEWER_COOKIE,&par::sendRecvEvent);
    PetscLogEventRegister("Par-gather", PETSC_VIEWER_COOKIE,&par::gatherEvent);
    PetscLogEventRegister("a2avs-Wait", PETSC_VIEWER_COOKIE,&par::a2avWaitEvent);
    PetscLogEventRegister("Par-a2avs", PETSC_VIEWER_COOKIE,&par::all2AllvSparseEvent);
    PetscLogEventRegister("Par-a2avd", PETSC_VIEWER_COOKIE,&par::all2AllvDenseEvent);
    PetscLogEventRegister("Par-ag", PETSC_VIEWER_COOKIE,&par::allGatherEvent);
    PetscLogEventRegister("Par-agv", PETSC_VIEWER_COOKIE,&par::allGathervEvent);
    PetscLogEventRegister("Par-a2a", PETSC_VIEWER_COOKIE,&par::all2AllEvent);
    PetscLogEventRegister("Par-aRed", PETSC_VIEWER_COOKIE,&par::allReduceEvent);
    PetscLogEventRegister("Par-red", PETSC_VIEWER_COOKIE,&par::reduceEvent);
    PetscLogEventRegister("Par-bcast",	PETSC_VIEWER_COOKIE,&par::bcastEvent);
    PetscLogEventRegister("Par-scan", PETSC_VIEWER_COOKIE,&par::scanEvent);
    PetscLogEventRegister("Par-ConCat",PETSC_VIEWER_COOKIE,&par::concatEvent);
    PetscLogEventRegister("nhBlk",PETSC_VIEWER_COOKIE,&pickNhBlocksEvent);
    PetscLogEventRegister("mergeComboBal",PETSC_VIEWER_COOKIE,&mergeComboBalEvent);
    PetscLogEventRegister("bal-mergeRecv",PETSC_VIEWER_COOKIE,&mergeRecvKeysBalEvent);
    PetscLogEventRegister("finalBalMerge",PETSC_VIEWER_COOKIE,&finalBalMergeEvent);
    PetscLogEventRegister("bal-comm1",PETSC_VIEWER_COOKIE,&prepBalComm1MssgEvent);
    PetscLogEventRegister("bal-comm2",PETSC_VIEWER_COOKIE,&prepBalComm2MssgEvent);
    PetscLogEventRegister("bal-wList",PETSC_VIEWER_COOKIE,&prepBalWlistEvent);
    PetscLogEventRegister("mergeOctrees",PETSC_VIEWER_COOKIE,&mergeOctreesEvent);
    PetscLogEventRegister("rg2o",PETSC_VIEWER_COOKIE,&rg2oEvent);
    PetscLogEventRegister("p2o",PETSC_VIEWER_COOKIE,&p2oEvent);
    PetscLogEventRegister("p2oSeq",PETSC_VIEWER_COOKIE,&p2oSeqEvent);
    PetscLogEventRegister("p2oLoc",PETSC_VIEWER_COOKIE,&p2oLocalEvent);
    PetscLogEventRegister("ComboRipple",PETSC_VIEWER_COOKIE,&comboRippleEvent);
    PetscLogEventRegister("Coarsen",PETSC_VIEWER_COOKIE,&coarsenEvent);
    PetscLogEventRegister("CoarsenSeq",PETSC_VIEWER_COOKIE,&coarsenSeqEvent);
    PetscLogEventRegister("ezCoarse",PETSC_VIEWER_COOKIE,&simpleCoarsenEvent);
    PetscLogEventRegister("bal",PETSC_VIEWER_COOKIE,&balanceEvent);
    PetscLogEventRegister("DA-apriori",PETSC_VIEWER_COOKIE,&DAaprioriCommEvent);
    PetscLogEventRegister("balSubtree", PETSC_VIEWER_COOKIE,&balSubtreeEvent);
    PetscLogEventRegister("conSubtree", PETSC_VIEWER_COOKIE,&completeSubtreeEvent);
    PetscLogEventRegister("n2o",PETSC_VIEWER_COOKIE,&n2oEvent);
    PetscLogEventRegister("n2oSeq",PETSC_VIEWER_COOKIE,&n2oSeqEvent);
    PetscLogEventRegister("compReg",PETSC_VIEWER_COOKIE,&completeRegionEvent);
    PetscLogEventRegister("pickBnd",PETSC_VIEWER_COOKIE,&pickBndEvent);
    PetscLogEventRegister("bPart1",PETSC_VIEWER_COOKIE,&blockPart1Event);
    PetscLogEventRegister("bPart2",PETSC_VIEWER_COOKIE,&blockPart2Event);
    PetscLogEventRegister("bPart3",PETSC_VIEWER_COOKIE,&blockPart3Event);
    PetscLogEventRegister("DAbPart1",PETSC_VIEWER_COOKIE,&DAbPart1Event);
    PetscLogEventRegister("DAbPart2",PETSC_VIEWER_COOKIE,&DAbPart2Event);
    PetscLogEventRegister("DAbPart3",PETSC_VIEWER_COOKIE,&DAbPart3Event);
    PetscLogEventRegister("conBal",PETSC_VIEWER_COOKIE,&conBalEvent);
    PetscLogEventRegister("rippleBal",PETSC_VIEWER_COOKIE,&rippleBalEvent);
    PetscLogEventRegister("ptrRipple",PETSC_VIEWER_COOKIE,&ptrRippleBalEvent);
    PetscLogEventRegister("parRipple1",PETSC_VIEWER_COOKIE,&parRippleType1Event);
    PetscLogEventRegister("parRipple2",PETSC_VIEWER_COOKIE,&parRippleType2Event);
    PetscLogEventRegister("parRipple3",PETSC_VIEWER_COOKIE,&parRippleType3Event);
    PetscLogEventRegister("BuildDA",PETSC_VIEWER_COOKIE,&buildDaEvent);
    PetscLogEventRegister("SetMatValues",PETSC_VIEWER_COOKIE,&setMatValuesEvent);
    PetscLogEventRegister("FlagNodes",PETSC_VIEWER_COOKIE,&markHangingEvent);
    PetscLogEventRegister("BuildNodeList",PETSC_VIEWER_COOKIE,&buildNlistEvent);
    PetscLogEventRegister("Nlist-Comm",PETSC_VIEWER_COOKIE,&buildNlistCommEvent);
    PetscLogEventRegister("Add Bdy",PETSC_VIEWER_COOKIE,&addBdyEvent);
    PetscLogEventRegister("AddBdySiblings",PETSC_VIEWER_COOKIE,&addBdySiblingsEvent);
    PetscLogEventRegister("Pick Ghosts",PETSC_VIEWER_COOKIE,&pickGhostsEvent);
    PetscLogEventRegister("SetDA",PETSC_VIEWER_COOKIE,&setDaEvent);
    PetscLogEventRegister("SetUp",PETSC_VIEWER_COOKIE,&setUpEvent);
    PetscLogEventRegister("crRp1",PETSC_VIEWER_COOKIE,&createRp1Event);
    PetscLogEventRegister("crRp2",PETSC_VIEWER_COOKIE,&createRp2Event);
    PetscLogEventRegister("SetKSP",PETSC_VIEWER_COOKIE,&setKspEvent);
    PetscLogEventRegister("Restrict",PETSC_VIEWER_COOKIE,&restrictEvent);
    PetscLogEventRegister("Restrict-Dummy",PETSC_VIEWER_COOKIE,&dummyRestrictEvent);
    PetscLogEventRegister("Prolong",PETSC_VIEWER_COOKIE,&prolongEvent);
    PetscLogEventRegister("MG Scatter",PETSC_VIEWER_COOKIE,&scatterEvent);
    PetscLogEventRegister("R-Gh-N-Begin",PETSC_VIEWER_COOKIE,&readFromGhostNodesBeginEvent);
    PetscLogEventRegister("R-Gh-N-End",PETSC_VIEWER_COOKIE,&readFromGhostNodesEndEvent);
    PetscLogEventRegister("R-Gh-E-Begin",PETSC_VIEWER_COOKIE,&readFromGhostElemsBeginEvent);
    PetscLogEventRegister("R-Gh-E-End",PETSC_VIEWER_COOKIE,&readFromGhostElemsEndEvent);
    PetscLogEventRegister("W-Gh-N-Begin",PETSC_VIEWER_COOKIE,&writeToGhostNodesBeginEvent);
    PetscLogEventRegister("W-Gh-N-End",PETSC_VIEWER_COOKIE,&writeToGhostNodesEndEvent);
    PetscLogEventRegister("W-Gh-E-Begin",PETSC_VIEWER_COOKIE,&writeToGhostElemsBeginEvent);
    PetscLogEventRegister("W-Gh-E-End",PETSC_VIEWER_COOKIE,&writeToGhostElemsEndEvent);
    PetscLogEventRegister("DA-Init",PETSC_VIEWER_COOKIE,&daInitEvent);
    PetscLogEventRegister("DA-Final",PETSC_VIEWER_COOKIE,&daFinalEvent);
    PetscLogEventRegister("DAMG-Init",PETSC_VIEWER_COOKIE,&damgInitEvent);
    PetscLogEventRegister("DAMG-Final",PETSC_VIEWER_COOKIE,&damgFinalEvent);
    
    PetscLogEventRegister("PC-KSP-Shell-Setup",PETSC_VIEWER_COOKIE,&pcKspShellSetupEvent);
    PetscLogEventRegister("PC-KSP-Shell-Apply",PETSC_VIEWER_COOKIE,&pcKspShellApplyEvent);
    PetscLogEventRegister("PC-KSP-Shell-Destroy",PETSC_VIEWER_COOKIE,&pcKspShellDestroyEvent);

    PetscLogEventRegister("Bal-Comm",PETSC_VIEWER_COOKIE,&balCommEvent);
    PetscLogEventRegister("Bal-SplitComm",PETSC_VIEWER_COOKIE,&balSplitCommEvent);
    PetscLogEventRegister("Bal-Scatter",PETSC_VIEWER_COOKIE,&balScatterEvent);
    PetscLogEventRegister("Bal-Bpart1",PETSC_VIEWER_COOKIE,&balBpart1Event);
    PetscLogEventRegister("Bal-Bpart2",PETSC_VIEWER_COOKIE,&balBpart2Event);

    PetscLogEventRegister("SetDA-stg1",PETSC_VIEWER_COOKIE,&setDAstage1Event);
    PetscLogEventRegister("SetDA-stg2",PETSC_VIEWER_COOKIE,&setDAstage2Event);
    PetscLogEventRegister("SetDA-stg3",PETSC_VIEWER_COOKIE,&setDAstage3Event);
    PetscLogEventRegister("SetDA-stg4",PETSC_VIEWER_COOKIE,&setDAstage4Event);
    PetscLogEventRegister("SetDA-stg5",PETSC_VIEWER_COOKIE,&setDAstage5Event);
    PetscLogEventRegister("SetDA-stg6",PETSC_VIEWER_COOKIE,&setDAstage6Event);

    PetscLogEventRegister("buildDA-stg1",PETSC_VIEWER_COOKIE,&buildDAstage1Event);
    PetscLogEventRegister("buildDA-stg2",PETSC_VIEWER_COOKIE,&buildDAstage2Event);
    PetscLogEventRegister("buildDA-stg3",PETSC_VIEWER_COOKIE,&buildDAstage3Event);
    PetscLogEventRegister("buildDA-stg4",PETSC_VIEWER_COOKIE,&buildDAstage4Event);
    PetscLogEventRegister("buildDA-stg5",PETSC_VIEWER_COOKIE,&buildDAstage5Event);
    PetscLogEventRegister("buildDA-stg6",PETSC_VIEWER_COOKIE,&buildDAstage6Event);
    PetscLogEventRegister("buildDA-stg7",PETSC_VIEWER_COOKIE,&buildDAstage7Event);
    PetscLogEventRegister("buildDA-stg8",PETSC_VIEWER_COOKIE,&buildDAstage8Event);
    PetscLogEventRegister("buildDA-stg9",PETSC_VIEWER_COOKIE,&buildDAstage9Event);

    PetscLogEventRegister("FLN-stg1",PETSC_VIEWER_COOKIE,&FLNstage1Event);
    PetscLogEventRegister("FLN-stg2",PETSC_VIEWER_COOKIE,&FLNstage2Event);
    PetscLogEventRegister("FLN-stg3",PETSC_VIEWER_COOKIE,&FLNstage3Event);
    PetscLogEventRegister("FLN-stg4",PETSC_VIEWER_COOKIE,&FLNstage4Event);
    PetscLogEventRegister("FLN-stg5",PETSC_VIEWER_COOKIE,&FLNstage5Event);
    PetscLogEventRegister("FLN-stg6",PETSC_VIEWER_COOKIE,&FLNstage6Event);
    PetscLogEventRegister("FLN-stg7",PETSC_VIEWER_COOKIE,&FLNstage7Event);
    PetscLogEventRegister("FLN-stg8",PETSC_VIEWER_COOKIE,&FLNstage8Event);
    PetscLogEventRegister("FLN-stg9",PETSC_VIEWER_COOKIE,&FLNstage9Event);
    PetscLogEventRegister("FLN-stg10",PETSC_VIEWER_COOKIE,&FLNstage10Event);
    PetscLogEventRegister("FLN-stg11",PETSC_VIEWER_COOKIE,&FLNstage11Event);
#endif

    PetscFunctionReturn(0);
  }

}//end namespace

