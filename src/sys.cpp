
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

    // PCRegister("blockDiag","${DENDRO_DIR}/lib/libPC","ot::PCCreate_BlockDiag", ot::PCCreate_BlockDiag);
    // PCRegister("blockDiag", ot::PCCreate_BlockDiag);

#ifdef PETSC_USE_LOG
    PetscClassId classid;
    PetscClassIdRegister("Dendro",&classid);

    PetscLogEventRegister("splitComm-2way", classid, &par::splitComm2wayEvent);
    PetscLogEventRegister("splitComm",classid,&par::splitCommEvent);
    PetscLogEventRegister("partW",classid,&par::partwEvent);
    PetscLogEventRegister("sort",classid,&par::sortEvent);
    PetscLogEventRegister("remdup",classid,&par::remdupEvent);
    PetscLogEventRegister("ParallelSearch", classid,&par::searchEvent);
    PetscLogEventRegister("Par-Scatter", classid,&par::parScatterEvent);
    PetscLogEventRegister("Par-sendRecv", classid,&par::sendRecvEvent);
    PetscLogEventRegister("Par-gather", classid,&par::gatherEvent);
    PetscLogEventRegister("a2avs-Wait", classid,&par::a2avWaitEvent);
    PetscLogEventRegister("Par-a2avs", classid,&par::all2AllvSparseEvent);
    PetscLogEventRegister("Par-a2avd", classid,&par::all2AllvDenseEvent);
    PetscLogEventRegister("Par-ag", classid,&par::allGatherEvent);
    PetscLogEventRegister("Par-agv", classid,&par::allGathervEvent);
    PetscLogEventRegister("Par-a2a", classid,&par::all2AllEvent);
    PetscLogEventRegister("Par-aRed", classid,&par::allReduceEvent);
    PetscLogEventRegister("Par-red", classid,&par::reduceEvent);
    PetscLogEventRegister("Par-bcast",	classid,&par::bcastEvent);
    PetscLogEventRegister("Par-scan", classid,&par::scanEvent);
    PetscLogEventRegister("Par-ConCat",classid,&par::concatEvent);
    PetscLogEventRegister("nhBlk",classid,&pickNhBlocksEvent);
    PetscLogEventRegister("mergeComboBal",classid,&mergeComboBalEvent);
    PetscLogEventRegister("bal-mergeRecv",classid,&mergeRecvKeysBalEvent);
    PetscLogEventRegister("finalBalMerge",classid,&finalBalMergeEvent);
    PetscLogEventRegister("bal-comm1",classid,&prepBalComm1MssgEvent);
    PetscLogEventRegister("bal-comm2",classid,&prepBalComm2MssgEvent);
    PetscLogEventRegister("bal-wList",classid,&prepBalWlistEvent);
    PetscLogEventRegister("mergeOctrees",classid,&mergeOctreesEvent);
    PetscLogEventRegister("rg2o",classid,&rg2oEvent);
    PetscLogEventRegister("p2o",classid,&p2oEvent);
    PetscLogEventRegister("p2oSeq",classid,&p2oSeqEvent);
    PetscLogEventRegister("p2oLoc",classid,&p2oLocalEvent);
    PetscLogEventRegister("ComboRipple",classid,&comboRippleEvent);
    PetscLogEventRegister("Coarsen",classid,&coarsenEvent);
    PetscLogEventRegister("CoarsenSeq",classid,&coarsenSeqEvent);
    PetscLogEventRegister("ezCoarse",classid,&simpleCoarsenEvent);
    PetscLogEventRegister("bal",classid,&balanceEvent);
    PetscLogEventRegister("DA-apriori",classid,&DAaprioriCommEvent);
    PetscLogEventRegister("balSubtree", classid,&balSubtreeEvent);
    PetscLogEventRegister("conSubtree", classid,&completeSubtreeEvent);
    PetscLogEventRegister("n2o",classid,&n2oEvent);
    PetscLogEventRegister("n2oSeq",classid,&n2oSeqEvent);
    PetscLogEventRegister("compReg",classid,&completeRegionEvent);
    PetscLogEventRegister("pickBnd",classid,&pickBndEvent);
    PetscLogEventRegister("bPart1",classid,&blockPart1Event);
    PetscLogEventRegister("bPart2",classid,&blockPart2Event);
    PetscLogEventRegister("bPart3",classid,&blockPart3Event);
    PetscLogEventRegister("DAbPart1",classid,&DAbPart1Event);
    PetscLogEventRegister("DAbPart2",classid,&DAbPart2Event);
    PetscLogEventRegister("DAbPart3",classid,&DAbPart3Event);
    PetscLogEventRegister("conBal",classid,&conBalEvent);
    PetscLogEventRegister("rippleBal",classid,&rippleBalEvent);
    PetscLogEventRegister("ptrRipple",classid,&ptrRippleBalEvent);
    PetscLogEventRegister("parRipple1",classid,&parRippleType1Event);
    PetscLogEventRegister("parRipple2",classid,&parRippleType2Event);
    PetscLogEventRegister("parRipple3",classid,&parRippleType3Event);
    PetscLogEventRegister("BuildDA",classid,&buildDaEvent);
    PetscLogEventRegister("SetMatValues",classid,&setMatValuesEvent);
    PetscLogEventRegister("FlagNodes",classid,&markHangingEvent);
    PetscLogEventRegister("BuildNodeList",classid,&buildNlistEvent);
    PetscLogEventRegister("Nlist-Comm",classid,&buildNlistCommEvent);
    PetscLogEventRegister("Add Bdy",classid,&addBdyEvent);
    PetscLogEventRegister("AddBdySiblings",classid,&addBdySiblingsEvent);
    PetscLogEventRegister("Pick Ghosts",classid,&pickGhostsEvent);
    PetscLogEventRegister("SetDA",classid,&setDaEvent);
    PetscLogEventRegister("SetUp",classid,&setUpEvent);
    PetscLogEventRegister("crRp1",classid,&createRp1Event);
    PetscLogEventRegister("crRp2",classid,&createRp2Event);
    PetscLogEventRegister("SetKSP",classid,&setKspEvent);
    PetscLogEventRegister("Restrict",classid,&restrictEvent);
    PetscLogEventRegister("Restrict-Dummy",classid,&dummyRestrictEvent);
    PetscLogEventRegister("Prolong",classid,&prolongEvent);
    PetscLogEventRegister("MG Scatter",classid,&scatterEvent);
    PetscLogEventRegister("R-Gh-N-Begin",classid,&readFromGhostNodesBeginEvent);
    PetscLogEventRegister("R-Gh-N-End",classid,&readFromGhostNodesEndEvent);
    PetscLogEventRegister("R-Gh-E-Begin",classid,&readFromGhostElemsBeginEvent);
    PetscLogEventRegister("R-Gh-E-End",classid,&readFromGhostElemsEndEvent);
    PetscLogEventRegister("W-Gh-N-Begin",classid,&writeToGhostNodesBeginEvent);
    PetscLogEventRegister("W-Gh-N-End",classid,&writeToGhostNodesEndEvent);
    PetscLogEventRegister("W-Gh-E-Begin",classid,&writeToGhostElemsBeginEvent);
    PetscLogEventRegister("W-Gh-E-End",classid,&writeToGhostElemsEndEvent);
    PetscLogEventRegister("DA-Init",classid,&daInitEvent);
    PetscLogEventRegister("DA-Final",classid,&daFinalEvent);
    PetscLogEventRegister("DAMG-Init",classid,&damgInitEvent);
    PetscLogEventRegister("DAMG-Final",classid,&damgFinalEvent);
    
    PetscLogEventRegister("PC-KSP-Shell-Setup",classid,&pcKspShellSetupEvent);
    PetscLogEventRegister("PC-KSP-Shell-Apply",classid,&pcKspShellApplyEvent);
    PetscLogEventRegister("PC-KSP-Shell-Destroy",classid,&pcKspShellDestroyEvent);

    PetscLogEventRegister("Bal-Comm",classid,&balCommEvent);
    PetscLogEventRegister("Bal-SplitComm",classid,&balSplitCommEvent);
    PetscLogEventRegister("Bal-Scatter",classid,&balScatterEvent);
    PetscLogEventRegister("Bal-Bpart1",classid,&balBpart1Event);
    PetscLogEventRegister("Bal-Bpart2",classid,&balBpart2Event);

    PetscLogEventRegister("SetDA-stg1",classid,&setDAstage1Event);
    PetscLogEventRegister("SetDA-stg2",classid,&setDAstage2Event);
    PetscLogEventRegister("SetDA-stg3",classid,&setDAstage3Event);
    PetscLogEventRegister("SetDA-stg4",classid,&setDAstage4Event);
    PetscLogEventRegister("SetDA-stg5",classid,&setDAstage5Event);
    PetscLogEventRegister("SetDA-stg6",classid,&setDAstage6Event);

    PetscLogEventRegister("buildDA-stg1",classid,&buildDAstage1Event);
    PetscLogEventRegister("buildDA-stg2",classid,&buildDAstage2Event);
    PetscLogEventRegister("buildDA-stg3",classid,&buildDAstage3Event);
    PetscLogEventRegister("buildDA-stg4",classid,&buildDAstage4Event);
    PetscLogEventRegister("buildDA-stg5",classid,&buildDAstage5Event);
    PetscLogEventRegister("buildDA-stg6",classid,&buildDAstage6Event);
    PetscLogEventRegister("buildDA-stg7",classid,&buildDAstage7Event);
    PetscLogEventRegister("buildDA-stg8",classid,&buildDAstage8Event);
    PetscLogEventRegister("buildDA-stg9",classid,&buildDAstage9Event);

    PetscLogEventRegister("FLN-stg1",classid,&FLNstage1Event);
    PetscLogEventRegister("FLN-stg2",classid,&FLNstage2Event);
    PetscLogEventRegister("FLN-stg3",classid,&FLNstage3Event);
    PetscLogEventRegister("FLN-stg4",classid,&FLNstage4Event);
    PetscLogEventRegister("FLN-stg5",classid,&FLNstage5Event);
    PetscLogEventRegister("FLN-stg6",classid,&FLNstage6Event);
    PetscLogEventRegister("FLN-stg7",classid,&FLNstage7Event);
    PetscLogEventRegister("FLN-stg8",classid,&FLNstage8Event);
    PetscLogEventRegister("FLN-stg9",classid,&FLNstage9Event);
    PetscLogEventRegister("FLN-stg10",classid,&FLNstage10Event);
    PetscLogEventRegister("FLN-stg11",classid,&FLNstage11Event);
#endif

    PetscFunctionReturn(0);
  }

}//end namespace

