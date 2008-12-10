
/**
  @file RestrictMatVec.C
  @brief Restriction MatVec
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  */

#include "petsc.h"
#include "petscmat.h"
#include "omg.h"
#include "oda.h"

#ifndef iC
#define iC(fun) {CHKERRQ(fun);}
#endif

#ifdef __DEBUG__
#ifndef __DEBUG_MG__
#define __DEBUG_MG__
#endif
#endif

namespace ot {

  extern double **** RmatType2Stencil;
  extern double ***** RmatType1Stencil;
  extern unsigned short**** VtxMap1; 
  extern unsigned short***** VtxMap2; 
  extern unsigned short***** VtxMap3; 
  extern unsigned short****** VtxMap4; 

  PetscErrorCode  addRestrictMatVec(Mat R, Vec v1, Vec v2, Vec v3)	
  {
    PetscScalar one = 1.0;
    PetscFunctionBegin;
    if((v2!=v3) && (v1!=v3)) {
      //Note This will fail only if v2==v3 or v1 ==v3!(i.e they are identical copies pointing to the same memory location)
      iC(MatMult(R, v1, v3));//v3 = R*v1
      iC(VecAXPY(v3,one,v2));//v3 = v3+ v2=v2 + R*v1
    }else {
      //This is less efficient but failproof.
      TransferOpData *data;
      MatShellGetContext( R, (void **)&data);
      Vec tmp = data->addRtmp;
      if(tmp == NULL) {
        VecDuplicate(v3,&tmp);
        data->addRtmp = tmp;
      }
      iC(MatMult(R, v1, tmp));//tmp=R*v1;
      iC(VecWAXPY(v3,one,v2,tmp));//v3 = (1*v2)+tmp=v2 + R*v1
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode   restrictMatVecType2(Mat R, Vec f, Vec c) {
    TransferOpData *data;
    PetscFunctionBegin;
    iC(MatShellGetContext( R, (void **)&data));
    MPI_Comm comm = data->comm;
    Vec tmp = data->tmp;		
    PetscInt tmpSz;
    PetscInt fSz;
    iC(VecGetLocalSize(tmp,&tmpSz));
    iC(VecGetLocalSize(f,&fSz));
    scatterValues(f, tmp, fSz, tmpSz, data->sendSzR,
        data->sendOffR, data->recvSzR, data->recvOffR, comm);
    restrictMatVecType1(R, tmp, c);
    PetscFunctionReturn(0);
  }

#define ITLB_SET_VALUE_NO_SUPPRESSED_DOFS {\
  carr[cidx+l] += (Rval*farr[fidx+l]);\
}

#define ITLB_SET_VALUE_SUPPRESSED_DOFS {\
  if(!( (suppressedDOFf && suppressedDOFf[fidx+l]) ||\
        (suppressedDOFc && suppressedDOFc[cidx+l]) )) {\
    carr[cidx+l] += (Rval*farr[fidx+l]);\
  }\
}

#define INTERGRID_TRANSFER_LOOP_BLOCK(ITLB_SET_VALUE) {\
  /*The fine loop is always Writable, but the coarse loop*/\
  /*could be Independent or W_Dependent. Hence the fine counter must*/\
  /*be incremented properly to align with the coarse.*/\
  Point Cpt = dac->getCurrentOffset();\
  while(daf->getCurrentOffset() != Cpt) {\
    if(daf->isLUTcompressed()) {\
      daf->updateQuotientCounter();\
    }\
    daf->next<ot::DA_FLAGS::WRITABLE>();\
  } \
  unsigned char chnMask = dac->getHangingNodeIndex(dac->curr());\
  unsigned char cNumCoarse = dac->getChildNumber();\
  unsigned int cIndices[8];\
  dac->getNodeIndices(cIndices);\
  unsigned char ctype = 0;\
  GET_ETYPE_BLOCK(ctype,chnMask,cNumCoarse)\
  if(daf->getLevel(daf->curr()) == dac->getLevel(dac->curr())) {\
    /*The coarse and fine elements are the same,*/\
    /*so cNumCoarse = cNumFine. This is type-2*/\
    double** type2RmatPtr = RmatType2Stencil[cNumCoarse][ctype];\
    unsigned char fhnMask = daf->getHangingNodeIndex(daf->curr());\
    unsigned int fIndices[8];\
    daf->getNodeIndices(fIndices);\
    for(unsigned char fCtr = 0; fCtr < 8; fCtr++) {\
      if(!(fhnMask & (1 << fCtr))) {\
        ot::FineTouchedStatus* fineTouchedStatusPtr = (&(fineTouchedFlagsArr[fIndices[fCtr]]));\
        unsigned int fidx = fIndices[fCtr]*dof;\
        for(unsigned char cCtr = 0; cCtr < 8; cCtr++) {\
          /*Read fineTouchedFlagsArr Directly*/\
          if( (fineTouchedStatusPtr->flags[fCtr]) & (1 << cCtr) ) {\
            unsigned int cidx = cIndices[cCtr]*dof;\
            double Rval = type2RmatPtr[cCtr][fCtr];\
            for(unsigned int l = 0; l < dof; l++) {\
              ITLB_SET_VALUE\
            }\
          }\
        }\
      }\
    }\
    daf->next<ot::DA_FLAGS::WRITABLE>();\
  }else {\
    for(unsigned char cNumFine = 0; cNumFine < 8; cNumFine++) {\
      /*The coarse and fine elements are NOT the same. This is type-1.*/\
      /*Loop over each of the 8 children of the coarse element.*/\
      /*These are the underlying fine elements.*/\
      double** type1RmatPtr = RmatType1Stencil[cNumCoarse][cNumFine][ctype];\
      unsigned char fhnMask = daf->getHangingNodeIndex(daf->curr());\
      unsigned int fIndices[8];\
      daf->getNodeIndices(fIndices);\
      for(unsigned int fCtr = 0; fCtr < 8; fCtr++) {\
        if(!(fhnMask & (1 << fCtr))) {\
          ot::FineTouchedStatus* fineTouchedStatusPtr = (&(fineTouchedFlagsArr[fIndices[fCtr]]));\
          unsigned int fidx = fIndices[fCtr]*dof;\
          for(unsigned int cCtr = 0; cCtr < 8; cCtr++) {\
            /*Read fineTouchedFlagsArr Directly*/\
            if( (fineTouchedStatusPtr->flags[fCtr]) & (1 << cCtr) ) {\
              unsigned int cidx = cIndices[cCtr]*dof;\
              double Rval = type1RmatPtr[cCtr][fCtr];\
              for(unsigned int l = 0; l < dof; l++) {\
                ITLB_SET_VALUE\
              }\
            }\
          }\
        }\
      }\
      daf->next<ot::DA_FLAGS::WRITABLE>();\
    }\
  }\
}

#define INTERGRID_TRANSFER_LOOP_BLOCK_DUMMY {\
  /*The fine loop is always Writable, but the coarse loop*/\
  /*could be Independent or W_Dependent. Hence the fine counter must*/\
  /*be incremented properly to align with the coarse.*/\
  Point Cpt = dac->getCurrentOffset();\
  while(daf->getCurrentOffset() != Cpt) {\
    if(daf->isLUTcompressed()) {\
      daf->updateQuotientCounter();\
    }\
    daf->next<ot::DA_FLAGS::WRITABLE>();\
  }\
  unsigned char chnMask = dac->getHangingNodeIndex(dac->curr());\
  unsigned char cNumCoarse = dac->getChildNumber();\
  unsigned int cIndices[8];\
  dac->getNodeIndices(cIndices);\
  unsigned char ctype = 0;\
  GET_ETYPE_BLOCK(ctype,chnMask,cNumCoarse)\
  if(daf->getLevel(daf->curr()) == dac->getLevel(dac->curr())) {\
    /*The coarse and fine elements are the same,*/\
    /*so cNumCoarse = cNumFine. This is type-2*/\
    unsigned char fhnMask = daf->getHangingNodeIndex(daf->curr());\
    unsigned int fIndices[8];\
    daf->getNodeIndices(fIndices);\
    for(unsigned char fCtr = 0; fCtr < 8; fCtr++) {\
      if(!(fhnMask & (1 << fCtr))) {\
        unsigned char thisElemLev = daf->getLevel(daf->curr());\
        unsigned char refElemLev = daf->getLevel(fIndices[fCtr]);\
        ot::FineTouchedDummyStatus* fineTouchedDummyStatusPtr = \
        (&(fineTouchedDummyFlagsArr[fIndices[fCtr]]));\
        if(thisElemLev == refElemLev) {\
          /*Set the dummy flag for this fine element and fine node pair*/\
          /*Since, I deal with this one byte at a time. It is consistent*/\
          /* across all computers. Endian issues do not enter. Hence, I*/\
          /* choose to manipulate chars instead of shorts. Besides, the*/\
          /* size of short is not guaranteed to be 2 bytes on all*/\
          /* machines.*/\
          /* A char is always 1 byte.*/\
          /*1 bit: Entry set or not (flags[2*fCtr][0]) */\
          /*1 bit: scalingCtr (flags[2*fCtr][1]) */\
          /*2 bits: stencilType (flags[2*fCtr][3,2]) */\
          /*3 bits: cNumFine (flags[2*fCtr][6,5,4]) */\
          /*3 bits: cNumCoarse (flags[2*fCtr+1][2,1,0]) */\
          /*5 bits: ctype (flags[2*fCtr+1][7,6,5,4,3]) */\
          fineTouchedDummyStatusPtr->flags[fCtr<<1] = 1;\
          fineTouchedDummyStatusPtr->flags[(fCtr<<1)+1] = cNumCoarse;\
          fineTouchedDummyStatusPtr->flags[(fCtr<<1)+1] |= (ctype<<3);\
        }else {\
          unsigned char scalingCtr = 0;\
          if(thisElemLev < refElemLev) {\
            scalingCtr = 1;\
          }\
          fineTouchedDummyStatusPtr->flags[fCtr<<1] = 1;\
          fineTouchedDummyStatusPtr->flags[fCtr<<1] |= (scalingCtr<<1);\
          fineTouchedDummyStatusPtr->flags[fCtr<<1] |= (2<<2);\
          fineTouchedDummyStatusPtr->flags[(fCtr<<1)+1] = cNumCoarse;\
          fineTouchedDummyStatusPtr->flags[(fCtr<<1)+1] |= (ctype<<3);\
        }\
      }\
    }\
    daf->next<ot::DA_FLAGS::WRITABLE>();\
  }else {\
    for(unsigned char cNumFine = 0; cNumFine < 8; cNumFine++) {\
      /*The coarse and fine elements are NOT the same. This is type-1.*/\
      /*Loop over each of the 8 children of the coarse element.*/\
      /*These are the underlying fine elements.*/\
      unsigned char fhnMask = daf->getHangingNodeIndex(daf->curr());\
      unsigned int fIndices[8];\
      daf->getNodeIndices(fIndices);\
      for(unsigned int fCtr = 0; fCtr < 8; fCtr++) {\
        if(!(fhnMask & (1 << fCtr))) {\
          unsigned char thisElemLev = daf->getLevel(daf->curr());\
          unsigned char refElemLev = daf->getLevel(fIndices[fCtr]);\
          ot::FineTouchedDummyStatus* fineTouchedDummyStatusPtr = \
          (&(fineTouchedDummyFlagsArr[fIndices[fCtr]]));\
          if(thisElemLev == refElemLev) {\
            fineTouchedDummyStatusPtr->flags[fCtr<<1] = 1;\
            fineTouchedDummyStatusPtr->flags[fCtr<<1] |= (1 << 2);\
            fineTouchedDummyStatusPtr->flags[fCtr<<1] |= (cNumFine<<4);\
            fineTouchedDummyStatusPtr->flags[(fCtr<<1)+1] = cNumCoarse;\
            fineTouchedDummyStatusPtr->flags[(fCtr<<1)+1] |= (ctype<<3);\
          }else {\
            unsigned char scalingCtr = 0;\
            if(thisElemLev < refElemLev) {\
              scalingCtr = 1;\
            }\
            fineTouchedDummyStatusPtr->flags[fCtr<<1] = 1;\
            fineTouchedDummyStatusPtr->flags[fCtr<<1] |= (scalingCtr<<1);\
            fineTouchedDummyStatusPtr->flags[fCtr<<1] |= (3<<2);\
            fineTouchedDummyStatusPtr->flags[fCtr<<1] |= (cNumFine<<4);\
            fineTouchedDummyStatusPtr->flags[(fCtr<<1)+1] = cNumCoarse;\
            fineTouchedDummyStatusPtr->flags[(fCtr<<1)+1] |= (ctype<<3);\
          }\
        }\
      }\
      daf->next<ot::DA_FLAGS::WRITABLE>();\
    }\
  }\
}

#define ITLB_DUMMY_FCTR_BLOCK1 {\
  if(fineTouchedDummyStatusPtr->flags[dummyFctr<<1]) {\
    /*12_10 = (00001100)_2*/\
    unsigned char dummyStencilType = \
    ((12 & (fineTouchedDummyStatusPtr->flags[dummyFctr<<1]))>>2);\
    /*7_10 = (00000111)_2*/\
    unsigned char dummyCnumCoarse = \
    (7 & (fineTouchedDummyStatusPtr->flags[(dummyFctr<<1)+1]));\
    /*248_10 = (11111000)_2*/\
    unsigned char dummyCtype = \
    ((248 & (fineTouchedDummyStatusPtr->flags[(dummyFctr<<1)+1]))>>3);\
    switch(dummyStencilType) {\
      case 0: {\
                dummyMapPtrs[dummyFctr] = \
                VtxMap1[dummyFctr][dummyCnumCoarse][dummyCtype];\
                break;\
              }\
      case 1: {\
                /*112_10 = (01110000)_2*/\
                unsigned char dummyCnumFine = \
                ((112 & (fineTouchedDummyStatusPtr->flags[dummyFctr<<1]))>>4);\
                dummyMapPtrs[dummyFctr] = \
                VtxMap2[dummyFctr][dummyCnumFine][dummyCnumCoarse][dummyCtype];\
                break;\
              }\
      case 2: {\
                /*2_10 = (00000010)_2*/\
                unsigned char dummyScalingCtr = \
                ((2 & (fineTouchedDummyStatusPtr->flags[dummyFctr<<1]))>>1);\
                dummyMapPtrs[dummyFctr] = \
                VtxMap3[dummyFctr-1][dummyScalingCtr][dummyCnumCoarse][dummyCtype];\
                break;\
              }\
      case 3: {\
                unsigned char dummyScalingCtr = \
                ((2 & (fineTouchedDummyStatusPtr->flags[dummyFctr<<1]))>>1);\
                unsigned char dummyCnumFine = \
                ((112 & (fineTouchedDummyStatusPtr->flags[dummyFctr<<1]))>>4);\
                dummyMapPtrs[dummyFctr] = \
                VtxMap4[dummyFctr-1][dummyScalingCtr][dummyCnumFine][dummyCnumCoarse][dummyCtype];\
                break;\
              }\
      default: {\
                 assert(false);\
               }\
    }\
  }\
}

#define ITLB_DUMMY_FCTR_BLOCK2 {\
  if(dummyMapPtrs[dummyFctr]) {\
    for(unsigned char dummyCctr = 0; dummyCctr < 8; dummyCctr++) {\
      if(coarseVtxId == dummyMapPtrs[dummyFctr][dummyCctr]) {\
        skipThisEntry = true;\
        break;\
      }\
    }\
    if(skipThisEntry) {\
      break;\
    }\
  }\
}

#define ITLB_DUMMY_FINAL_SET_VALUE(nodeNum,idx) {\
  if(!(fhnMask & (1 << nodeNum))) {\
    ot::FineTouchedDummyStatus* fineTouchedDummyStatusPtr = (&(fineTouchedDummyFlagsArr[idx]));\
    typedef unsigned short* ushPtr;\
    ushPtr dummyMapPtrs[8];\
    for(unsigned char dummyFctr = 0; dummyFctr < 8; dummyFctr++) {\
      dummyMapPtrs[dummyFctr] = NULL;\
      ITLB_DUMMY_FCTR_BLOCK1\
    }\
    ot::FineTouchedStatus* fineTouchedStatusPtr = (&(fineTouchedFlagsArr[idx]));\
    for(unsigned char fCtr = 0; fCtr < 8; fCtr++) {\
      fineTouchedStatusPtr->flags[fCtr] = 0;\
      /*Handle negative boundaries*/\
      if(dummyMapPtrs[fCtr]) {\
        for(unsigned char cCtr = 0; cCtr < 8; cCtr++) {\
          /*Can't guarantee that each of the 8 elements*/\
          /*surrounding this fine node has a different cNum*/\
          /*But, they will have a different fCtr*/\
          unsigned short coarseVtxId = dummyMapPtrs[fCtr][cCtr];\
          bool skipThisEntry = false;\
          for(unsigned char dummyFctr = 0; dummyFctr < fCtr; dummyFctr++) {\
            ITLB_DUMMY_FCTR_BLOCK2\
          }\
          if(!skipThisEntry) {\
            fineTouchedStatusPtr->flags[fCtr] |= (1 << cCtr);\
          }\
        }\
      }\
    }\
  }\
}

#define INTERGRID_TRANSFER_LOOP_BLOCK_DUMMY_FINAL_W {\
  /*To avoid redundant writes, only the element whose anchor is */\
  /*the regular fine grid node writes for all the 8 elements */\
  /*surrounding this node. Note, that some of these elements */\
  /* may be owned by other processors. So all processors take */\
  /* care of all the 8 elements surrounding the node they own.*/\
  /* The only other problem is with the positive boundary nodes. */\
  /* The element whose anchor is this positive boundary node is only */\
  /* a pseudo-element and will never be visited while looping through the */\
  /* elements. To make things worse, we can have situtations where a */\
  /* positive boundary node is owned by one processor and all the */\
  /* true  elements that share this node are on different processors. */\
  /* Thus a writable loop will never suffice to take care of this */\
  /* scenario. Hence, this WRITABLE loop will handle all nodes except */\
  /* positive boundaries and a separate ALL loop will handle positive */\
  /* boundary nodes alone. NOTE, that unlike the other loops in the */\
  /* restriction/prolongation, these two loops are not simultaneous */\
  /* loops through both the fine and coarse grids. Looping through the */\
  /* fine mesh will suffice.*/\
  unsigned char fhnMask = daf->getHangingNodeIndex(daf->curr());\
  ITLB_DUMMY_FINAL_SET_VALUE(0,daf->curr())\
}

#define INTERGRID_TRANSFER_LOOP_BLOCK_DUMMY_FINAL_A {\
  unsigned char fBndFlag;\
  bool fIsBnd = daf->isBoundaryOctant(&fBndFlag);\
  fIsBnd = (fBndFlag > ot::TreeNode::NEG_POS_DEMARCATION);\
  if(fIsBnd) {\
    unsigned char fhnMask = daf->getHangingNodeIndex(daf->curr());\
    unsigned int fIndices[8];\
    daf->getNodeIndices(fIndices);\
    if(fBndFlag & ot::TreeNode::X_POS_BDY) {\
      ITLB_DUMMY_FINAL_SET_VALUE(1,fIndices[1])\
    }\
    if(fBndFlag & ot::TreeNode::Y_POS_BDY) {\
      ITLB_DUMMY_FINAL_SET_VALUE(2,fIndices[2])\
    }\
    if(fBndFlag & ot::TreeNode::Z_POS_BDY) {\
      ITLB_DUMMY_FINAL_SET_VALUE(4,fIndices[4])\
    }\
    if( (fBndFlag & (ot::TreeNode::X_POS_BDY + ot::TreeNode::Y_POS_BDY))\
        == (ot::TreeNode::X_POS_BDY + ot::TreeNode::Y_POS_BDY) ) {\
      ITLB_DUMMY_FINAL_SET_VALUE(3,fIndices[3])\
    }\
    if( (fBndFlag & (ot::TreeNode::Y_POS_BDY + ot::TreeNode::Z_POS_BDY))\
        == (ot::TreeNode::Y_POS_BDY + ot::TreeNode::Z_POS_BDY) ) {\
      ITLB_DUMMY_FINAL_SET_VALUE(6,fIndices[6])\
    }\
    if( (fBndFlag & (ot::TreeNode::Z_POS_BDY + ot::TreeNode::X_POS_BDY))\
        == (ot::TreeNode::Z_POS_BDY + ot::TreeNode::X_POS_BDY) ) {\
      ITLB_DUMMY_FINAL_SET_VALUE(5,fIndices[5])\
    }\
    if( (fBndFlag & (ot::TreeNode::X_POS_BDY + ot::TreeNode::Y_POS_BDY + ot::TreeNode::Z_POS_BDY))\
        == (ot::TreeNode::X_POS_BDY + ot::TreeNode::Y_POS_BDY + ot::TreeNode::Z_POS_BDY) ) {\
      ITLB_DUMMY_FINAL_SET_VALUE(7,fIndices[7])\
    }\
  }else {\
    if(daf->isLUTcompressed()) {\
      daf->updateQuotientCounter();\
    }\
  }\
}

PetscErrorCode dummyRestrictMatVecType1(TransferOpData *data) {

  PROF_MG_RESTRICT_DUMMY_BEGIN

    ot::DA * dac = data->dac;
  ot::DA * daf = data->daf;

  ot::FineTouchedDummyStatus* fineTouchedDummyFlagsArr;
  std::vector<ot::FineTouchedDummyStatus > fineTouchedDummyFlags;

  daf->createVector<ot::FineTouchedDummyStatus >(fineTouchedDummyFlags, false, false, 1);
  daf->vecGetBuffer<ot::FineTouchedDummyStatus >(fineTouchedDummyFlags,
      fineTouchedDummyFlagsArr, false, false, false, 1);//writable 

  if(dac->iAmActive()) {
    for(dac->init<ot::DA_FLAGS::W_DEPENDENT>(), daf->init<ot::DA_FLAGS::WRITABLE>();
        dac->curr() < dac->end<ot::DA_FLAGS::W_DEPENDENT>(); dac->next<ot::DA_FLAGS::W_DEPENDENT>()) {
      INTERGRID_TRANSFER_LOOP_BLOCK_DUMMY;	
    }//end dependent loop
  }

  if(daf->iAmActive()) {
    daf->WriteToGhostsBegin<ot::FineTouchedDummyStatus>(fineTouchedDummyFlagsArr, 1);
  }

  if(dac->iAmActive()) {
    //Note: If Coarse is Independent, then the corresponding Fine is also independent.
    //Hence, overlapping comm with comp is possible.		
    for(dac->init<ot::DA_FLAGS::INDEPENDENT>(), daf->init<ot::DA_FLAGS::WRITABLE>();
        dac->curr() < dac->end<ot::DA_FLAGS::INDEPENDENT>(); dac->next<ot::DA_FLAGS::INDEPENDENT>()) {
      INTERGRID_TRANSFER_LOOP_BLOCK_DUMMY;	
    }//end Independent loop (overlapping with write to coarse ghosts) 
  }

  if(daf->iAmActive()) {
    daf->WriteToGhostsEnd<ot::FineTouchedDummyStatus >(fineTouchedDummyFlagsArr, 1);
  }

  //Take care of the discrepancies across processors.
  //It is not sufficient to loop over the dependent elements alone to set the
  //status. i.e. setting status for the independent elements can not be
  //combined with computing dummystatus. This is because by definition, a
  //dependent element is one which has atleast 1 writable and 1 ghost node. But
  //we could have cases where the owner of the node is not a dependent element,
  //but this node is shared with ghost elements.

  ot::FineTouchedStatus* fineTouchedFlagsArr;
  std::vector<ot::FineTouchedStatus >* fineTouchedFlags = data->fineTouchedFlags;

  daf->vecGetBuffer<ot::FineTouchedStatus >(*fineTouchedFlags, fineTouchedFlagsArr,
      false, false, false, 1);//writable 

  if(daf->iAmActive()) {
    for(daf->init<ot::DA_FLAGS::WRITABLE>(); daf->curr() < daf->end<ot::DA_FLAGS::WRITABLE>(); daf->next<ot::DA_FLAGS::WRITABLE>()) {
      INTERGRID_TRANSFER_LOOP_BLOCK_DUMMY_FINAL_W;	
    }//end  W loop
    for(daf->init<ot::DA_FLAGS::ALL>(); daf->curr() < daf->end<ot::DA_FLAGS::ALL>(); daf->next<ot::DA_FLAGS::ALL>()) {
      INTERGRID_TRANSFER_LOOP_BLOCK_DUMMY_FINAL_A;	
    }//end  A loop
  }

  daf->vecRestoreBuffer<ot::FineTouchedStatus >(*fineTouchedFlags, fineTouchedFlagsArr, 
      false, false, false, 1);//writable 

  //THIS IS A HACK FOR EFFICIENCY PURPOSES. Although, the buffer was modified
  //there is no need to write the changes back to the vector. This is because
  //the dummy vector is only temporary.   
  daf->vecRestoreBuffer<ot::FineTouchedDummyStatus >(fineTouchedDummyFlags, fineTouchedDummyFlagsArr, 
      false, false, true, 1);//READ-ONLY

  fineTouchedDummyFlags.clear();

  PROF_MG_RESTRICT_DUMMY_END
}//restrict-3

PetscErrorCode restrictMatVecType1(Mat R, Vec f, Vec c) {

  PROF_MG_RESTRICT_BEGIN

    TransferOpData *data;
  iC(MatShellGetContext( R, (void **)&data));

  unsigned int dof = data->dof;

  //Overlap 75% of the independent computation with the first communication and
  //25% with the second communication. In the first communication, we exchange
  // fine grid ghosts. In the second, we exchange coarse grid ghosts (1/4 of
  //fine grid, assuming uniform refinement). So, the
  //first comm. is more expensive.
  unsigned int fop = 75;

  unsigned char* suppressedDOFc = data->suppressedDOFc;
  unsigned char* suppressedDOFf = data->suppressedDOFf;

  ot::DA * dac = data->dac;
  ot::DA * daf = data->daf;

  PetscInt cSz;
  iC(VecGetLocalSize(c,&cSz));

  unsigned int fopCnt = (fop*cSz)/(100*dof);

  //unsigned int fopCnt = data->minIndependentSize;

  std::vector<ot::FineTouchedStatus >* fineTouchedFlags = data->fineTouchedFlags;
  ot::FineTouchedStatus* fineTouchedFlagsArr;

  PetscScalar *farr = NULL;
  PetscScalar *carr = NULL;

  daf->vecGetBuffer(f,farr,false,false,true,dof);//Read-only
  daf->vecGetBuffer<ot::FineTouchedStatus >(*fineTouchedFlags, 
      fineTouchedFlagsArr, false, false, true, 1);//read-only 

  if(daf->iAmActive()) {
    daf->ReadFromGhostsBegin<PetscScalar>(farr, dof);
    daf->ReadFromGhostsBegin<ot::FineTouchedStatus>(fineTouchedFlagsArr, 1);
  }

  VecZeroEntries(c);
  dac->vecGetBuffer(c,carr,false,false,false,dof);//Writable

  if(dac->iAmActive()) {
    //Note: If Coarse is Independent, then the corresponding Fine is also independent.
    //Hence, overlapping comm with comp is possible.		
    //Order of the test condition is important. We want to store the info before checking loopCtr.		 
    unsigned int loopCtr = 0;
    if(suppressedDOFc || suppressedDOFf) {
      for(dac->init<ot::DA_FLAGS::INDEPENDENT>(), daf->init<ot::DA_FLAGS::WRITABLE>();
          ( (daf->currWithInfo() == daf->currWithInfo()) && 
            (dac->currWithInfo() < dac->end<ot::DA_FLAGS::INDEPENDENT>()) && (loopCtr < fopCnt) );
          dac->next<ot::DA_FLAGS::INDEPENDENT>(), loopCtr++) {
        INTERGRID_TRANSFER_LOOP_BLOCK(ITLB_SET_VALUE_SUPPRESSED_DOFS);	
      }//end Independent loop (overlapping with read from fine ghosts)
    } else {
      for(dac->init<ot::DA_FLAGS::INDEPENDENT>(), daf->init<ot::DA_FLAGS::WRITABLE>();
          ( (daf->currWithInfo() == daf->currWithInfo()) && 
            (dac->currWithInfo() < dac->end<ot::DA_FLAGS::INDEPENDENT>()) && (loopCtr < fopCnt) );
          dac->next<ot::DA_FLAGS::INDEPENDENT>(), loopCtr++) {
        INTERGRID_TRANSFER_LOOP_BLOCK(ITLB_SET_VALUE_NO_SUPPRESSED_DOFS);	
      }//end Independent loop (overlapping with read from fine ghosts)
    }
  }

  if(daf->iAmActive()) {
    daf->ReadFromGhostsEnd<PetscScalar>(farr);
    daf->ReadFromGhostsEnd<ot::FineTouchedStatus>(fineTouchedFlagsArr);
  }

  if(dac->iAmActive()) {
    if(suppressedDOFc || suppressedDOFf) {
      for(dac->init<ot::DA_FLAGS::W_DEPENDENT>(), daf->init<ot::DA_FLAGS::WRITABLE>();
          dac->curr() < dac->end<ot::DA_FLAGS::W_DEPENDENT>(); dac->next<ot::DA_FLAGS::W_DEPENDENT>()) {
        INTERGRID_TRANSFER_LOOP_BLOCK(ITLB_SET_VALUE_SUPPRESSED_DOFS);	
      }//end dependent loop
    } else {
      for(dac->init<ot::DA_FLAGS::W_DEPENDENT>(), daf->init<ot::DA_FLAGS::WRITABLE>();
          dac->curr() < dac->end<ot::DA_FLAGS::W_DEPENDENT>(); dac->next<ot::DA_FLAGS::W_DEPENDENT>()) {
        INTERGRID_TRANSFER_LOOP_BLOCK(ITLB_SET_VALUE_NO_SUPPRESSED_DOFS);	
      }//end dependent loop
    }
  }

  if(dac->iAmActive()) {
    dac->WriteToGhostsBegin<PetscScalar>(carr,  dof);
  }

  if(dac->iAmActive()) {
    //Continue Independent loop from where we left off.
    if(suppressedDOFc || suppressedDOFf) {
      for(dac->init<ot::DA_FLAGS::FROM_STORED>(), daf->init<ot::DA_FLAGS::FROM_STORED>();
          dac->curr() < dac->end<ot::DA_FLAGS::INDEPENDENT>(); dac->next<ot::DA_FLAGS::INDEPENDENT>()) {
        INTERGRID_TRANSFER_LOOP_BLOCK(ITLB_SET_VALUE_SUPPRESSED_DOFS);	
      }//end Independent loop (overlapping with write to coarse ghosts) 
    } else {
      for(dac->init<ot::DA_FLAGS::FROM_STORED>(), daf->init<ot::DA_FLAGS::FROM_STORED>();
          dac->curr() < dac->end<ot::DA_FLAGS::INDEPENDENT>(); dac->next<ot::DA_FLAGS::INDEPENDENT>()) {
        INTERGRID_TRANSFER_LOOP_BLOCK(ITLB_SET_VALUE_NO_SUPPRESSED_DOFS);	
      }//end Independent loop (overlapping with write to coarse ghosts) 
    }
  }

  if(dac->iAmActive()) {
    dac->WriteToGhostsEnd<PetscScalar>(carr, dof);
  }

  daf->vecRestoreBuffer(f,farr,false,false,true,dof);//Read-only
  dac->vecRestoreBuffer(c,carr,false,false,false,dof);//Writable  
  daf->vecRestoreBuffer<ot::FineTouchedStatus >(*fineTouchedFlags, 
      fineTouchedFlagsArr, false, false, true, 1);//read-only 

#ifdef PETSC_USE_LOG
  PetscLogFlops(128*dof*(daf->getElementSize()));
#endif

  PROF_MG_RESTRICT_END
}//restrict-3

#undef ITLB_SET_VALUE_NO_SUPPRESSED_DOFS
#undef ITLB_SET_VALUE_SUPPRESSED_DOFS
#undef INTERGRID_TRANSFER_LOOP_BLOCK
#undef INTERGRID_TRANSFER_LOOP_BLOCK_DUMMY
#undef ITLB_DUMMY_FCTR_BLOCK1 
#undef ITLB_DUMMY_FCTR_BLOCK2 
#undef ITLB_DUMMY_FINAL_SET_VALUE 
#undef INTERGRID_TRANSFER_LOOP_BLOCK_DUMMY_FINAL_W
#undef INTERGRID_TRANSFER_LOOP_BLOCK_DUMMY_FINAL_A

}//end namespace

