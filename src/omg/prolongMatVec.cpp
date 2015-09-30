
/**
  @file ProlongMatVec.C
  @brief Prolongation MatVec
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

  PetscErrorCode addProlongMatVec(Mat R, Vec v1, Vec v2, Vec v3)
  {
    PetscScalar one = 1.0;
    PetscFunctionBegin;
    if((v2!=v3) && (v1!=v3)) {
      //Note This will fail only if v2==v3 or v1 ==v3!
      //(i.e they are identical copies pointing to the same memory location)
      iC(MatMultTranspose(R, v1, v3));//v3 = R*v1
      iC(VecAXPY(v3,one,v2));//v3 = v3+ v2=v2 + R*v1
    }else {
      //This is less efficient but failproof.
      TransferOpData *data;			
      MatShellGetContext( R, (void **)&data);
      Vec tmp = data->addPtmp;
      if(tmp == NULL) {
        VecDuplicate(v3,&tmp);
        data->addPtmp = tmp;
      }
      iC(MatMultTranspose(R, v1, tmp));//tmp=R'*v1;
      iC(VecWAXPY(v3,one,v2,tmp));//v3 = (1*v2)+tmp=v2 + R'*v1
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode prolongMatVecType2(Mat R, Vec c, Vec f) {		
    TransferOpData *data;			
    PetscFunctionBegin;
    iC(MatShellGetContext( R, (void **)&data));
    MPI_Comm comm = data->comm;
    Vec tmp = data->tmp;				
    PetscInt tmpSz;
    PetscInt fSz;
    iC(VecGetLocalSize(tmp,&tmpSz));
    iC(VecGetLocalSize(f,&fSz));
    prolongMatVecType1(R,c,tmp);		
    scatterValues(tmp, f, tmpSz, fSz, data->sendSzP, data->sendOffP, 
        data->recvSzP, data->recvOffP, comm);
    PetscFunctionReturn(0);
  }

#define ITLB_SET_VALUE_NO_SUPPRESSED_DOFS {\
  farr[fidx+l] += (Rval*carr[cidx+l]);\
}

#define ITLB_SET_VALUE_SUPPRESSED_DOFS {\
  if(!( (suppressedDOFc && suppressedDOFc[cidx+l]) ||\
        (suppressedDOFf && suppressedDOFf[fidx+l]) )) {\
    farr[fidx+l] += (Rval*carr[cidx+l]);\
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
    double** type2RmatPtr =  RmatType2Stencil[cNumCoarse][ctype];\
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
      double** type1RmatPtr =  RmatType1Stencil[cNumCoarse][cNumFine][ctype];\
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

PetscErrorCode prolongMatVecType1(Mat R, Vec c, Vec f) {

  PROF_MG_PROLONG_BEGIN 

    TransferOpData *data;
  iC(MatShellGetContext(R, (void **)&data));

  unsigned int dof = data->dof;

  //Overlap 60% of independent computation with the first communication and
  //40% with the second communication. Since the message size for the fist 
  //communication will roughly be 1.25 times that for the second communication
  //This is because in the first communication, we exchange both coarse and
  //fine grid ghosts. In the second, we only exchange fine grid ghosts.
  //I estimate the number of fine grid ghosts to be 4 times that of the coarse
  //grid ghosts (1 surface, uniform refinement)
  unsigned int fop = 60;

  ot::DA * dac = data->dac;
  ot::DA * daf = data->daf;	

  unsigned char* suppressedDOFc = data->suppressedDOFc;
  unsigned char* suppressedDOFf = data->suppressedDOFf;

  PetscInt cSz;
  iC(VecGetLocalSize(c,&cSz));

  unsigned int fopCnt = (fop*cSz)/(100*dof);

  PetscScalar *farr = NULL;
  PetscScalar *carr = NULL;

  std::vector<ot::FineTouchedStatus >* fineTouchedFlags = data->fineTouchedFlags;
  ot::FineTouchedStatus* fineTouchedFlagsArr;

  dac->vecGetBuffer(c, carr, false, false, true, dof);//Read-only
  daf->vecGetBuffer<ot::FineTouchedStatus >(*fineTouchedFlags,
      fineTouchedFlagsArr, false, false, true, 1);//read-only 

  if(dac->iAmActive()) {
    dac->ReadFromGhostsBegin<PetscScalar>(carr, dof);		
  }

  if(daf->iAmActive()) {
    daf->ReadFromGhostsBegin<ot::FineTouchedStatus>(fineTouchedFlagsArr, 1);
  }

  VecZeroEntries(f);
  daf->vecGetBuffer(f, farr, false, false, false, dof);//Writable

  //Note: If Coarse is Independent, then the corresponding Fine is also independent.
  //Hence, overlapping comm with comp is possible.		
  //Order of the test condition is important. We want to store
  //the info before checking loopCtr.	

  if(dac->iAmActive()) {
    unsigned int loopCtr = 0;
    if(suppressedDOFc || suppressedDOFf) {
      for(dac->init<ot::DA_FLAGS::INDEPENDENT>(), daf->init<ot::DA_FLAGS::WRITABLE>();
          ( (daf->currWithInfo() == daf->currWithInfo()) && 
            (dac->currWithInfo() < dac->end<ot::DA_FLAGS::INDEPENDENT>()) &&
            (loopCtr < fopCnt) );
          dac->next<ot::DA_FLAGS::INDEPENDENT>(), loopCtr++) {
        INTERGRID_TRANSFER_LOOP_BLOCK(ITLB_SET_VALUE_SUPPRESSED_DOFS)
      }//end Independent loop (overlapping with read from coarse ghosts)
    } else {
      for(dac->init<ot::DA_FLAGS::INDEPENDENT>(), daf->init<ot::DA_FLAGS::WRITABLE>();
          ( (daf->currWithInfo() == daf->currWithInfo()) && 
            (dac->currWithInfo() < dac->end<ot::DA_FLAGS::INDEPENDENT>()) &&
            (loopCtr < fopCnt) );
          dac->next<ot::DA_FLAGS::INDEPENDENT>(), loopCtr++) {
        INTERGRID_TRANSFER_LOOP_BLOCK(ITLB_SET_VALUE_NO_SUPPRESSED_DOFS)
      }//end Independent loop (overlapping with read from coarse ghosts)
    }
  }

  if(dac->iAmActive()) {
    dac->ReadFromGhostsEnd<PetscScalar>(carr);
  }

  if(daf->iAmActive()) {
    daf->ReadFromGhostsEnd<ot::FineTouchedStatus>(fineTouchedFlagsArr);
  }

  if(dac->iAmActive()) {
    if(suppressedDOFc || suppressedDOFf) {
      for(dac->init<ot::DA_FLAGS::W_DEPENDENT>(), daf->init<ot::DA_FLAGS::WRITABLE>();
          dac->curr() < dac->end<ot::DA_FLAGS::W_DEPENDENT>(); dac->next<ot::DA_FLAGS::W_DEPENDENT>()) {
        INTERGRID_TRANSFER_LOOP_BLOCK(ITLB_SET_VALUE_SUPPRESSED_DOFS)
      }//end dependent loop
    } else {
      for(dac->init<ot::DA_FLAGS::W_DEPENDENT>(), daf->init<ot::DA_FLAGS::WRITABLE>();
          dac->curr() < dac->end<ot::DA_FLAGS::W_DEPENDENT>(); dac->next<ot::DA_FLAGS::W_DEPENDENT>()) {
        INTERGRID_TRANSFER_LOOP_BLOCK(ITLB_SET_VALUE_NO_SUPPRESSED_DOFS)
      }//end dependent loop
    }
  }

  if(daf->iAmActive()) {
    daf->WriteToGhostsBegin<PetscScalar>(farr, dof);
  }

  if(dac->iAmActive()) {
    //Continue Independent loop from where we left off.
    if(suppressedDOFc || suppressedDOFf) {
      for(dac->init<ot::DA_FLAGS::FROM_STORED>(), daf->init<ot::DA_FLAGS::FROM_STORED>();
          dac->curr() < dac->end<ot::DA_FLAGS::INDEPENDENT>(); dac->next<ot::DA_FLAGS::INDEPENDENT>()) {
        INTERGRID_TRANSFER_LOOP_BLOCK(ITLB_SET_VALUE_SUPPRESSED_DOFS)
      }//end Independent loop (overlapping with write to fine ghosts) 
    } else {
      for(dac->init<ot::DA_FLAGS::FROM_STORED>(), daf->init<ot::DA_FLAGS::FROM_STORED>();
          dac->curr() < dac->end<ot::DA_FLAGS::INDEPENDENT>(); dac->next<ot::DA_FLAGS::INDEPENDENT>()) {
        INTERGRID_TRANSFER_LOOP_BLOCK(ITLB_SET_VALUE_NO_SUPPRESSED_DOFS)
      }//end Independent loop (overlapping with write to fine ghosts) 
    }
  }

  if(daf->iAmActive()) {
    daf->WriteToGhostsEnd<PetscScalar>(farr, dof);
  }

  daf->vecRestoreBuffer(f, farr, false, false, false, dof);//Writable 
  dac->vecRestoreBuffer(c, carr, false, false, true, dof);//Read-only
  daf->vecRestoreBuffer<ot::FineTouchedStatus >(*fineTouchedFlags, fineTouchedFlagsArr, false, false, true, 1);//read-only 

#ifdef PETSC_USE_LOG
  PetscLogFlops(128*dof*(daf->getElementSize()));
#endif

  PROF_MG_PROLONG_END 
}//end Prolong-3

#undef ITLB_SET_VALUE_NO_SUPPRESSED_DOFS
#undef ITLB_SET_VALUE_SUPPRESSED_DOFS
#undef INTERGRID_TRANSFER_LOOP_BLOCK

}//end namespace


