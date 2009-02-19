
/**
 * @file odaUtils.txx
 * @author		Rahul S. Sampath, rahul.sampath@gmail.com 
 * @brief A list of non-member templated functions for the ot::DA class.
 **/ 

#include <cassert>
#include "cnumEtypes.h"

namespace ot { 

  template<unsigned char cNum>
    unsigned char getElemType(unsigned char hnMask) {
      unsigned char eType = 0;
      switch(hnMask) {
        case  ot::cNumEtype<cNum>::ET_N: {
                                           eType = 0;
                                           break;
                                         }
        case  ot::cNumEtype<cNum>::ET_Y: {
                                           eType = 1;
                                           break;
                                         }
        case  ot::cNumEtype<cNum>::ET_X: {
                                           eType = 2;
                                           break;
                                         }
        case  ot::cNumEtype<cNum>::ET_XY: {
                                            eType = 3;
                                            break;
                                          }
        case  ot::cNumEtype<cNum>::ET_Z: {
                                           eType = 4;
                                           break;
                                         }
        case  ot::cNumEtype<cNum>::ET_ZY: {
                                            eType = 5;
                                            break;
                                          }
        case  ot::cNumEtype<cNum>::ET_ZX: {
                                            eType = 6;
                                            break;
                                          }
        case  ot::cNumEtype<cNum>::ET_ZXY: {
                                             eType = 7;
                                             break;
                                           }
        case  ot::cNumEtype<cNum>::ET_XY_XY: {
                                               eType = 8;
                                               break;
                                             }
        case  ot::cNumEtype<cNum>::ET_XY_ZXY: {
                                                eType = 9;
                                                break;
                                              }
        case  ot::cNumEtype<cNum>::ET_YZ_ZY: {
                                               eType = 10;
                                               break;
                                             }
        case  ot::cNumEtype<cNum>::ET_YZ_ZXY: {
                                                eType = 11;
                                                break;
                                              }
        case  ot::cNumEtype<cNum>::ET_YZ_XY_ZXY: {
                                                   eType = 12;
                                                   break;
                                                 }
        case  ot::cNumEtype<cNum>::ET_ZX_ZX: {
                                               eType = 13;
                                               break;
                                             }
        case  ot::cNumEtype<cNum>::ET_ZX_ZXY: {
                                                eType = 14;
                                                break;
                                              }
        case  ot::cNumEtype<cNum>::ET_ZX_XY_ZXY: {
                                                   eType = 15;
                                                   break;
                                                 }
        case  ot::cNumEtype<cNum>::ET_ZX_YZ_ZXY: {
                                                   eType = 16;
                                                   break;
                                                 }
        case  ot::cNumEtype<cNum>::ET_ZX_YZ_XY_ZXY: {
                                                      eType = 17;
                                                      break;
                                                    }
        default:  assert(false);
      }
      return eType; 
    }//end fn;

  template <typename T>
    void injectNodalVector(ot::DA* dac, ot::DA* daf, unsigned int dof,
        std::vector<T>& fVec, std::vector<T>& cVec, void (*setZero)(T&)) {

      dac->createVector<T>(cVec, false, false, dof);

      for(int i = 0; i < cVec.size(); i++) {
        (*setZero)(cVec[i]);
      }

      T* cArr = NULL;
      T* fArr = NULL;
      dac->vecGetBuffer<T>(cVec, cArr, false, false, false, dof);
      daf->vecGetBuffer<T>(fVec, fArr, false, false, true, dof);

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

}//end namespace


