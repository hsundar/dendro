
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

}//end namespace


