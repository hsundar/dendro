
/**
 * @file cnumEtypes.h
 * @brief 		The set of templated classes for performing
 child-number based operations.
 * @author		Rahul S. Sampath, rahul.sampath@gmail.com
 * 
 */ 

#ifndef __CNUM_ETYPE_H__
#define __CNUM_ETYPE_H__

namespace ot  {

  template<unsigned char cNum>
    class cNumEtype; 

  //Cnum-0
  template<>
    class cNumEtype<0> {
      public:
        enum eTypes {
          ET_N                   = 0,
          ET_Y                   = 4 ,
          ET_X                   = 2 ,
          ET_XY                  = 6 ,
          ET_Z                  = 16 ,
          ET_ZY                  = 20 ,
          ET_ZX                  = 18 ,
          ET_ZXY                 = 22 ,
          ET_XY_XY               = 14 ,
          ET_XY_ZXY              = 30 ,
          ET_YZ_ZY               = 84 ,
          ET_YZ_ZXY              = 86 ,
          ET_YZ_XY_ZXY           = 94 ,
          ET_ZX_ZX               = 50 ,
          ET_ZX_ZXY              = 54 ,
          ET_ZX_XY_ZXY           = 62 ,
          ET_ZX_YZ_ZXY           = 118 ,
          ET_ZX_YZ_XY_ZXY        = 126 
        };
    };

  //Cnum-1
  template<>
    class cNumEtype<1> {
      public:
        enum eTypes {
          ET_N                   = 0 , 
          ET_Y                   = 1 ,
          ET_X                   = 8 ,
          ET_XY                  = 9 ,
          ET_Z                  = 32 ,
          ET_ZY                  = 33 ,
          ET_ZX                  = 40 ,
          ET_ZXY                 = 41 ,
          ET_XY_XY               = 13 ,
          ET_XY_ZXY              = 45 ,
          ET_YZ_ZY               = 49 ,
          ET_YZ_ZXY              = 57 ,
          ET_YZ_XY_ZXY           = 61 ,
          ET_ZX_ZX               = 168 ,
          ET_ZX_ZXY              = 169 ,
          ET_ZX_XY_ZXY           = 173 ,
          ET_ZX_YZ_ZXY           = 185 ,
          ET_ZX_YZ_XY_ZXY        = 189 
        };
    };


  //Cnum-2
  template<>
    class cNumEtype<2> {
      public:
        enum eTypes {
          ET_N                   = 0 , 
          ET_Y                   = 8 ,
          ET_X                   = 1 ,
          ET_XY                  = 9 ,
          ET_Z                  = 64 ,
          ET_ZY                  = 72 ,
          ET_ZX                  = 65 ,
          ET_ZXY                 = 73 ,
          ET_XY_XY               = 11 ,
          ET_XY_ZXY              = 75 ,
          ET_YZ_ZY               = 200 ,
          ET_YZ_ZXY              = 201 ,
          ET_YZ_XY_ZXY           = 203 ,
          ET_ZX_ZX               = 81 ,
          ET_ZX_ZXY              = 89 ,
          ET_ZX_XY_ZXY           = 91 ,
          ET_ZX_YZ_ZXY           = 217 ,
          ET_ZX_YZ_XY_ZXY        = 219 
        };
    };


  //Cnum-3
  template<>
    class cNumEtype<3> {
      public:
        enum eTypes {
          ET_N                   = 0 , 
          ET_Y                   = 2 ,
          ET_X                   = 4 ,
          ET_XY                  = 6 ,
          ET_Z                  = 128 ,
          ET_ZY                  = 130 ,
          ET_ZX                  = 132 ,
          ET_ZXY                 = 134 ,
          ET_XY_XY               = 7 ,
          ET_XY_ZXY              = 135 ,
          ET_YZ_ZY               = 162 ,
          ET_YZ_ZXY              = 166 ,
          ET_YZ_XY_ZXY           = 167 ,
          ET_ZX_ZX               = 196 ,
          ET_ZX_ZXY              = 198 ,
          ET_ZX_XY_ZXY           = 199 ,
          ET_ZX_YZ_ZXY           = 230 ,
          ET_ZX_YZ_XY_ZXY        = 231 
        };
    };


  //Cnum-4
  template<>
    class cNumEtype<4> {
      public:
        enum eTypes {
          ET_N                   = 0 , 
          ET_Y                   = 1 ,
          ET_X                   = 32 ,
          ET_XY                  = 33 ,
          ET_Z                  = 64 ,
          ET_ZY                  = 65 ,
          ET_ZX                  = 96 ,
          ET_ZXY                 = 97 ,
          ET_XY_XY               = 35 ,
          ET_XY_ZXY              = 99 ,
          ET_YZ_ZY               = 69 ,
          ET_YZ_ZXY              = 101 ,
          ET_YZ_XY_ZXY           = 103 ,
          ET_ZX_ZX               = 224 ,
          ET_ZX_ZXY              = 225 ,
          ET_ZX_XY_ZXY           = 227 ,
          ET_ZX_YZ_ZXY           = 229 ,
          ET_ZX_YZ_XY_ZXY        = 231 
        };
    };


  //Cnum-5
  template<>
    class cNumEtype<5> {
      public:
        enum eTypes {
          ET_N                   = 0 , 
          ET_Y                   = 2 ,
          ET_X                   = 128 ,
          ET_XY                  = 130 ,
          ET_Z                  = 16 ,
          ET_ZY                  = 18 ,
          ET_ZX                  = 144 ,
          ET_ZXY                 = 146 ,
          ET_XY_XY               = 138 ,
          ET_XY_ZXY              = 154 ,
          ET_YZ_ZY               = 19 ,
          ET_YZ_ZXY              = 147 ,
          ET_YZ_XY_ZXY           = 155 ,
          ET_ZX_ZX               = 208 ,
          ET_ZX_ZXY              = 210 ,
          ET_ZX_XY_ZXY           = 218 ,
          ET_ZX_YZ_ZXY           = 211 ,
          ET_ZX_YZ_XY_ZXY        = 219 
        };
    };


  //Cnum-6
  template<>
    class cNumEtype<6> {
      public:
        enum eTypes {
          ET_N                  = 0 , 
          ET_Y                  = 4 ,
          ET_X                  = 16 ,
          ET_XY                 = 20 ,
          ET_Z                 = 128 ,
          ET_ZY                 = 132 ,
          ET_ZX                 = 144 ,
          ET_ZXY                = 148 ,
          ET_XY_XY              = 21 ,
          ET_XY_ZXY             = 149 ,
          ET_YZ_ZY              = 140 ,
          ET_YZ_ZXY             = 156 ,
          ET_YZ_XY_ZXY          = 157 ,
          ET_ZX_ZX              = 176 ,
          ET_ZX_ZXY             = 180 ,
          ET_ZX_XY_ZXY          = 181 ,
          ET_ZX_YZ_ZXY          = 188 ,
          ET_ZX_YZ_XY_ZXY       = 189 
        };
    };

  //Cnum-7
  template<>
    class cNumEtype<7> {
      public:
        enum eTypes {
          ET_N                   = 0 , 
          ET_Y                   = 8 ,
          ET_X                   = 64 ,
          ET_XY                  = 72 ,
          ET_Z                  = 32 ,
          ET_ZY                  = 40 ,
          ET_ZX                  = 96 ,
          ET_ZXY                 = 104 ,
          ET_XY_XY               = 76 ,
          ET_XY_ZXY              = 108 ,
          ET_YZ_ZY               = 42 ,
          ET_YZ_ZXY              = 106 ,
          ET_YZ_XY_ZXY           = 110 ,
          ET_ZX_ZX               = 112 ,
          ET_ZX_ZXY              = 120 ,
          ET_ZX_XY_ZXY           = 124 ,
          ET_ZX_YZ_ZXY           = 122 ,
          ET_ZX_YZ_XY_ZXY        = 126 
        };
    };

} //end namespace

#endif

