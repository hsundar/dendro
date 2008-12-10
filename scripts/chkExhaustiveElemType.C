
#include<cassert>

enum exhaustiveElemType {
  N = 0,
  Y = 1,
  X = 2,
  XY = 3,
  Z = 4,
  ZY = 5,
  ZX = 6,
  ZXY = 7,
  XY_XY = 11,
  XY_ZXY = 15,
  YZ_ZY = 21,
  YZ_ZXY = 23,
  YZ_XY_ZXY = 31,
  ZX_ZX = 38,
  ZX_ZXY = 39,
  ZX_XY_ZXY = 47,
  ZX_YZ_ZXY = 55,
  ZX_YZ_XY_ZXY = 63
};

//flags MUST be in the order: 2,1,4,3,6,5
char getExhaustiveElemType(bool flags[6]) {
  char type = 0;
  for(int i=0;i<6;i++){
    type = ( (flags[i]) ? ( type | (1 << i) ) : type );
  }
  return type;
}//end fn2

int main() {
  //5,6,3,4,1,2
  bool v[18][6] = {
    {0,0,0,0,0,0},
    {0,0,0,0,0,1},
    {0,0,0,0,1,0},
    {0,0,0,0,1,1},
    {0,0,0,1,0,0},
    {0,0,0,1,0,1},
    {0,0,0,1,1,0},
    {0,0,0,1,1,1},
    {0,0,1,0,1,1},
    {0,0,1,1,1,1},
    {0,1,0,1,0,1},
    {0,1,0,1,1,1},
    {0,1,1,1,1,1},
    {1,0,0,1,1,0},
    {1,0,0,1,1,1},
    {1,0,1,1,1,1},
    {1,1,0,1,1,1},
    {1,1,1,1,1,1},
  };

  char types[18];
  for(int i=0;i<18;i++) {
    bool flags[6] = {v[i][5], v[i][4], v[i][3], v[i][2], v[i][1],v[i][0]};
    types[i] = getExhaustiveElemType(flags);
  }
  
  assert( N == types[0]);
  assert( Y ==  types[1]);
  assert( X ==  types[2]);
  assert( XY ==  types[3]);
  assert( Z ==  types[4]);
  assert( ZY ==  types[5]);
  assert( ZX ==  types[6]);
  assert( ZXY ==  types[7]);
  assert( XY_XY ==  types[8]);
  assert( XY_ZXY ==  types[9]);
  assert( YZ_ZY ==  types[10]);
  assert( YZ_ZXY ==  types[11]);
  assert( YZ_XY_ZXY ==  types[12]);
  assert( ZX_ZX ==  types[13]);
  assert( ZX_ZXY ==  types[14]);
  assert( ZX_XY_ZXY ==  types[15]);
  assert( ZX_YZ_ZXY ==  types[16]);
  assert( ZX_YZ_XY_ZXY ==  types[17]);

  assert( N != types[2]);
  assert( Y !=  types[0]);
  assert( X !=  types[3]);
  assert( XY !=  types[1]);
  assert( Z !=  types[6]);
  assert( ZY !=  types[4]);
  assert( ZX !=  types[8]);
  assert( ZXY !=  types[2]);
  assert( XY_XY !=  types[3]);
  assert( XY_ZXY !=  types[4]);
  assert( YZ_ZY !=  types[17]);
  assert( YZ_ZXY !=  types[10]);
  assert( YZ_XY_ZXY !=  types[2]);
  assert( ZX_ZX !=  types[3]);
  assert( ZX_ZXY !=  types[4]);
  assert( ZX_XY_ZXY !=  types[5]);
  assert( ZX_YZ_ZXY !=  types[6]);
  assert( ZX_YZ_XY_ZXY !=  types[7]);


}
