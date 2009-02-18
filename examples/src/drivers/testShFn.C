
#include<iostream>

/*
   Shape Functions for ChildNumber type=1 and hangingType=1
   */
#define phi0(x,y,z) (0.1875 - (0.1875*(x)) - (0.0625*(y)) - (0.1875*(z))\
 + (0.0625*(x)*(y)) + (0.0625*(y)*(z)) + (0.1875*(z)*(x)) - (0.0625*(x)*(y)*(z)))

#define phi1(x,y,z) (0.125 + (0.125*(x)) - (0.125*(y)) - (0.125*(z))\
 - (0.125*(x)*(y)) + (0.125*(y)*(z)) - (0.125*(z)*(x)) + (0.125*(x)*(y)*(z)))

#define phi2(x,y,z) (0.0625 - (0.0625*(x)) + (0.0625*(y)) - (0.0625*(z))\
 - (0.0625*(x)*(y)) - (0.0625*(y)*(z)) + (0.0625*(z)*(x)) + (0.0625*(x)*(y)*(z)))

#define phi3(x,y,z) (0.125 + (0.125*(x)) + (0.125*(y)) - (0.125*(z))\
 + (0.125*(x)*(y)) - (0.125*(y)*(z)) - (0.125*(z)*(x)) - (0.125*(x)*(y)*(z)))

#define phi4(x,y,z) (0.125 - (0.125*(x)) - (0.125*(y)) + (0.125*(z))\
 + (0.125*(x)*(y)) - (0.125*(y)*(z)) - (0.125*(z)*(x)) + (0.125*(x)*(y)*(z)))

#define phi5(x,y,z) (0.125 + (0.125*(x)) - (0.125*(y)) + (0.125*(z))\
 - (0.125*(x)*(y)) - (0.125*(y)*(z)) + (0.125*(z)*(x)) - (0.125*(x)*(y)*(z)))

#define phi6(x,y,z) (0.125 - (0.125*(x)) + (0.125*(y)) + (0.125*(z))\
 - (0.125*(x)*(y)) + (0.125*(y)*(z)) - (0.125*(z)*(x)) - (0.125*(x)*(y)*(z)))

#define phi7(x,y,z) (0.125 + (0.125*(x)) + (0.125*(y)) + (0.125*(z))\
 + (0.125*(x)*(y)) + (0.125*(y)*(z)) + (0.125*(z)*(x)) + (0.125*(x)*(y)*(z)))

int main() {

  int xyz[][3] = {
    {-1,-1,-1},
    {1,-1,-1},
    {-1,1,-1},
    {1,1,-1},
    {-1,-1,1},
    {1,-1,1},
    {-1,1,1},
    {1,1,1},
    {-1,3,-1}};

  for(int i = 0;i < 9; i++) {
    std::cout<<"phi0 @ "<<i<<": "<<(phi0(xyz[i][0],xyz[i][1],xyz[i][2]))<<std::endl;
  }
  std::cout<<std::endl;

  for(int i = 0;i < 9; i++) {
    std::cout<<"phi1 @ "<<i<<": "<<(phi1(xyz[i][0],xyz[i][1],xyz[i][2]))<<std::endl;
  }
  std::cout<<std::endl;

  for(int i = 0;i < 9; i++) {
    std::cout<<"phi2 @ "<<i<<": "<<(phi2(xyz[i][0],xyz[i][1],xyz[i][2]))<<std::endl;
  }
  std::cout<<std::endl;

  for(int i = 0;i < 9; i++) {
    std::cout<<"phi3 @ "<<i<<": "<<(phi3(xyz[i][0],xyz[i][1],xyz[i][2]))<<std::endl;
  }
  std::cout<<std::endl;

  for(int i = 0;i < 9; i++) {
    std::cout<<"phi4 @ "<<i<<": "<<(phi4(xyz[i][0],xyz[i][1],xyz[i][2]))<<std::endl;
  }
  std::cout<<std::endl;
  
  for(int i = 0;i < 9; i++) {
    std::cout<<"phi5 @ "<<i<<": "<<(phi5(xyz[i][0],xyz[i][1],xyz[i][2]))<<std::endl;
  }
  std::cout<<std::endl;
  
  for(int i = 0;i < 9; i++) {
    std::cout<<"phi6 @ "<<i<<": "<<(phi6(xyz[i][0],xyz[i][1],xyz[i][2]))<<std::endl;
  }
  std::cout<<std::endl;
  
  for(int i = 0;i < 9; i++) {
    std::cout<<"phi7 @ "<<i<<": "<<(phi7(xyz[i][0],xyz[i][1],xyz[i][2]))<<std::endl;
  }
  std::cout<<std::endl;

}

