
/**
@file binUtils.C
@brief A set of functions for fast binary operations
@author Rahul S. Sampath, rahul.sampath@gmail.com
@author Hari Sundar, hsundar@gmail.com
*/

#include <vector>

namespace binOp {

  unsigned int binLength(unsigned int num) {
    unsigned int len =1;
    while(num > 1) {
        num = (num>>1);
        len++;
    }
    return len;
  }//end function

  unsigned int fastLog2(unsigned int num) {
    unsigned int len = 0;
    while(num > 0) {
        num = (num>>1);
        len++;
    }
    return len;
  }//end function

  int toBin(unsigned int num, unsigned int binLen,  std::vector<bool>& numBin) {
    numBin = std::vector<bool>(binLen);
    for(unsigned int i=0;i<binLen;i++) {
      numBin[i]=0;
    }//end for
    unsigned int pos = binLen -1;
    while(num > 0) {
      numBin[pos] = (num%2);
      num = num/2;  
      pos--;
    }  //end while  
    return 1;
  }//end function

  unsigned int binToDec(unsigned int* numBin, unsigned int binLen) {
    unsigned int res = 0;
    for(unsigned int i = 0; i< binLen; i++) {
      res = (2*res) + numBin[i];
    }
    return res;
  }//end function


bool isPowerOfTwo(unsigned int n) {
  return !(n & (n - 1)) && n;
}

// compute the next highest power of 2 of 32-bit v
int getNextHighestPowerOfTwo(unsigned int n) {
  unsigned int v = n;
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return v;
}

// compute the prev highest power of 2 of 32-bit v
int getPrevHighestPowerOfTwo(unsigned int n) {
  unsigned int v = n;
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return  v >> 1;
}

}//end namespace

