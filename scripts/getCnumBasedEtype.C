
int createChildNumBasedVtxLUT(unsigned char**&lut ) ;
unsigned char getEtypeInv(unsigned char inMask) ;
int mapFlagsInv(unsigned char **lut, int childNum,unsigned char & mask);
int destroyChildNumBasedVtxLUT(unsigned char**&a) ;

int int2str(int n,char*s) ;
char int2char(int d) ;
enum exhaustiveElemType {
  //Order: 654321
  //YZ ZX Z XY Y X
  ET_N = 0,
  ET_Y = 2,
  ET_X = 1,
  ET_XY = 3,
  ET_Z = 8,
  ET_ZY = 10,
  ET_ZX = 9,
  ET_ZXY = 11,
  ET_XY_XY = 7,
  ET_XY_ZXY = 15,
  ET_YZ_ZY = 42,
  ET_YZ_ZXY = 43,
  ET_YZ_XY_ZXY = 47,
  ET_ZX_ZX = 25,
  ET_ZX_ZXY = 27,
  ET_ZX_XY_ZXY = 31,
  ET_ZX_YZ_ZXY = 59,
  ET_ZX_YZ_XY_ZXY = 63
};

int main() {
  unsigned char**lut;
  createChildNumBasedVtxLUT(lut);
  for(int childNum =0; childNum<8; childNum++) {
    unsigned char trueMask;
    char filename[20];
    char tmpStr[20];
    strcpy(filename,"cNumEtype_\0");
    int2str(childNum,tmpStr);
    strcat(filename,tmpStr);
    strcat(filename,".h\0");

    std::ofstream outfile (filename);

    trueMask = getEtypeInv( ET_Y ) ;
    mapFlagsInv(lut, childNum, trueMask) ;
    outfile<<" = "<<((int)trueMask)<<" ,\n";

    trueMask = getEtypeInv( ET_X ) ;
    mapFlagsInv(lut, childNum, trueMask) ;
    outfile<<" = "<<((int)trueMask)<<" ,\n";

    trueMask = getEtypeInv( ET_XY ) ;
    mapFlagsInv(lut, childNum, trueMask) ;
    outfile<<" = "<<((int)trueMask)<<" ,\n";

    trueMask = getEtypeInv( ET_Z ) ;
    mapFlagsInv(lut, childNum, trueMask) ;
    outfile<<" = "<<((int)trueMask)<<" ,\n";

    trueMask = getEtypeInv( ET_ZY ) ;
    mapFlagsInv(lut, childNum, trueMask) ;
    outfile<<" = "<<((int)trueMask)<<" ,\n";

    trueMask = getEtypeInv( ET_ZX ) ;
    mapFlagsInv(lut, childNum, trueMask) ;
    outfile<<" = "<<((int)trueMask)<<" ,\n";

    trueMask = getEtypeInv( ET_ZXY ) ;
    mapFlagsInv(lut, childNum, trueMask) ;
    outfile<<" = "<<((int)trueMask)<<" ,\n";

    trueMask = getEtypeInv( ET_XY_XY ) ;
    mapFlagsInv(lut, childNum, trueMask) ;
    outfile<<" = "<<((int)trueMask)<<" ,\n";

    trueMask = getEtypeInv( ET_XY_ZXY ) ;
    mapFlagsInv(lut, childNum, trueMask) ;
    outfile<<" = "<<((int)trueMask)<<" ,\n";

    trueMask = getEtypeInv( ET_YZ_ZY ) ;
    mapFlagsInv(lut, childNum, trueMask) ;
    outfile<<" = "<<((int)trueMask)<<" ,\n";

    trueMask = getEtypeInv( ET_YZ_ZXY ) ;
    mapFlagsInv(lut, childNum, trueMask) ;
    outfile<<" = "<<((int)trueMask)<<" ,\n";

    trueMask = getEtypeInv( ET_YZ_XY_ZXY ) ;
    mapFlagsInv(lut, childNum, trueMask) ;
    outfile<<" = "<<((int)trueMask)<<" ,\n";

    trueMask = getEtypeInv( ET_ZX_ZX ) ;
    mapFlagsInv(lut, childNum, trueMask) ;
    outfile<<" = "<<((int)trueMask)<<" ,\n";

    trueMask = getEtypeInv( ET_ZX_ZXY ) ;
    mapFlagsInv(lut, childNum, trueMask) ;
    outfile<<" = "<<((int)trueMask)<<" ,\n";

    trueMask = getEtypeInv( ET_ZX_XY_ZXY ) ;
    mapFlagsInv(lut, childNum, trueMask) ;
    outfile<<" = "<<((int)trueMask)<<" ,\n";

    trueMask = getEtypeInv( ET_ZX_YZ_ZXY ) ;
    mapFlagsInv(lut, childNum, trueMask) ;
    outfile<<" = "<<((int)trueMask)<<" ,\n";

    trueMask = getEtypeInv( ET_ZX_YZ_XY_ZXY ) ;
    mapFlagsInv(lut, childNum, trueMask) ;
    outfile<<" = "<<((int)trueMask)<<" ,\n";

    outfile.close();
  }
  destroyChildNumBasedVtxLUT(lut);
}

int int2str(int n,char*s) {
  int tmpd[20];
  int i=0;
  int j;   
  if (n==0) {
    strcpy(s,"0\0");
  } else {
    while (n>0) {
      tmpd[i]= (n%10);
      n= (int)(n/10);
      i++;
    }
    for (j=i-1;j>=0;j--) {
      s[i-j-1]=int2char(tmpd[j]);
    }
    s[i]='\0';
  }
  return 1;
}//end function

char int2char(int d) {
  switch (d) {
    case 0: return '0';
    case 1: return '1';
    case 2: return '2';
    case 3: return '3';
    case 4: return '4';
    case 5: return '5';
    case 6: return '6';
    case 7: return '7';
    case 8: return '8';
    case 9: return '9';
    default: return '\0';
  }
}


//inMask: XX654321
//OutMask: 76543210
//0 and 7 will always be 0.
unsigned char getEtypeInv(unsigned char inMask) {
  unsigned char type = (inMask<<1);
  return type;
}//end fn2

int mapFlagsInv(unsigned char **lut, int childNum,unsigned char & mask) {
  unsigned char tmpFlags = 0;
  for(int i=0;i<8;i++) {
    tmpFlags = ( tmpFlags | ( ( (1 << i) & mask ) ? (1 << (lut[childNum][i])) : 0 ) );
  }
  mask = tmpFlags;
  return(0);
}//end function

int createChildNumBasedVtxLUT(unsigned char**&lut ) {
  //Note: It is not symmetric.
  unsigned char tmp[8][8]={
    {0,1,2,3,4,5,6,7},
    {1,3,0,2,5,7,4,6},
    {2,0,3,1,6,4,7,5},
    {3,2,1,0,7,6,5,4},
    {4,5,0,1,6,7,2,3},
    {5,7,1,3,4,6,0,2},
    {6,4,2,0,7,5,3,1},
    {7,6,3,2,5,4,1,0}
  };

  //Is Stored in  ROW_MAJOR Format.  
  typedef unsigned char* charPtr;
  lut = new charPtr[8];
  for(int i=0;i<8;i++) {
    lut[i] = new unsigned char[8]; 
    for(int j=0;j<8;j++) {
      lut[i][j] = tmp[i][j];
    }
  }
  return(0);
}//end function.

int destroyChildNumBasedVtxLUT(unsigned char**&a) {
  for(int i=0;i<8;i++) {
    delete [] a[i];
    a[i] = NULL;
  }//end i
  delete []a;
  a=NULL;
  return(0);
}//end fn2

