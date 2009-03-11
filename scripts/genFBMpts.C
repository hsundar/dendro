
#include <iostream>
#include <cstdio>
#include <vector>
#include <cassert>
#include <cmath>

#define __PI__ 3.14159265

void writePtsToFile( const double* pts, const unsigned int ptsLen, char* filename) {
  FILE*  outfile = fopen(filename,"w");
  if(ptsLen >0) {
    unsigned int numPts = ptsLen/3;
    fwrite(&numPts, sizeof(unsigned int),1,outfile);
    fwrite(pts, sizeof(double),ptsLen,outfile);
  }
  fclose(outfile);
}//end function

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "<exe> numAng1 numAng2 fbmR" << std::endl;
    return -1;
  }

  int numAng1 = atoi(argv[1]);
  int numAng2 = atoi(argv[2]);
  double fbmR = atof(argv[3]);

  double angInc1 = 2*__PI__/static_cast<double>(numAng1);
  double angInc2 = __PI__/static_cast<double>(numAng2);

  std::vector<double> outPts;

  for(double ang1 = 0; ang1 < (2*__PI__); ang1 += angInc1) {
    for(double ang2 = 0; ang2 < __PI__; ang2 += angInc2) {
      double x = 0.5 + fbmR*cos(ang1)*sin(ang2);
      double y = 0.5 + fbmR*sin(ang1)*sin(ang2);
      double z = 0.5 + fbmR*cos(ang2);
      outPts.push_back(x);
      outPts.push_back(y);
      outPts.push_back(z);
    }
  }

  std::cout<<"Writing "<<(outPts.size()/3)<<" points."<<std::endl;

  writePtsToFile(&(*outPts.begin()), outPts.size(), "fbmInp.pts");

  outPts.clear();

  return 0;
}



