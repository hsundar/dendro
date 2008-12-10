
#include <iostream>
#include <cstdio>
#include <vector>
#include <cassert>
#include <cmath>

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define POWFIVE(x) (CUBE(x)*SQR(x))

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
    std::cerr << "<exe> outFile hx hy" << std::endl;
    return -1;
  }

  double hx = atof(argv[2]);
  double hy = atof(argv[3]);

  std::vector<double> outPts;

  for(double x = -3.0; x < 3.0; x += hx) {
    for(double y = -3.0; y < 3.0; y += hy) {
      double z = (3.0*SQR(1.0-x)*exp(-(SQR(x) + SQR(y + 1.0)))) -
        (10.0*((x/5.0) - CUBE(x) - POWFIVE(y))*exp(-(SQR(x) + SQR(y)))) - 
        ((1.0/3.0)*exp(-(SQR(x + 1.0) + SQR(y))));
      outPts.push_back(x);
      outPts.push_back(y);
      outPts.push_back(z);
    }
  }

  double maxXYZ[3];
  double minXYZ[3];
  for(int j = 0 ; j < 3; j++) {
    maxXYZ[j] = outPts[j];
    minXYZ[j] = outPts[j];
  }

  unsigned int numPts = outPts.size()/3;

  for(unsigned int i = 0; i < numPts; i++) {
    for(int j = 0; j < 3; j++) {
      if(outPts[3*i + j] > maxXYZ[j]) {
        maxXYZ[j] = outPts[3*i + j];
      }
      if(outPts[3*i + j] < minXYZ[j]) {
        minXYZ[j] = outPts[3*i + j];
      }
    }
  }

  for(unsigned int i = 0; i < numPts; i++) {
    for(int j = 0; j < 3; j++) {
      outPts[3*i +j] = (outPts[3*i + j] - minXYZ[j])/(maxXYZ[j] - minXYZ[j] + 0.001);
    }
  }

  std::cout<<"Writing "<<numPts<<" points."<<std::endl;

  writePtsToFile(&(*outPts.begin()), outPts.size(), argv[1]);

  outPts.clear();

  return 0;
}



