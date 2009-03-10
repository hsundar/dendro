
#include "mpi.h"
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
  MPI_Init(&argc, &argv);

  if (argc < 4) {
    std::cerr << "<exe> numAng1 numAng2 fbmR" << std::endl;
    return -1;
  }

  int numAng1 = atoi(argv[1]);
  int numAng2 = atoi(argv[2]);
  double fbmR = atof(argv[3]);

  int rank, npes;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);

  int rootP = static_cast<int>(sqrt(npes));
  assert( (rootP*rootP) == npes );

  int iIdx = rank/rootP;
  int jIdx = rank%rootP;

  double ang1Start = static_cast<double>(iIdx)*2.0*__PI__/static_cast<double>(rootP);
  double ang1End = static_cast<double>(iIdx + 1)*2.0*__PI__/static_cast<double>(rootP);

  double ang2Start = static_cast<double>(jIdx)*__PI__/static_cast<double>(rootP);
  double ang2End = static_cast<double>(jIdx + 1)*__PI__/static_cast<double>(rootP);

  //The gather is for testing only
#ifdef __DEBUG__
  double* allAng1Start = new double[npes];
  double* allAng2Start = new double[npes];

  MPI_Gather(&ang1Start, 1, MPI_DOUBLE, allAng1Start, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&ang2Start, 1, MPI_DOUBLE, allAng2Start, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if(!rank) {
    for(int i = 0; i < npes; i++) {
      std::cout<<"ang1Start["<<i<<"] = "<<allAng1Start[i]<<std::endl;
      std::cout<<"ang2Start["<<i<<"] = "<<allAng2Start[i]<<std::endl;
    }
  }

  delete [] allAng1Start;
  delete [] allAng2Start;
#endif 

  double angInc1 = (ang1End - ang1Start)/static_cast<double>(numAng1);
  double angInc2 = (ang2End - ang2Start)/static_cast<double>(numAng2);

  std::vector<double> outPts;

  for(double ang1 = ang1Start; ang1 < ang1End; ang1 += angInc1) {
    for(double ang2 = ang2Start; ang2 < ang2End; ang2 += angInc2) {
      double x = 0.5 + fbmR*cos(ang1)*sin(ang2);
      double y = 0.5 + fbmR*sin(ang1)*sin(ang2);
      double z = 0.5 + fbmR*cos(ang2);
      outPts.push_back(x);
      outPts.push_back(y);
      outPts.push_back(z);
    }
  }

  char pFile[256];
  sprintf(pFile, "fbmInp%d_%d.pts", rank, npes);
  writePtsToFile(&(*outPts.begin()), outPts.size(), pFile);

  outPts.clear();

  MPI_Finalize();
}



