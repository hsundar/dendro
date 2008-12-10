
#include "mpi.h"
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>

float uniform();

double gaussian(double mean = 0.5, double std_dev = 0.1);

int main(int argc, char**argv) {

  MPI_Init(&argc, &argv);

  int rank;
  int npes;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);

  unsigned int numPts = atoi(argv[1]);
  double mean = atof(argv[2]);

  srand(0);

  for(unsigned int i = 0; i < (3*numPts*rank); i++) {
    int dummy = rand();
  }

  char fname[256];
  sprintf(fname, "%s%d_%d.pts" , argv[3], rank, npes);

  if(!rank) {
    std::cout<<"Generating "<<numPts<<" points (per processor) with a Log-Normal"<<
      " distribution with mean = "<<mean<<std::endl;
    std::cout<<"P0 writing output to file: "<<fname<<std::endl;
  }

  double meanForGaussian = (log(mean) - 0.5);

  double* xyz = new double[3*numPts];

  for(unsigned int i = 0; i < (3*numPts); i++) {
    xyz[i] = exp(gaussian(meanForGaussian, 1.0));
  }

  double* xyzDup = new double[3*numPts];
  for(unsigned int i = 0; i < (3*numPts); i++) {
    xyzDup[(3*numPts)-i-1] = (1.0 - xyz[i]);
  }

  unsigned int totalPts = (2*numPts);

  FILE* outfile;
  outfile = fopen(fname,"wb");
  fwrite(&totalPts, sizeof(unsigned int), 1, outfile);
  fwrite(xyz, sizeof(double), (3*numPts), outfile);
  fwrite(xyzDup, sizeof(double), (3*numPts), outfile);
  fclose(outfile);

  delete [] xyz;
  delete [] xyzDup;

  MPI_Finalize();
}

float uniform() {
  return float(rand()) / RAND_MAX; // [0,1)
}

double gaussian(double mean, double std_deviation) {
  static double t1 = 0, t2=0;
  double x1, x2, x3, r;

  // reuse previous calculations
  if(t1) {
    const double tmp = t1;
    t1 = 0;
    return mean + std_deviation * tmp;
  }
  if(t2) {
    const double tmp = t2;
    t2 = 0;
    return mean + std_deviation * tmp;
  }

  // pick randomly a point inside the unit disk
  do {
    x1 = 2 * uniform() - 1;
    x2 = 2 * uniform() - 1;
    x3 = 2 * uniform() - 1;
    r = x1 * x1 + x2 * x2 + x3*x3;
  } while(r >= 1);

  // Box-Muller transform
  r = sqrt(-2.0 * log(r) / r);

  // save for next call
  t1 = (r * x2);
  t2 = (r * x3);

  return mean + (std_deviation * r * x1);
}


