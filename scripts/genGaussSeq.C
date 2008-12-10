
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>

float uniform();

double gaussian(double mean = 0.5, double std_dev = 0.1);

int main(int argc, char**argv) {
  unsigned int numPts = atoi(argv[1]);
  double mean = atof(argv[2]);

  std::cout<<"Generating "<<numPts<<" points with a Normal"<<
    " distribution with mean = "<<mean<<std::endl;

  std::cout<<"Writing output to file: "<<argv[3]<<std::endl;

  double* xyz = new double[3*numPts];

  srand(0);
  for(unsigned int i = 0; i < (3*numPts); i++) {
    xyz[i] = gaussian(mean, 1.0);
  }

  FILE* outfile;
  outfile = fopen(argv[3],"wb");
  fwrite(&numPts, sizeof(unsigned int), 1, outfile);
  fwrite(xyz, sizeof(double), (3*numPts), outfile);
  fclose(outfile);

  delete [] xyz;

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


