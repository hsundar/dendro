
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>

float uniform();

double gaussian(double mean = 0.5, double std_dev = 0.1);

int main(int argc, char**argv) {
  if(argc < 4) {
    std::cout<<"<exe> numPts mean outFile"<<std::endl;
    exit(0);
  }
  unsigned int numPts = atoi(argv[1]);
  double mean = atof(argv[2]);

  std::cout<<"Generating "<<numPts<<" points with a Log-Normal"<<
    " distribution with mean = "<<mean<<std::endl;

  std::cout<<"Writing output to file: "<<argv[3]<<std::endl;

  double meanForGaussian = (log(mean) - 0.5);

  double* xyz = new double[3*numPts];

  srand(0);
  for(unsigned int i = 0; i < (3*numPts); i++) {
    xyz[i] = exp(gaussian(meanForGaussian, 1.0));
  }

  double maxXYZ[3];
  double minXYZ[3];

  maxXYZ[0] = xyz[0];
  minXYZ[0] = xyz[0];
  maxXYZ[1] = xyz[1];
  minXYZ[1] = xyz[1];
  maxXYZ[2] = xyz[2];
  minXYZ[2] = xyz[2];
  for(int i = 0; i < numPts; i++) {

    for(int j = 0; j < 3; j++) {
      if(xyz[3*i + j] > maxXYZ[j]) {
        maxXYZ[j] = xyz[3*i + j];
      }

      if(xyz[3*i + j] < minXYZ[j]) {
        minXYZ[j] = xyz[3*i + j];
      }

    }//end for j

  }//end for i

  //rescale to 0.0 to 1.0

  for(int i = 0; i < numPts; i++) {
    for(int j = 0; j < 3; j++) {
      xyz[3*i + j] = (xyz[3*i + j] - minXYZ[j])/((maxXYZ[j] - minXYZ[j]) + 0.01);
    }
  }

  double* xyzDup = new double[3*numPts];
  for(unsigned int i = 0; i < (3*numPts); i++) {
    xyzDup[(3*numPts)-i-1] = (0.5 - xyz[i]);
  }

  unsigned int totalPts = (2*numPts);

  maxXYZ[0] = xyz[0];
  minXYZ[0] = xyz[0];
  maxXYZ[1] = xyz[1];
  minXYZ[1] = xyz[1];
  maxXYZ[2] = xyz[2];
  minXYZ[2] = xyz[2];
  for(int i = 0; i < numPts; i++) {

    for(int j = 0; j < 3; j++) {
      if(xyz[3*i + j] > maxXYZ[j]) {
        maxXYZ[j] = xyz[3*i + j];
      }

      if(xyz[3*i + j] < minXYZ[j]) {
        minXYZ[j] = xyz[3*i + j];
      }

      if(xyzDup[3*i + j] > maxXYZ[j]) {
        maxXYZ[j] = xyzDup[3*i + j];
      }

      if(xyzDup[3*i + j] < minXYZ[j]) {
        minXYZ[j] = xyzDup[3*i + j];
      }

    }//end for j

  }//end for i

  //rescale to 0.0 to 1.0

  for(int i = 0; i < numPts; i++) {
    for(int j = 0; j < 3; j++) {
      xyz[3*i + j] = (xyz[3*i + j] - minXYZ[j])/((maxXYZ[j] - minXYZ[j]) + 0.01);
      xyzDup[3*i + j] = (xyzDup[3*i + j] - minXYZ[j])/((maxXYZ[j] - minXYZ[j]) + 0.01);
    }
  }

  FILE* outfile;
  outfile = fopen(argv[3],"wb");
  fwrite(&totalPts, sizeof(unsigned int), 1, outfile);
  fwrite(xyz, sizeof(double), (3*numPts), outfile);
  fwrite(xyzDup, sizeof(double), (3*numPts), outfile);
  fclose(outfile);

  delete [] xyz;
  delete [] xyzDup;

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


