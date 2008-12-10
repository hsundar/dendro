
#include <cstdio>
#include <cmath>
#include <iostream>

#define MY_PI 3.14159265

int main(int argc, char**argv) {
  if(argc < 3) {
    std::cerr<<"exe outFilePrefix numNodesPerDim"<<std::endl;
    exit(0);
  }
  char vtkFileName[256];
  sprintf(vtkFileName, "%s.vtk", argv[1]);

  int numNodesPerDim = atoi(argv[2]);
  double h = (1.0/static_cast<double>(numNodesPerDim));

  FILE * outFile = fopen(vtkFileName, "wa");

  fprintf(outFile, "# vtk DataFile Version 3.0\n"); 
  fprintf(outFile, "Material Property File\n");
  fprintf(outFile, "ASCII\n");
  fprintf(outFile, "DATASET STRUCTURED_POINTS\n");
  fprintf(outFile, "DIMENSIONS %d %d %d\n", numNodesPerDim,
      numNodesPerDim, numNodesPerDim);
  fprintf(outFile, "ORIGIN 0.0 0.0 0.0 \n");
  fprintf(outFile, "SPACING %lf %lf %lf\n", h, h, h);
  fprintf(outFile, "POINT_DATA %d\n",(numNodesPerDim*numNodesPerDim*numNodesPerDim));
  fprintf(outFile, "SCALARS matProp double 1\n");
  fprintf(outFile, "LOOKUP_TABLE default\n");

  double posFac = 2.0*MY_PI*h;
  for(int k = 0; k < numNodesPerDim; k++) {
    for(int j = 0; j < numNodesPerDim; j++) {
      for(int i = 0; i < numNodesPerDim; i++) {
        double xPos = static_cast<double>(i)*posFac;
        double yPos = static_cast<double>(j)*posFac;
        double zPos = static_cast<double>(k)*posFac;
        double val = 1.0 + 1.0e6*((cos(xPos)*cos(xPos)) +
            (cos(yPos)*cos(yPos)) +
            (cos(zPos)*cos(zPos)));
        fprintf(outFile, "%lf \n", val);
      }//end for i
    }//end for j
  }//end for k

  fclose(outFile);
}


