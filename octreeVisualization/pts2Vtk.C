
#include <cstdio>
#include <iostream>

int main(int argc, char**argv) {
  if(argc < 2) {
    std::cerr<<"exe ptsFilePrefix"<<std::endl;
    exit(0);
  }
  char ptsFileName[256];
  char vtkFileName[256];
  sprintf(ptsFileName, "%s.pts", argv[1]);
  sprintf(vtkFileName, "%s.vtk", argv[1]);

  FILE* inFile = fopen(ptsFileName, "rb");

  unsigned int numPts;
  fread(&numPts, sizeof(unsigned int), 1, inFile);

  double* xyz = new double[3*numPts];
  fread(xyz, sizeof(double), (3*numPts), inFile);

  fclose(inFile);

  FILE * outFile = fopen(vtkFileName, "wa");

  fprintf(outFile, "# vtk DataFile Version 3.0\n"); 
  fprintf(outFile, "Octree field file\n");
  fprintf(outFile, "ASCII\n");
  fprintf(outFile, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(outFile, "POINTS %u float\n", numPts);

  for(unsigned int i = 0; i < numPts; i++) {
    fprintf(outFile, "%f %f %f\n", static_cast<float>(xyz[3*i]),
        static_cast<float>(xyz[(3*i)+1]), static_cast<float>(xyz[(3*i)+2]));
  }

  fprintf(outFile, "\nCELLS %u %u\n", numPts, (2*numPts));

  for(unsigned int i = 0; i < numPts; i++) {
    fprintf(outFile, "1 %u\n", i);
  }

  fprintf(outFile, "\nCELL_TYPES %u\n", numPts);

  for(unsigned int i = 0; i < numPts; i++) {
    fprintf(outFile, "1\n");
  }

  fprintf(outFile, "\nCELL_DATA %u\n", numPts);
  fprintf(outFile, "SCALARS scalars float\n");
  fprintf(outFile, "LOOKUP_TABLE default\n");

  for(unsigned int i = 0; i < numPts; i++) {
    fprintf(outFile, "1.0\n");
  }

  delete [] xyz;

  fclose(outFile);
}


