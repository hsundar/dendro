#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;
int dim=2;

void readPtsFromFile( double* & pts, unsigned int* ptsLen, char* filename) {
  FILE* infile = fopen(filename,"rb");
  unsigned int temp;
  fread(&temp, sizeof(unsigned int),1,infile);
  *ptsLen = dim*temp;
  pts = new double[dim*temp];
  std::cout << temp << " points" << std::endl;

  fread(pts, sizeof(double),dim*temp,infile);

  fclose(infile);
}//end function

void writePtsToFile( const double* pts, const unsigned int ptsLen, char* filename) {
 FILE*  outfile = fopen(filename,"wb");
  if(ptsLen >0) {
    unsigned int numPts = ptsLen/dim;
    std::cout << "numpts: " << numPts << std::endl;
    fwrite(&numPts, sizeof(unsigned int),1,outfile);
    fwrite(pts, sizeof(double),ptsLen,outfile);
  }
  fclose(outfile);
}//end function


int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << argv[0] << " infile numProcs outFilePrefix dimention" << std::endl;
    return -1;
  }

  double *points;
  unsigned int numPts;
  unsigned int p = atoi(argv[2]);
  dim=atoi(argv[4]);
  // Read in the points ...
  readPtsFromFile(points, &numPts, argv[1]);

  numPts /= dim;
  // Now split and write ...
  unsigned int n = numPts/p;
  unsigned int n1;

  for (int i=0; i<p; i++) {
    char fname[256];
    sprintf(fname, "%s%d_%d.pts", argv[3], i, p);
    
    if ( i == p-1 )
      n1 = dim*(numPts - n*i);
    else 
      n1 = dim*n;

    writePtsToFile( points + dim*i*n, n1, fname);
  }
  return 0;
}

