#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;

void readPtsFromFile( double* & pts, unsigned int* ptsLen, char* filename) {
  FILE* infile = fopen(filename,"rb");
  unsigned int temp;
  fread(&temp, sizeof(unsigned int),1,infile);
  *ptsLen = 3*temp;
  pts = new double[3*temp];
  std::cout << temp << " points" << std::endl;

  fread(pts, sizeof(double),3*temp,infile);

  fclose(infile);
}//end function

void writePtsToFile( const double* pts, const unsigned int ptsLen, char* filename) {
 FILE*  outfile = fopen(filename,"wb");
  if(ptsLen >0) {
    unsigned int numPts = ptsLen/3;
    fwrite(&numPts, sizeof(unsigned int),1,outfile);
    fwrite(pts, sizeof(double),ptsLen,outfile);
  }
  fclose(outfile);
}//end function


int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << argv[0] << " infile numProcs outFilePrefix" << std::endl;
    return -1;
  }

  double *points;
  unsigned int numPts;
  unsigned int p = atoi(argv[2]);

  // Read in the points ...
  readPtsFromFile(points, &numPts, argv[1]);

  numPts /= 3;
  // Now split and write ...
  unsigned int n = numPts/p;
  unsigned int n1;

  for (int i=0; i<p; i++) {
    char fname[256];
    sprintf(fname, "%s%d_%d.pts", argv[3], i, p);
    
    if ( i == p-1 )
      n1 = 3*(numPts - n*i);
    else 
      n1 = 3*n;

    writePtsToFile( points + 3*i*n, n1, fname);
  }
  return 0;
}

