#include <iostream>
#include <cstdio>
#include <vector>
#include <cassert>

using namespace std;

void readPtsFromFile( double* & pts, unsigned int* ptsLen, char* filename) {
  FILE* infile = fopen(filename,"r");
  unsigned int temp, temp2;
  fread(&temp, sizeof(unsigned int),1,infile);
  *ptsLen = 3*temp;
  pts = new double[3*temp];
  std::cout << temp << " points" << std::endl;

  fread(pts, sizeof(double),3*temp,infile);

  fclose(infile);
}//end function

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
  if (argc < 3) {
    std::cerr << "<exe> inFile halfWidth(<= 0.16) outFile" << std::endl;
    return -1;
  }

  double *points;
  std::vector<double> outPts;
  unsigned int numPts;
  double halfWidth  = atof(argv[2]);
  assert(halfWidth > 0);
  assert(halfWidth <= 0.16);
  double beg = 0.5 - halfWidth;
  double end = 0.5 + halfWidth;

  // Read in the points ...
  readPtsFromFile(points, &numPts, argv[1]);

  numPts /= 3;

  // Now modify and write ...
    
  for (int i=0; i < numPts; i++) {   
    double px = points[3*i];
    double py = points[3*i + 1];
    double pz = points[3*i + 2];	
    if(px > beg && px < end && py > beg && py < end && pz > beg && pz < end ) {
	//-1,-1,-1 
	outPts.push_back(px-0.33);
	outPts.push_back(py-0.33);
	outPts.push_back(pz-0.33);

	//+1,+1,+1
	outPts.push_back(px+0.33);
	outPts.push_back(py+0.33);
	outPts.push_back(pz+0.33);
	
    }
  }

  delete [] points;

  std::cout<<"Writing "<<(outPts.size()/3)<<" points."<<std::endl;

  writePtsToFile(&(*outPts.begin()), outPts.size(), argv[3]);

  outPts.clear();

  return 0;
}

