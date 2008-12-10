#include <iostream>
#include <cstdio>

using namespace std;

void readHeaderFromFile( unsigned int & numPts, char* filename) {
	FILE* infile = fopen(filename,"rb");

	fread(&numPts, sizeof(unsigned int),1,infile);

	fclose(infile);
}//end function

void readPtsFromFile( double* & pts, unsigned int* ptsLen, char* filename) {
	FILE* infile = fopen(filename,"rb");
	unsigned int temp;
	fread(&temp, sizeof(unsigned int),1,infile);
	*ptsLen = 3*temp;
	pts = new double[3*temp];

	fread(pts, sizeof(double),3*temp,infile);

	fclose(infile);
}//end function

void writeHeaderToFile( const unsigned int numPts, char* filename) {
	FILE*  outfile = fopen(filename,"wb");
	if(numPts > 0) {
		fwrite(&numPts, sizeof(unsigned int),1,outfile);
		std::cout << " Total points: "<<numPts << std::endl;
	}
	fclose(outfile);
}//end function

void appendPtsToFile( const double* pts, const unsigned int ptsLen, char* filename) {
	FILE*  outfile = fopen(filename,"ab");
	if(ptsLen >0) {
		fwrite(pts, sizeof(double),ptsLen,outfile);
	}
	fclose(outfile);
}//end function


int main(int argc, char **argv) {
	if (argc < 3) {
		std::cerr << "mergePoints InFilePrefix numProcs outFile" << std::endl;
		return -1;
	}

	int p = atoi(argv[2]);

	char fname[256];

	unsigned int totalPts = 0;
	// Read in the header ...
	for(int i = 0; i < p ; i++) {
		unsigned int numPts = 0;
		sprintf(fname, "%s%d_%d.pts", argv[1], i, p);
		readHeaderFromFile( numPts, fname); 
		totalPts += numPts;
	}

	//Write the header
	writeHeaderToFile( totalPts, argv[3]);

	//Read and write points

	for (int i = 0; i < p; i++) {
		sprintf(fname, "%s%d_%d.pts", argv[1], i, p);

		double *points = NULL;
		unsigned int ptsLen;

		// Read in the points ...
		readPtsFromFile(points, &ptsLen, fname);

		appendPtsToFile( points, ptsLen, argv[3]);

		delete [] points;

	}

	return 0;
}

