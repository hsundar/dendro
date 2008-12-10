
#include<iostream>
#include<vector>
#include<cstdio>
#include<cstdlib>

int writePtsToFile(char* filename, std::vector<double>& pts);

int readNodesFromFile (char* filename, std::vector<unsigned int> & octants,
		unsigned int & dim, unsigned int & maxDepth );

int main(int argc, char** argv) {

	if(argc < 3) {
		std::cout<<"<exe> octFile ptsFile"<<std::endl;
		exit(0);
	}

	unsigned int dim, maxDepth;
	std::vector<unsigned int> octants;

	readNodesFromFile(argv[1], octants, dim, maxDepth);

	std::vector<double> pts;
	unsigned int numPts = (octants.size())/4;
	pts.resize(3*numPts);

	for(int i = 0; i < numPts; i++) {
		unsigned int lev = octants[4*i + 3];
		double h = static_cast<double>(1u << (maxDepth - lev))/static_cast<double>(1u << maxDepth);
		for(int j = 0; j < 3; j++) {
			unsigned int coord = octants[4*i + j];
			double val = static_cast<double>(coord)/static_cast<double>(1u << maxDepth);
			pts[3*i + j] = (val + (0.5*h));
		}
	}//end for i 

	writePtsToFile(argv[2], pts);

}

int writePtsToFile(char* filename, std::vector<double>& pts) {
	FILE* outfile = fopen(filename,"wb");
	unsigned int ptsLen = pts.size();
	double * ptsTemp = NULL;
	if(!pts.empty()) {
		ptsTemp = (&(*(pts.begin())));
	}

	if (ptsLen >0) {
		unsigned int numPts = ptsLen/3;
		fwrite(&numPts,sizeof(unsigned int),1,outfile);
		fwrite(ptsTemp, sizeof(double),ptsLen,outfile);
	}
	fclose(outfile);
	return 1;
}//end function

int readNodesFromFile (char* filename, std::vector<unsigned int> & octants,
		unsigned int & dim, unsigned int & maxDepth ) {
	FILE* infile = fopen(filename,"r");
	unsigned int numNode;
	fscanf(infile,"%u", &dim);
	fscanf(infile,"%u", &maxDepth); 
	fscanf(infile,"%u", &numNode);
	octants.resize(4*numNode) ;   

	for (unsigned int i = 0; i < octants.size(); i++) {
		fscanf(infile,"%u", &(octants[i]));
	}
	fclose(infile);
	return 1;
}//end function


