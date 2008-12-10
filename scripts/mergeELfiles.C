#include <iostream>
#include <cstdio>
#include <cassert>

using namespace std;

void readElemListHeaderFromFile(unsigned int & maxD, unsigned int & numElems, char* filename) {
	FILE* infile = fopen(filename,"r");

	fscanf(infile, "%u", &maxD); 
	fscanf(infile, "%u", &numElems); 

	fclose(infile);
}//end function

void readElemListFromFile(unsigned int** & elemList, unsigned int & numElems, char* filename) {
	FILE* infile = fopen(filename,"r");
	unsigned int  maxD;
	fscanf(infile, "%u", &maxD); 
	fscanf(infile, "%u", &numElems); 

	typedef unsigned int* uiPtr;

	elemList = new uiPtr [numElems];
	for(unsigned int i = 0; i < numElems; i++) {
		elemList[i] = new unsigned int [5];
		for(unsigned int j = 0; j < 5; j++) {
			fscanf(infile, "%u", &(elemList[i][j]));	
		}	   	
	}

	fclose(infile);
}//end function

void writeElemListHeaderToFile( const unsigned int maxD, const unsigned int numElems, char* filename) {
	FILE*  outfile = fopen(filename,"w");

	fprintf(outfile, "%u %u\n", maxD,numElems);

	fclose(outfile);
}//end function

void appendElemListToFile(unsigned** elemList, const unsigned int numElems, char* filename) {
	FILE*  outfile = fopen(filename,"a");
	for(unsigned int i = 0; i < numElems; i++) {
		for(unsigned int j = 0; j < 5; j++) {
			fprintf(outfile, " %u ", elemList[i][j]);
		}
		fprintf(outfile, "\n");
	}
	fclose(outfile);
}//end function


int main(int argc, char **argv) {
	if (argc < 3) {
		std::cerr << "mergeELfiles InFilePrefix numProcs outFile" << std::endl;
		return -1;
	}

	int p = atoi(argv[2]);

	char fname[256];

	unsigned int totalElems = 0;
	unsigned int maxD = 0;
	// Read in the header ...
	for(int i = 0; i < p ; i++) {
		unsigned int numElems = 0;
		unsigned int tmpMaxD = 0;
		sprintf(fname, "%s_EL_%d_%d.out", argv[1], i, p);
		readElemListHeaderFromFile(tmpMaxD, numElems, fname); 

		if(!i) {
		   maxD = tmpMaxD;	
		} else {
		  assert( (maxD == tmpMaxD) || (numElems == 0) );
		}

		totalElems += numElems;
	}

	//Write the header
	writeElemListHeaderToFile(maxD, totalElems, argv[3]);

	//Read and write 

	for (int i = 0; i < p; i++) {
		sprintf(fname, "%s_EL_%d_%d.out", argv[1], i, p);

		unsigned int **elemList = NULL;
		unsigned int numElems = 0;

		readElemListFromFile(elemList, numElems, fname);

		appendElemListToFile(elemList, numElems, argv[3]);

		for(unsigned int j = 0; j < numElems; j++) {
			delete [] (elemList[j]);
		}
		delete [] elemList;

	}

	return 0;
}

