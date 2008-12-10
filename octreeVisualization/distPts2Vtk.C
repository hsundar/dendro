
#include <cstdio>
#include <vector>
#include <iostream>

int main(int argc, char**argv) {

	if(argc < 3) {
		std::cerr<<"exe ptsFilePrefix numProcs <suffix>"<<std::endl;
		exit(0);
	}

	int numProcs = atoi(argv[2]);

	int totalPts = 0; 
	std::vector<unsigned int> ptsSizes(numProcs);

	for(int pId = 0; pId < numProcs; pId++) {

		char ptsFileName[256];
		if(argc > 3) {
			sprintf(ptsFileName, "%s%d_%d%s.pts", argv[1], pId, numProcs, argv[3]);
		} else {
			sprintf(ptsFileName, "%s%d_%d.pts", argv[1], pId, numProcs);
		}

		FILE* inFile = fopen(ptsFileName, "rb");

		int numPts;
		fread(&numPts, sizeof(unsigned int), 1, inFile);

		totalPts += numPts;

		ptsSizes[pId] = numPts;

		fclose(inFile);

	}//end for pId

	char vtkFileName[256];
	sprintf(vtkFileName, "%s.vtk", argv[1]);

	FILE * outFile = fopen(vtkFileName, "wa");

	fprintf(outFile, "# vtk DataFile Version 3.0\n"); 
	fprintf(outFile, "Octree field file\n");
	fprintf(outFile, "ASCII\n");
	fprintf(outFile, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(outFile, "POINTS %u float\n", totalPts);

	for(int pId = 0; pId < numProcs; pId++) {

		char ptsFileName[256];
		if(argc > 3) {
			sprintf(ptsFileName, "%s%d_%d%s.pts", argv[1], pId, numProcs, argv[3]);
		} else {
			sprintf(ptsFileName, "%s%d_%d.pts", argv[1], pId, numProcs);
		}


		FILE* inFile = fopen(ptsFileName, "rb");

		unsigned int numPts;
		fread(&numPts, sizeof(unsigned int), 1, inFile);

		double* xyz = new double[3*numPts];
		fread(xyz, sizeof(double), (3*numPts), inFile);

		for(unsigned int i = 0; i < numPts; i++) {
			fprintf(outFile, "%f %f %f\n", static_cast<float>(xyz[3*i]),
					static_cast<float>(xyz[(3*i)+1]), static_cast<float>(xyz[(3*i)+2]));
		}//end for i

		delete [] xyz;

		fclose(inFile);

	}//end for pId

	fprintf(outFile, "\nCELLS %u %u\n", totalPts, (2*totalPts));

	for(unsigned int i = 0; i < totalPts; i++) {
		fprintf(outFile, "1 %u\n", i);
	}

	fprintf(outFile, "\nCELL_TYPES %u\n", totalPts);

	for(unsigned int i = 0; i < totalPts; i++) {
		fprintf(outFile, "1\n");
	}

	fprintf(outFile, "\nCELL_DATA %u\n", totalPts);
	fprintf(outFile, "SCALARS scalars float\n");
	fprintf(outFile, "LOOKUP_TABLE default\n");

	for(unsigned int pId = 0; pId < numProcs; pId++) {
		for(unsigned int i = 0; i < ptsSizes[pId]; i++) {
			fprintf(outFile, " %f \n", static_cast<float>(pId + 1));
		}//end for i
	}//end for pId

	fclose(outFile);
}


