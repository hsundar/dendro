static char help[] = " Visualization routine";

#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "oct.h"
#include "oda.h"
#include "omg.h"
#include "par.h"
#include "Point.h"
#include "parUtils.h"
#include "octUtils.h"
#include "TreeNode.h"
#include "dendro.h"


namespace par {
 template <>
   class Mpi_datatype< Point > {
     public:
     static MPI_Datatype value()
     {
       static bool         first = true;
       static MPI_Datatype datatype;

       if (first)
       {
         first = false;
         MPI_Type_contiguous(3, Mpi_datatype<double>::value(), &datatype);
         MPI_Type_commit(&datatype);
       }

       return datatype;
     }
   };
}

int main(int argc, char **argv)
{
  int rank, size;
  double sliceAt = 0.5;

  MPI_Init(&argc,&argv);

  // get the size and rank
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  FILE* fpNode;
  FILE* fpSolution;

  
  char fileNameNode[50];
  char fileNameSolution[50];

  // local files to read points and solution
  sprintf(fileNameNode,"nodes_%d.dat",rank);
  sprintf(fileNameSolution,"solution_%d.dat",rank);
  fpNode = fopen(fileNameNode,"r");
  fpSolution = fopen(fileNameSolution,"r");
  
  std::vector<float> nodeX;
  std::vector<float> nodeY;
  std::vector<float> nodeZ;
  std::vector<float> solution;
  // read nodes files and solution files
  int numNodes;
  int numSolution;
  fscanf(fpNode,"%d",&numNodes);
  fscanf(fpSolution,"%d",&numSolution);

  // assert that they are equal
  assert(numNodes == numSolution);
  if(numNodes != numSolution)
	 {
		std::cerr << rank << ":: sizes of the nodes and the solution are not compatible " << std::endl;
	 }

  // Get the Slice information
  for(int i = 0; i < numNodes; i++)
	 {
		float x,y,z;
		float value;
		fscanf(fpNode,"%f",&x);
		fscanf(fpNode,"%f",&y);
		fscanf(fpNode,"%f",&z);
		fscanf(fpSolution,"%f",&value);
		if( z < sliceAt + 0.05 && z > sliceAt - 0.05)
		  {
			 //			 Point pt(x,y,z);
			 nodeX.push_back(x);
			 nodeY.push_back(y);
			 nodeZ.push_back(z);
			 solution.push_back(value);
		  }
	 }
  fclose(fpNode);
  fclose(fpSolution);


  // Merge routines
  float *localSolution;
  float *localPointX;
  float *localPointY;
  float *localPointZ;
  int localSize = nodeX.size();
  DendroIntL globalSize;
  int scanSize;

  // local arrays
  //localSolution = new float[localSize];
  //localPointX = new float[localSize];
  //  localPointY = new float[localSize];
  //  localPointZ = new float[localSize];
  localSolution = (&(*solution.begin()));
  localPointX   = (&(*nodeX.begin()));
  localPointY   = (&(*nodeY.begin()));
  localPointZ   = (&(*nodeZ.begin()));

  int *recvCnts;
  recvCnts = new int[size];
  MPI_Gather(&localSize, 1, MPI_INT, recvCnts, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Convert vectors to int* and float*
  /*  for(int i = 0; i < localSize; i++)
	 {
		localSolution[i] = solution[i];
		localPointX[i]  = nodeX[i];
		localPointY[i]  = nodeY[i];
		localPointZ[i]  = nodeZ[i];
	 }
  nodeX.clear();
  nodeY.clear();
  nodeZ.clear();
  solution.clear();  */

  // get the global size
  DendroIntL localSizeCopy = localSize;
  par::Mpi_Reduce<DendroIntL>(&localSizeCopy, &globalSize, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  if(!rank)
	 {
	   std::cout << " Global size = " << globalSize << std::endl;
		for(int i = 0; i < size; i++)
		  {
			 std::cout << recvCnts[i] << " ";
		  }
		std::cout << std::endl;
	 }
  float *solvec;
  float *solpointX;
  float *solpointY;
  float *solpointZ;
  // allocate memory for points and solution
  if(!rank){
	 solvec = new float[globalSize];
	 solpointX = new float[globalSize];
	 solpointY = new float[globalSize];
	 solpointZ = new float[globalSize];
  }

  // evaluate displacements
  int *displs;
  if(!rank){
	 displs = new int[size];
	 for(int i = 0; i  <size; i++)
		{
		  displs[i] = 0;
		  for(int j = 0; j < i; j++)
			 {
				displs[i] += recvCnts[j];
			 }
		  std::cout << rank << " :: " << displs[i] << std::endl;
		}
  }

  // get the solution
  MPI_Gatherv(localSolution,localSize,MPI_FLOAT,solvec,recvCnts,displs,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Gatherv(localPointX,localSize,MPI_FLOAT,solpointX,recvCnts,displs,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Gatherv(localPointY,localSize,MPI_FLOAT,solpointY,recvCnts,displs,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Gatherv(localPointZ,localSize,MPI_FLOAT,solpointZ,recvCnts,displs,MPI_FLOAT,0,MPI_COMM_WORLD);

  if(!rank){
	 char fileNameVTK[50];
	 sprintf(fileNameVTK,"solAtSlice_%0.1f.vtk",sliceAt);
	 FILE* outfile;
	 outfile = fopen(fileNameVTK,"w");
	 fprintf(outfile,"# vtk DataFile Version 3.0\n");
	 fprintf(outfile,"Octree field file\n");
	 fprintf(outfile,"ASCII\n");
	 fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");
	 fprintf(outfile,"POINTS %d float\n",globalSize);
	 
	 for(int i = 0; i < globalSize; i++)
		{
		  float fx, fy, fz;
		  fx = solpointX[i]; 		  fy = solpointY[i]; 		  fz = solpointZ[i]; 
		  fprintf(outfile,"%f %f %f\n",fx,fy,fz);
		}
	 fprintf(outfile,"\n CELLS %d %d\n",globalSize,globalSize*2);
	 for(int i = 0; i < globalSize; i++)
		{
		  fprintf(outfile,"1 ");
		  fprintf(outfile,"%d\n",i);
		}
	 fprintf(outfile,"\nCELL_TYPES %d\n",globalSize);

	 for(int i = 0; i < globalSize; i++)
		{
		int eleven = 11;
		fprintf(outfile,"1 \n");
	 }
	 fprintf(outfile,"\nCELL_DATA %d\n",globalSize);
	 fprintf(outfile,"SCALARS scalars float\n");
	 fprintf(outfile,"LOOKUP_TABLE default\n");


	 for (int i =0; i< globalSize; i++) {
		float v = solvec[i];	
		fprintf(outfile,"%f \n",v);
	 }

	 fclose(outfile);
  }

  nodeX.clear();
  nodeY.clear();
  nodeZ.clear();
  solution.clear();
  MPI_Finalize();
}

