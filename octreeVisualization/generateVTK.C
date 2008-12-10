
#include <stdio.h>
#include <iostream>
#include <fstream>

int generateVTKUnstructuredASCII(char* fileNameIn, char* fileNameOut);

int main(int argc, char **argv)
{
  char outFileName[50];
  char inFileName[50];
  for(int i = 1; i < argc; i++){
	 strcpy(inFileName,argv[i]);
	 strcpy(outFileName,argv[i]);
	 strcat(inFileName,".ot");
	 strcat(outFileName,".vtk");
	 generateVTKUnstructuredASCII(inFileName,outFileName);
  }
  
}

int generateVTKUnstructuredASCII(char* fileNameIn, char* fileNameOut)
{
  
  FILE* infile  = fopen(fileNameIn,"r");
  FILE* outfile = fopen(fileNameOut,"w");
  FILE* outfileV = fopen("Vfile","w");

  unsigned int numNode;
  unsigned int dim, maxDepth;

  float coord[8][3] = {
	 {0.0,0.0,0.0},
	 {1.0,0.0,0.0},
	 {0.0,1.0,0.0},
	 {1.0,1.0,0.0},
	 {0.0,0.0,1.0},
	 {1.0,0.0,1.0},
	 {0.0,1.0,1.0},
	 {1.0,1.0,1.0}
  };
  fscanf(infile,"%u",&dim);
  fscanf(infile,"%u",&maxDepth); 
  fscanf(infile,"%u",&numNode);
  
  fprintf(outfile,"# vtk DataFile Version 3.0\n");
  fprintf(outfile,"Octree field file\n");
  fprintf(outfile,"ASCII\n");
  fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(outfile,"POINTS %d float\n",numNode*8);

  
  for (unsigned int i =0; i< numNode; i++) {
	 unsigned int x,y,z,d;
	 float fx, fy, fz,hx;
	 fscanf(infile,"%u",&x);
	 fscanf(infile,"%u",&y);
	 fscanf(infile,"%u",&z);
	 fscanf(infile,"%u",&d);
	 fx = ((float)x)/((float)(1u<<maxDepth));
	 fy = ((float)y)/((float)(1u<<maxDepth));
	 fz = ((float)z)/((float)(1u<<maxDepth));
	 hx = ((float)(1u<<(maxDepth-d)))/((float)(1u<<maxDepth));
	 
	 for(int j = 0; j < 8; j++)
		{
		  float fxn,fyn,fzn;
		  fxn = fx + coord[j][0]*hx;
		  fyn = fy + coord[j][1]*hx;
		  fzn = fz + coord[j][2]*hx;
		  fprintf(outfile,"%f %f %f \n",fxn,fyn,fzn);
		}
	 fprintf(outfileV,"%u ",d);
  }

  fclose(outfileV);
  fclose(infile);
  fprintf(outfile,"\nCELLS %d %d\n",numNode,numNode*9);


  for(int i = 0; i < numNode; i++)
	 {
		fprintf(outfile,"8 ");
		
		for(int j = 0; j < 8; j++)
		  {
			 int index = (8*i)+j;
			 fprintf(outfile,"%d ", index);
		  }
		fprintf(outfile,"\n");
	 }

  fprintf(outfile,"\nCELL_TYPES %d\n",numNode);

  for(int i = 0; i < numNode; i++)
	 {
		fprintf(outfile,"11 \n");
	 }


  fprintf(outfile,"\nCELL_DATA %d\n",numNode);
  fprintf(outfile,"SCALARS scalars unsigned_int\n");
  fprintf(outfile,"LOOKUP_TABLE default\n");

  outfileV = fopen("Vfile","r");
  for (unsigned int i =0; i< numNode; i++) {
	 unsigned int v;
	 fscanf(outfileV,"%u",&v);
	 fprintf(outfile,"%u \n",v);
  }

  fclose(outfileV);
  fclose(outfile);
  return 0;
}
