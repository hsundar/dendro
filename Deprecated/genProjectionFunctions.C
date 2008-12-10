
#include <cstdio>

void getRemainingNodes(int v1, int v2, int& v3, int& v4);

int main() {

  int faceHn[][3] = {
    {3, 5, 6},
    {2, 4, 7},
    {1, 4, 7},
    {0, 5, 6},
    {1, 2, 7},
    {0, 3, 6},
    {0, 3, 5},
    {1, 2, 4}
  };

  int edgeHn[][3] = {
    {1, 2, 4},
    {0, 3, 5},
    {0, 3, 6},
    {1, 2, 7},
    {0, 5, 6},
    {1, 4, 7},
    {2, 4, 7},
    {3, 5, 6},
  };

  FILE* fptr = fopen("tmpFile.txt","w");

  //Project In
  fprintf(fptr,"\n inline void projectIn(unsigned char cNum, unsigned char hnMask,\n");
  fprintf(fptr,"unsigned int* indices, unsigned int dof, double* vals, double* projections) {\n");
  fprintf(fptr,"//Note: index cNum is always used for projections. This is the vertex that\n");
  fprintf(fptr,"//is shared with the parent\n"); 
  fprintf(fptr,"//The other index used in the projection is the vertex that is hanging\n");
  fprintf(fptr,"//For face hanging the remaining vertices are simply the remaining vertices on the\n");
  fprintf(fptr,"//plane formed by the above two vertices\n");
  fprintf(fptr,"switch(cNum) {\n");
  for(int cNum = 0 ; cNum < 8; cNum++) {
    fprintf(fptr," case %d: {\n",cNum);
    fprintf(fptr," for(unsigned int d = 0; d < dof; d++) {\n");
    fprintf(fptr,"  projections[(dof*%d) + d] = vals[(dof*indices[%d]) + d];\n",
        cNum, cNum);
    fprintf(fptr,"  projections[(dof*%d) + d] = vals[(dof*indices[%d]) + d];\n",
        (7 - cNum), (7 - cNum));
    fprintf(fptr," } // end for d\n");
    fprintf(fptr,"  //Edge Hanging: %d, %d, %d\n", edgeHn[cNum][0],
        edgeHn[cNum][1], edgeHn[cNum][2]); 
    for(int i = 0; i < 3; i++) {
      fprintf(fptr,"if(hnMask & (1 << %d)) {\n", edgeHn[cNum][i]);
      fprintf(fptr,"for(unsigned int d = 0; d < dof; d++) {\n");
      fprintf(fptr,"projections[(dof*%d) + d] = ", edgeHn[cNum][i]);
      fprintf(fptr,"0.5*(vals[(dof*indices[%d]) + d] + vals[(dof*indices[%d]) + d]);\n",
          edgeHn[cNum][i], cNum);
      fprintf(fptr," } // end for d\n");
      fprintf(fptr," } else {\n");
      fprintf(fptr," for(unsigned int d = 0; d < dof; d++) {\n");
      fprintf(fptr," projections[(dof*%d) + d] = vals[(dof*indices[%d]) + d];\n",
          edgeHn[cNum][i], edgeHn[cNum][i]);
      fprintf(fptr," } // end for d\n");
      fprintf(fptr,"} //end if\n");
    }//end for i
    fprintf(fptr,"  //Face Hanging: %d, %d, %d\n", faceHn[cNum][0],
        faceHn[cNum][1], faceHn[cNum][2]); 
    for(int i = 0; i < 3; i++) {
      int v3, v4;
      getRemainingNodes(cNum, faceHn[cNum][i], v3, v4);
      fprintf(fptr,"if(hnMask & (1 << %d)) {\n", faceHn[cNum][i]);
      fprintf(fptr,"for(unsigned int d = 0; d < dof; d++) {\n");
      fprintf(fptr,"projections[(dof*%d) + d] = ", faceHn[cNum][i]);
      fprintf(fptr,"0.25*(vals[(dof*indices[%d]) + d] + \n", faceHn[cNum][i]);
      fprintf(fptr,"vals[(dof*indices[%d]) + d] + ", v3);
      fprintf(fptr,"vals[(dof*indices[%d]) + d] + \n", v4);
      fprintf(fptr,"vals[(dof*indices[%d]) + d]);\n", cNum);
      fprintf(fptr," } // end for d\n");
      fprintf(fptr," } else {\n");
      fprintf(fptr," for(unsigned int d = 0; d < dof; d++) {\n");
      fprintf(fptr," projections[(dof*%d) + d] = vals[(dof*indices[%d]) + d];\n",
          faceHn[cNum][i], faceHn[cNum][i]);
      fprintf(fptr," } // end for d\n");
      fprintf(fptr,"} //end if\n");
    }//end for i
    fprintf(fptr,"break;\n");
    fprintf(fptr,"}\n\n");
  }//end for cNum
  fprintf(fptr,"default: {\n");
  fprintf(fptr,"assert(false);\n");
  fprintf(fptr,"}\n");
  fprintf(fptr,"} //end switch\n");
  fprintf(fptr,"} //end function\n");

  //Project Out 
  fprintf(fptr,"\n inline void projectOut(unsigned char cNum, unsigned char hnMask,\n");
  fprintf(fptr," unsigned int dof, double* vals, double* projections) {\n");
  fprintf(fptr,"//Note: Refer to the comments in projectIn. This is simply the reverse operation. \n");
  fprintf(fptr,"switch(cNum) {\n");
  for(int cNum = 0 ; cNum < 8; cNum++) {
    fprintf(fptr," case %d: {\n",cNum);
    fprintf(fptr," for(unsigned int d = 0; d < dof; d++) {\n");
    fprintf(fptr," vals[(dof*%d) + d] = projections[(dof*%d) + d];\n",
        cNum, cNum);
    fprintf(fptr," vals[(dof*%d) + d] = projections[(dof*%d) + d];\n",
        (7 - cNum), (7 - cNum));
    fprintf(fptr," } // end for d\n");
    fprintf(fptr,"  //Edge Hanging: %d, %d, %d\n", edgeHn[cNum][0],
        edgeHn[cNum][1], edgeHn[cNum][2]); 
    for(int i = 0; i < 3; i++) {
      fprintf(fptr,"if(hnMask & (1 << %d)) {\n", edgeHn[cNum][i]);
      fprintf(fptr,"for(unsigned int d = 0; d < dof; d++) {\n");
      fprintf(fptr," vals[(dof*%d) + d] = ", edgeHn[cNum][i]);
      fprintf(fptr,"(2.0*projections[(dof*%d) + d]) - vals[(dof*%d) + d];\n",
          edgeHn[cNum][i], cNum);
      fprintf(fptr," } // end for d\n");
      fprintf(fptr," } else {\n");
      fprintf(fptr," for(unsigned int d = 0; d < dof; d++) {\n");
      fprintf(fptr," vals[(dof*%d) + d] = projections[(dof*%d) + d];\n",
          edgeHn[cNum][i], edgeHn[cNum][i]);
      fprintf(fptr," } // end for d\n");
      fprintf(fptr,"} //end if\n");
    } //end for i
    fprintf(fptr,"  //Face Hanging: %d, %d, %d\n", faceHn[cNum][0],
        faceHn[cNum][1], faceHn[cNum][2]); 
    for(int i = 0; i < 3; i++) {
      int v3, v4;
      getRemainingNodes(cNum, faceHn[cNum][i], v3, v4);
      fprintf(fptr,"if(hnMask & (1 << %d)) {\n", faceHn[cNum][i]);
      fprintf(fptr,"for(unsigned int d = 0; d < dof; d++) {\n");
      fprintf(fptr,"vals[(dof*%d) + d] = ", faceHn[cNum][i]);
      fprintf(fptr,"(4.0*projections[(dof*%d) + d]) -  \n", faceHn[cNum][i]);
      fprintf(fptr,"(vals[(dof*%d) + d] + ", v3);
      fprintf(fptr,"vals[(dof*%d) + d] + \n", v4);
      fprintf(fptr,"vals[(dof*%d) + d]);\n", cNum);
      fprintf(fptr," } // end for d\n");
      fprintf(fptr," } else {\n");
      fprintf(fptr," for(unsigned int d = 0; d < dof; d++) {\n");
      fprintf(fptr," vals[(dof*%d) + d] = projections[(dof*%d) + d];\n",
          faceHn[cNum][i], faceHn[cNum][i]);
      fprintf(fptr," } // end for d\n");
      fprintf(fptr,"} //end if\n");
    } //end for i
    fprintf(fptr,"break;\n");
    fprintf(fptr,"}\n\n");
  }//end for cNum
  fprintf(fptr,"default: {\n");
  fprintf(fptr,"assert(false);\n");
  fprintf(fptr,"}\n");
  fprintf(fptr,"} //end switch\n");
  fprintf(fptr,"} //end function\n");

  fclose(fptr);
}

void getRemainingNodes(int v1, int v2, int& v3, int& v4) {

  int surfaces[][4] = {
    {0, 1, 2, 3},
    {4, 5, 6, 7},
    {1, 3, 5, 7},
    {0, 2, 4, 6},
    {0, 1, 4, 5},
    {2, 3, 6, 7}
  };

  int row = 0;
  for(int i = 0; i < 6; i++) {
    bool foundV1 = false;
    bool foundV2 = false;
    for(int j = 0; j < 4; j++) {
      if(surfaces[i][j] == v1) {
        foundV1 = true;
      }
      if(surfaces[i][j] == v2) {
        foundV2 = true;
      }
    } //end for j
    if(foundV1 && foundV2) {
      row = i;
      break;
    }
  } //end for i

  bool foundV3 = false;
  for(int j = 0; j < 4; j++) {
    if( (surfaces[row][j] != v1) && (surfaces[row][j] != v2) ) {
      if(foundV3) {
        v4 = surfaces[row][j];
      } else {
        v3 = surfaces[row][j];
        foundV3 = true;
      }
    }
  } //end for j

} //end function


