
#include <cstdio>
int main() {
  FILE* ifpL = fopen("LapType2Stencils.inp","r");
  FILE* ifpG = fopen("GDType2Stencils.inp","r");
  FILE* ofp = fopen("tmpStencils.txt","w");

  fprintf(ofp,"\nvoid InitRegularLaplacianStencil(double** stencil) {\n");

  for(int i = 0; i < 8; i++) {
    for(int j = 0; j < 8; j++) {
      double val;
      fscanf(ifpL,"%lf",&val);
      fprintf(ofp," stencil[%d][%d] = %lf; \n", i, j, val);
    }
  }

  fprintf(ofp,"\n}\n");
  fprintf(ofp,"\nvoid InitRegularGradDivStencil(double** stencil) {\n");

  for(int i = 0; i < 24; i++) {
    for(int j = 0; j < 24; j++) {
      double val;
      fscanf(ifpG,"%lf",&val);
      fprintf(ofp," stencil[%d][%d] = %lf; \n", i, j, val);
    }
  }

  fprintf(ofp,"\n}\n");

  fclose(ofp);
  fclose(ifpL);
  fclose(ifpG);

}

