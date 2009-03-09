
#include "mpi.h"
#include "petsc.h"
#include "sys.h"
#include "TreeNode.h"
#include "parUtils.h"
#include "omg.h"
#include "oda.h"
#include "omgJac.h"
#include "handleStencils.h"
#include "colors.h"
#include "externVars.h"
#include "dendro.h"

static char help[] = "FBM";

//Don't want time to be synchronized. Need to check load imbalance.
#ifdef MPI_WTIME_IS_GLOBAL
#undef MPI_WTIME_IS_GLOBAL
#endif

#ifdef PETSC_USE_LOG
//user-defined variables
int Jac1DiagEvent;
int Jac1MultEvent;
int Jac1FinestDiagEvent;
int Jac1FinestMultEvent;

int Jac2DiagEvent;
int Jac2MultEvent;
int Jac2FinestDiagEvent;
int Jac2FinestMultEvent;

int Jac3DiagEvent;
int Jac3MultEvent;
int Jac3FinestDiagEvent;
int Jac3FinestMultEvent;
#endif

double***** LaplacianType1Stencil; 
double**** LaplacianType2Stencil; 
double***** MassType1Stencil; 
double**** MassType2Stencil; 
double****** ShapeFnStencil;

int main(int argc, char ** argv ) {	
  int npes, rank;
  unsigned int maxNumPts= 1;
  unsigned int dim=3;
  unsigned int maxDepth=30;
  bool compressLut=false;
  double mgLoadFac = 2.0;
  bool incCorner = 1;  

  PetscInitialize(&argc,&argv,"optionsFBM",help);
  ot::RegisterEvents();

  ot::DAMG_Initialize(MPI_COMM_WORLD);

  MPI_Comm_size(MPI_COMM_WORLD,&npes);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  char pFile[256];
  sprintf(pFile, "fbmInp%d_%d.pts",rank, npes);

  std::vector<double> pts;
  ot::readPtsFromFile(pFile, pts);

  double gSize[3];
  gSize[0] = 1.0;
  gSize[1] = 1.0;
  gSize[2] = 1.0;

  //Construction
  std::vector<ot::TreeNode> linOct;
  ot::points2Octree(pts, gSize, linOct, dim, maxDepth, maxNumPts, MPI_COMM_WORLD);
  pts.clear();

  //Balancing...
  std::vector<ot::TreeNode> balOct;
  ot::balanceOctree(linOct, balOct, dim, maxDepth, incCorner, MPI_COMM_WORLD, NULL, NULL);
  linOct.clear();

  //Solve ...

  ot::DAMG        *damg;    
  int       nlevels = 1; //number of multigrid levels
  unsigned int       dof =1;// degrees of freedom per node  

  PetscInt nlevelsPetscInt = nlevels;
  PetscOptionsGetInt(0, "-nlevels", &nlevelsPetscInt, 0);
  nlevels = nlevelsPetscInt;

  if(!rank) {
    std::cout<<"nlevels initial: "<<nlevels<<std::endl;
  }

  // Note: The user context for all levels will be set separately later.
  ot::DAMGCreateAndSetDA(PETSC_COMM_WORLD, nlevels, NULL, &damg,
      balOct, dof, mgLoadFac, compressLut, incCorner);

  if(!rank) {
    std::cout<<"nlevels final: "<<nlevels<<std::endl;
  }

  if(!rank) {
    std::cout << "Created DA for all levels."<< std::endl;
  }

  ot::PrintDAMG(damg);

  createLmatType2(LaplacianType2Stencil);
  createMmatType2(MassType2Stencil);
  createShFnMat(ShapeFnStencil);

  //Global Function Handles for using KSP_Shell (will be used @ the coarsest grid if not all
  //processors are active on the coarsest grid)
  ot::getPrivateMatricesForKSP_Shell = getPrivateMatricesForKSP_Shell_Jac1;

  ot::DAMGSetKSP(damg, CreateJacobian1, ComputeJacobian1, ComputeFBM_RHS);

  ot::DAMGSolve(damg);

  destroyLmatType2(LaplacianType2Stencil);
  destroyMmatType2(MassType2Stencil);
  destroyShFnMat(ShapeFnStencil);

  DAMGDestroy(damg);

  balOct.clear();

  ot::DAMG_Finalize();

  PetscFinalize();
}//end function

