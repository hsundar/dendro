
#include "mpi.h"
#include "petsc.h"
#include "sys.h"
#include "TreeNode.h"
#include "parUtils.h"
#include "omg.h"
#include "oda.h"
#include "omgJac.h"
#include "handleStencils.h"
#include <cstdlib>
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

  PetscInt       numRefinements = 0;

  PetscOptionsGetInt(0,"-numRefinements",&numRefinements,0);
  for(int i = 0; i < numRefinements; i++) {
    std::vector<ot::TreeNode> tmpOct = balOct;
    balOct.clear();
    ot::refineOctree(tmpOct, balOct); 
  }

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

  MPI_Barrier(MPI_COMM_WORLD);
  double setupStartTime = MPI_Wtime();
  
  // Note: The user context for all levels will be set separately later.
  ot::DAMGCreateAndSetDA(PETSC_COMM_WORLD, nlevels, NULL, &damg,
      balOct, dof, mgLoadFac, compressLut, incCorner);

  MPI_Barrier(MPI_COMM_WORLD);
  double setupEndTime = MPI_Wtime();

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

  ot::DAMGCreateSuppressedDOFs(damg);

  SetDirichletJacContexts(damg);

  //Global Function Handles for using KSP_Shell (will be used @ the coarsest grid if not all
  //processors are active on the coarsest grid)
  ot::getPrivateMatricesForKSP_Shell = getPrivateMatricesForKSP_Shell_DirichletJac;

  ot::DAMGSetKSP(damg, CreateDirichletJacobian, ComputeDirichletJacobian, ComputeFBM_RHS);

  MPI_Barrier(MPI_COMM_WORLD);
  double solveStartTime = MPI_Wtime();
  
  ot::DAMGSolve(damg);
  
  MPI_Barrier(MPI_COMM_WORLD);
  double solveEndTime = MPI_Wtime();

  Vec solTrue;
  VecDuplicate(DAMGGetx(damg), &solTrue);
  SetSolutionFBM(damg[nlevels - 1], solTrue);
  EnforceZeroFBM(damg[nlevels - 1], DAMGGetx(damg));

  VecAXPY(solTrue, -1.0, DAMGGetx(damg));

  PetscReal maxNormErr;
  VecNorm(solTrue, NORM_INFINITY, &maxNormErr);

  double l2err = ComputeFBMerror(damg[nlevels - 1], DAMGGetx(damg));

  VecDestroy(solTrue);

  if(!rank) {
    std::cout<<" Total Setup Time: "<<(setupEndTime - setupStartTime)<<std::endl;
    std::cout<<" Total Solve Time: "<<(solveEndTime - solveStartTime)<<std::endl;
    std::cout<<" maxNormErr (Pointwise): "<<maxNormErr<<std::endl;
    std::cout<<" l2err: "<<l2err<<std::endl;
  }

  destroyLmatType2(LaplacianType2Stencil);
  destroyMmatType2(MassType2Stencil);
  destroyShFnMat(ShapeFnStencil);

  DestroyDirichletJacContexts(damg);

  DAMGDestroy(damg);

  balOct.clear();

  ot::DAMG_Finalize();

  PetscFinalize();
}//end function

