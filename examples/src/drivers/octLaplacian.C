
/**
  @file octLaplacian.C
  @author Hari Sundar, hsundar@gmail.com
  */

#include "mpi.h"
#include "petsc.h"
#include "sys.h"
#include "petscmat.h"
#include "petscvec.h"
#include "petscksp.h"
#include "oda.h"
#include "massMatrix.h"
#include "stiffnessMatrix.h"
#include "externVars.h"

static char help[] = "Driver for a variable coefficient Laplacian problem on octrees";

#ifdef MPI_WTIME_IS_GLOBAL
#undef MPI_WTIME_IS_GLOBAL
#endif

PetscErrorCode LaplacianMatMult(Mat J , Vec in, Vec out);
PetscErrorCode InitializeData(ot::DA* da, Vec nu, Vec inc);

typedef struct {
  massMatrix* Mass;
  stiffnessMatrix* Stiffness;
} AppCtx;

int main(int argc, char **argv)
{
  unsigned int dof = 1;
  double nuval = 1.0;

  Vec nu;   /* Coefficient vector for the laplacian */
  Vec rhs, inc, sol;

  int size, rank;

  // Parameters for the balancing algorithm.
  // Refer to manual for details ... 
  bool incCorner = 1; // balance across corners = true  
  unsigned int maxNumPts= 1; // maximum number of points per octant
  unsigned int dim=3; // spatial dimensions 
  unsigned int maxDepth=30; // maximum depth of the octree, has to be <= 30

  char filePrefix[PETSC_MAX_PATH_LEN];
  char filename[PETSC_MAX_PATH_LEN];

  // The global domain size
  double gSize[3];
  gSize[0] = 1.; gSize[1] = 1.; gSize[2] = 1.;

  std::vector<ot::TreeNode> linOct, balOct;
  std::vector<double> pts;

  PetscInitialize(&argc, &argv, "options", NULL);
  ot::RegisterEvents();
  ot::DA_Initialize(MPI_COMM_WORLD);

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(argc < 3) {
    std::cerr << "Usage: " << argv[0] << "-pfx file_prefix" << std::endl;
    return -1;
  }

  CHKERRQ ( PetscOptionsGetString(PETSC_NULL, "-pfx", filePrefix, PETSC_MAX_PATH_LEN-1, PETSC_NULL));

  /*********************************************************************** */
  // READ: Read pts from file ...
  /*********************************************************************** */
  sprintf(filename, "%s.%d.pts", filePrefix, rank); 
  ot::readPtsFromFile(filename, pts);

  MPI_Barrier(MPI_COMM_WORLD);	

  /*********************************************************************** */
  // CONSTRUCT: Construct the linear octree from the points ...
  /*********************************************************************** */
  ot::points2Octree(pts, gSize, linOct, dim, maxDepth, maxNumPts, MPI_COMM_WORLD);

  // The points are not needed anymore, and can be cleared to free memory.
  pts.clear();

  /*********************************************************************** */
  // BALANCE: Balance the linear octree to enforce the 2:1 balance conditions.
  /*********************************************************************** */
  ot::balanceOctree (linOct, balOct, dim, maxDepth, incCorner, MPI_COMM_WORLD, NULL, NULL);

  // The linear octree (unbalanced) can be cleared to free memory.
  linOct.clear();

  // If desired the octree can be written to a file using the supplied routine ...
  //  ot::writeNodesToFile("filename.rank", balOct);

  /*********************************************************************** */
  // MESH : Construct the octree-based Distruted Array.
  /*********************************************************************** */
  assert(!(balOct.empty()));
  ot::DA da(balOct,MPI_COMM_WORLD, MPI_COMM_WORLD);
  balOct.clear();

  MPI_Barrier(MPI_COMM_WORLD);

  if (!rank)
    std::cout <<"Finshed Meshing" << std::endl;

  /*********************************************************************** */
  // Laplacian
  /*********************************************************************** */

  PetscOptionsGetScalar(0,"-nu",&nuval,0);

  // Use the ot::DA to create PETSc Vec objects
  da.createVector(inc, false, false, dof);
  da.createVector(nu, false, false, dof);
  da.createVector(rhs, false, false, dof);
  da.createVector(sol, false, false, dof);

  // Simple PETSc based initialization
  CHKERRQ( VecSet(nu, nuval ) );
  CHKERRQ( VecSet(inc, 1.0) );   

  // We can use elemental and nodal loops using the DA to initialize elemental 
  // and nodal properties.
  InitializeData(&da, nu, inc);

  /* ********************************************************************** */
  // create the Matrix objects
  massMatrix *Mass = new massMatrix(ot::fem::feMat::OCT); // Mass Matrix
  stiffnessMatrix *Stiffness = new stiffnessMatrix(ot::fem::feMat::OCT); // Stiffness matrix

  // set Matrix parameters
  Mass->setProblemDimensions(1.0,1.0,1.0);
  Mass->setDA(&da);
  Mass->setDof(dof);

  Stiffness->setProblemDimensions(1.0,1.0,1.0);
  Stiffness->setDA(&da);
  Stiffness->setNuVec(nu);
  Stiffness->setDof(dof);

  // The RHS
  Mass->MatVec(inc, rhs);

  // Create the user context
  AppCtx user;
  user.Stiffness = Stiffness;
  user.Mass = Mass;

  unsigned int Nx = da.getNodeSize();
  
  Mat J;
  MatCreateShell(PETSC_COMM_WORLD,dof*Nx, dof*Nx, PETSC_DETERMINE,PETSC_DETERMINE,&user,&J);
  MatShellSetOperation(J,MATOP_MULT,(void(*)(void))LaplacianMatMult);

  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD,&ksp);
  KSPSetOperators(ksp,J,J,SAME_NONZERO_PATTERN);
  KSPSetType(ksp,KSPCG);
  KSPSetFromOptions(ksp);

  // Solve 
  KSPSolve(ksp, rhs, sol);

  ot::DA_Finalize();
  PetscFinalize();

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "LaplacianMatMult"
PetscErrorCode LaplacianMatMult(Mat J , Vec in, Vec out)
{
  PetscFunctionBegin;
  AppCtx *data;

  MatShellGetContext( J, (void **)&data);  

  VecZeroEntries(out);
  data->Mass->MatVec(in,out);
  data->Stiffness->MatVec(in,out,-1.0);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "InitializeData"
PetscErrorCode InitializeData(ot::DA* da, Vec nu, Vec inc) {
  PetscScalar *inarray;
  PetscScalar *nuarray;

  unsigned int dof = 1;

  // Get pointers to the data buffers 
  da->vecGetBuffer(inc, inarray,false,false,false,dof);
  da->vecGetBuffer(nu, nuarray,false,false,false,dof);
  
  unsigned int maxD = da->getMaxDepth();
  unsigned int balOctmaxD = maxD - 1;
  unsigned int Nx = da->getNodeSize();

  double facsol = 1.0;

  for( da->init<ot::DA_FLAGS::ALL>(), da->init<ot::DA_FLAGS::WRITABLE>(); da->curr() < da->end<ot::DA_FLAGS::ALL>(); da->next<ot::DA_FLAGS::ALL>()) {
    
    Point pt;
    pt = da->getCurrentOffset();
    unsigned levelhere = da->getLevel(da->curr()) - 1;
    double hxOct = (double)((double)(1u << (balOctmaxD - levelhere))/(double)(1u << balOctmaxD));
    
    double x = (double)(pt.xint())/((double)(1u << (maxD-1)));
    double y = (double)(pt.yint())/((double)(1u << (maxD-1)));
    double z = (double)(pt.zint())/((double)(1u << (maxD-1)));

    int xindx = (static_cast<int>(x+y+z))%2;
    unsigned int indices[8];
    da->getNodeIndices(indices); 
    double coord[8][3] = {
      {0.0,0.0,0.0},
      {1.0,0.0,0.0},
      {0.0,1.0,0.0},
      {1.0,1.0,0.0},
      {0.0,0.0,1.0},
      {1.0,0.0,1.0},
      {0.0,1.0,1.0},
      {1.0,1.0,1.0}
    };

    unsigned char hn = da->getHangingNodeIndex(da->curr());

    for(int i = 0; i < 8; i++) {
      if (!(hn & (1 << i))) {
        double xhere, yhere, zhere;
        xhere = x + coord[i][0]*hxOct ; 
        yhere = y + coord[i][1]*hxOct; 
        zhere = z + coord[i][2]*hxOct; 
        
        nuarray[dof*indices[i]] = cos(facsol*M_PI*xhere)*cos(facsol*M_PI*yhere)*cos(facsol*M_PI*zhere);
        inarray[dof*indices[i]] = (1.0 + 3.0*facsol*facsol*M_PI*M_PI)*cos(facsol*M_PI*xhere)*cos(facsol*M_PI*yhere)*cos(facsol*M_PI*zhere);
      }
    }
  }

  da->vecRestoreBuffer(inc,inarray,false,false,false,dof);
  da->vecRestoreBuffer(nu,nuarray,false,false,false,dof);

  PetscFunctionReturn(0);
}

