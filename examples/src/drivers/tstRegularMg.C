
#include "mpi.h"
#include "petsc.h"
#include "sys.h"
#include "TreeNode.h"
#include "omg.h"
#include "oda.h"
#include "omgJac.h"
#include <cstdlib>
#include "handleStencils.h"
#include "colors.h"
#include "externVars.h"

static char help[] = "Regular Grid MG";

//Don't want time to be synchronized. Need to check load imbalance.
#ifdef MPI_WTIME_IS_GLOBAL
#undef MPI_WTIME_IS_GLOBAL
#endif

//user-defined variables

#ifdef PETSC_USE_LOG
int Jac1MultEvent;
int Jac1DiagEvent;
int Jac1FinestMultEvent;
int Jac1FinestDiagEvent;
int matPropEvent;

int Jac2MultEvent;
int Jac2DiagEvent;
int Jac2FinestMultEvent;
int Jac2FinestDiagEvent;

int Jac3MultEvent;
int Jac3DiagEvent;
int Jac3FinestMultEvent;
int Jac3FinestDiagEvent;
#endif

double***** LaplacianType1Stencil; 
double**** LaplacianType2Stencil; 
double***** MassType1Stencil; 
double**** MassType2Stencil; 
double****** ShapeFnStencil;

int main(int argc, char ** argv ) {	
  int size, rank;
  bool incCorner = 1;  
  unsigned int dim=3;
  unsigned int maxDepth=29;
  bool compressLut=true;
  std::vector<ot::TreeNode> balOct;
  double mgLoadFac = 2.0;
  unsigned int regLev = 2;

  PetscInitialize(&argc, &argv, "options", help);
  ot::RegisterEvents();

  ot::DAMG_Initialize(MPI_COMM_WORLD);

#ifdef PETSC_USE_LOG
  PetscClassId classid;
  PetscClassIdRegister("Dendro",&classid);
  
  PetscLogEventRegister("matProp",classid, &matPropEvent);
  PetscLogEventRegister("ODAmatDiag",classid, &Jac1DiagEvent);
  PetscLogEventRegister("ODAmatMult",classid, &Jac1MultEvent);
  PetscLogEventRegister("ODAmatDiagFinest",classid, &Jac1FinestDiagEvent);
  PetscLogEventRegister("ODAmatMultFinest",classid, &Jac1FinestMultEvent);

  PetscLogEventRegister("OMGmatDiag-2",classid, &Jac2DiagEvent);
  PetscLogEventRegister("OMGmatMult-2",classid, &Jac2MultEvent);
  PetscLogEventRegister("OMGmatDiagFinest-2",classid, &Jac2FinestDiagEvent);
  PetscLogEventRegister("OMGmatMultFinest-2",classid, &Jac2FinestMultEvent);

  PetscLogEventRegister("OMGmatDiag-3",classid, &Jac3DiagEvent);
  PetscLogEventRegister("OMGmatMult-3",classid, &Jac3MultEvent);
  PetscLogEventRegister("OMGmatDiagFinest-3",classid, &Jac3FinestDiagEvent);
  PetscLogEventRegister("OMGmatMultFinest-3",classid, &Jac3FinestMultEvent);

  int stages[1];
  PetscLogStageRegister("Solve",&stages[0]);  
#endif

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(argc > 1) {
    regLev = atoi(argv[1]);
  }
  if(argc > 2) {
    maxDepth = atoi(argv[2]);
  }
  if(argc > 3) {
    dim = atoi(argv[3]);
  }
  if(argc > 4) { incCorner = (bool)(atoi(argv[4]));}  
  if(argc > 5) { compressLut = (bool)(atoi(argv[5]));}
  if(argc > 6) { mgLoadFac = atof(argv[6]); }

#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[0]);
#endif
  MPI_Barrier(MPI_COMM_WORLD);	

  ot::DAMG        *damg;    
  int       nlevels = 1; //number of multigrid levels
  unsigned int       dof = 1; // degrees of freedom per node  

  createRegularOctree(balOct, regLev, dim, maxDepth, MPI_COMM_WORLD);

  PetscInt nlevelsPetscInt = nlevels; //To keep the compilers happy when using 64-bit indices 
  PetscOptionsGetInt(0, "-nlevels", &nlevelsPetscInt, 0);
  nlevels = nlevelsPetscInt;
  
  // Note: The user context for all levels will be set separately later.
  MPI_Barrier(MPI_COMM_WORLD);	
  
  if(!rank) {
    std::cout<<" nlevels initial: "<<nlevels<<std::endl;
  }
  
  ot::DAMGCreateAndSetDA(PETSC_COMM_WORLD, nlevels, NULL, &damg, 
      balOct, dof, mgLoadFac, compressLut, incCorner);

  if(!rank) {
    std::cout<<" nlevels final: "<<nlevels<<std::endl;
  }
  
  MPI_Barrier(MPI_COMM_WORLD);	
  
  if(!rank) {
    std::cout << "Created DA for all levels."<< std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  ot::PrintDAMG(damg);
  MPI_Barrier(MPI_COMM_WORLD);

  for(int i=0;i<nlevels;i++) {
    bool isRegOct = isRegularGrid(damg[i]->da);
    if(!rank) {
    std::cout<<"Level "<<i<<" is regular? "<<isRegOct<<std::endl;
    }
  }//end for i


  SetUserContexts(damg);

  if(!rank) {
    std::cout << "Set User Contexts all levels."<< std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  PetscInt       jacType = 1;
  PetscOptionsGetInt(0,"-jacType",&jacType,0);

  PetscInt       rhsType = 1;
  PetscOptionsGetInt(0,"-rhsType",&rhsType,0);

  createLmatType2(LaplacianType2Stencil);
  createMmatType2(MassType2Stencil);
  if(jacType == 3) {
    createLmatType1(LaplacianType1Stencil);
    createMmatType1(MassType1Stencil);
  }
 createShFnMat(ShapeFnStencil);

  if(!rank) {
    std::cout << "Created Stencils."<< std::endl;
  }

  //Function handles
  PetscErrorCode (*ComputeRHSHandle)(ot::DAMG damg,Vec rhs) = NULL;
  PetscErrorCode (*CreateJacobianHandle)(ot::DAMG damg,Mat *B) = NULL;
  PetscErrorCode (*ComputeJacobianHandle)(ot::DAMG damg,Mat J, Mat B) = NULL;

  if(rhsType == 0) {
    ComputeRHSHandle = ComputeRHS0;
  } else if (rhsType == 1) {
    ComputeRHSHandle = ComputeRHS1;
  } else if (rhsType == 2) {
    ComputeRHSHandle = ComputeRHS2;
  } else if (rhsType == 3) {
    ComputeRHSHandle = ComputeRHS3;
  } else if (rhsType == 4) {
    ComputeRHSHandle = ComputeRHS4;
  } else if (rhsType == 5) {
    ComputeRHSHandle = ComputeRHS5;
  } else if (rhsType == 6) {
    ComputeRHSHandle = ComputeRHS6;
  } else if (rhsType == 7) {
    ComputeRHSHandle = ComputeRHS7;
  } else if (rhsType == 8) {
    ComputeRHSHandle = ComputeRHS8;
  } else {
    assert(false);
  }

  if(jacType == 1) {
    CreateJacobianHandle = CreateJacobian1;
    ComputeJacobianHandle = ComputeJacobian1;
  } else if (jacType == 2) {
    CreateJacobianHandle = CreateJacobian2;
    ComputeJacobianHandle = ComputeJacobian2;
  } else if (jacType == 3) {
    CreateJacobianHandle = CreateJacobian3;
    ComputeJacobianHandle = ComputeJacobian3;
    //Skip the finest and the coarsest levels. For the other levels, J and B
    //must be different
    for(int i = 1; i < (nlevels-1); i++) {
      ot::DAMGCreateJMatrix(damg[i], CreateJacobianHandle);
    }
  } else {
    assert(false);
  }

  //Global Function Handles for using KSP_Shell (will be used @ the coarsest grid if not all
  //processors are active on the coarsest grid)
  if (jacType == 1) {
    ot::getPrivateMatricesForKSP_Shell = getPrivateMatricesForKSP_Shell_Jac1;
  } else if (jacType == 2) {
    ot::getPrivateMatricesForKSP_Shell = getPrivateMatricesForKSP_Shell_Jac2;
  } else if (jacType == 3) {
    ot::getPrivateMatricesForKSP_Shell = getPrivateMatricesForKSP_Shell_Jac3;
  } else {
    assert(false);
  }

  ot::DAMGSetKSP(damg, CreateJacobianHandle, ComputeJacobianHandle, ComputeRHSHandle);

  if(!rank) {
    std::cout<<"Solving u-Lu=f"<<std::endl;
  }

  iC(ot::DAMGSolve(damg));

  destroyLmatType2(LaplacianType2Stencil);
  destroyMmatType2(MassType2Stencil);
  if(jacType == 3) {
    destroyLmatType1(LaplacianType1Stencil);
    destroyMmatType1(MassType1Stencil);
  }
  destroyShFnMat(ShapeFnStencil);

  MPI_Barrier(MPI_COMM_WORLD);

  DestroyUserContexts(damg);

  if (!rank) {
    std::cout << GRN << "Destroyed User Contexts." << NRM << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  iC(DAMGDestroy(damg));

  if (!rank) {
    std::cout << GRN << "Destroyed DAMG" << NRM << std::endl;
  }

#ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif
  balOct.clear();

  if (!rank) {
    std::cout << GRN << "Finalizing ..." << NRM << std::endl;
  }
  ot::DAMG_Finalize();
  PetscFinalize();
}//end function

