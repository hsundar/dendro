
/**
  @file omgRhs.C
  @author Rahul Sampath, rahul.sampath@gmail.com
  */

#include "petsc.h"
#include "parUtils.h"
#include "seqUtils.h"
#include "omg.h"
#include "oda.h"
#include "odaJac.h"
#include "omgJac.h"
#include "nodeAndValues.h"

extern double****** ShapeFnStencil;
extern double**** ShapeFnCoeffs;

#define square(x) ((x)*(x))

#define EVAL_FN(x,y,z,result) {\
  result = 0.0;\
  if( ((x) >= 0.0) && ((y) >= 0.0) && ((z) >= 0.0) && \
      ((x) < 0.2) && ((y) < 0.2) && ((z) < 0.2) ) {\
    double dx = ((x) - 0.1);\
    double dy = ((y) - 0.1);\
    double dz = ((z) - 0.1);\
    result += exp(-((square(dx) + square(dy) + square(dz))/0.005));\
  }\
  if( ((x) >= 0.2) && ((y) >= 0.2) && ((z) >= 0.2) && \
      ((x) <= 0.4) && ((y) <= 0.4) && ((z) <= 0.4) ) {\
    double dx = ((x) - 0.3);\
    double dy = ((y) - 0.3);\
    double dz = ((z) - 0.3);\
    result += exp(-((square(dx) + square(dy) + square(dz))/0.005));\
  }\
  if( ((x) >= 0.5) && ((y) >= 0.5) && ((z) >= 0.5) && \
      ((x) <= 0.7) && ((y) <= 0.7) && ((z) <= 0.7) ) {\
    double dx = ((x) - 0.6);\
    double dy = ((y) - 0.6);\
    double dz = ((z) - 0.6);\
    result += exp(-((square(dx) + square(dy) + square(dz))/0.005));\
  }\
}

//Rhs = sum of 3 gaussian functions centered at different points.
//Rhs is assembled using 3-pt gauss-quadrature integration per dimension. This is exact to
//polynomials of degree 5 or less
//w1 =8/9, x1 =0, w2=w3 =5/9, x2 = +sqrt(3/5), x3 = -sqrt(3/5)
//int_a_to_b(f) = ((b-a)/2)*sum(wi*f((((b-a)/2)*xi) + ((b+a)/2)))
PetscErrorCode ComputeRHS2(ot::DAMG damg,Vec in) {
  PetscFunctionBegin;	 	 
  ot::DA* da = damg->da;
  PetscScalar *inarray;

  //In PETSc's debug mode they use an Allreduce where all processors (active
  //and inactive) participate. Hence, VecZeroEntries must be called by active
  //and inactive processors.
  VecZeroEntries(in);
  da->vecGetBuffer(in,inarray,false,false,false,1);

  unsigned int maxD;
  unsigned int balOctmaxD;

  double wts[3] = { (8.0/9.0), (5.0/9.0), (5.0/9.0) };
  double gPts[3] = { 0.0, sqrt((3.0/5.0)), -sqrt((3.0/5.0)) };

  if(da->iAmActive()) {
    maxD = da->getMaxDepth();
    balOctmaxD = maxD - 1;
    for(da->init<ot::DA_FLAGS::ALL>();
        da->curr() < da->end<ot::DA_FLAGS::ALL>();
        da->next<ot::DA_FLAGS::ALL>())  
    {
      Point pt;
      pt = da->getCurrentOffset();
      unsigned int idx = da->curr();
      unsigned levelhere = (da->getLevel(idx) - 1);
      double hxOct = (double)((double)(1u << (balOctmaxD - levelhere))/(double)(1u << balOctmaxD));
      double x = (double)(pt.xint())/((double)(1u << (maxD-1)));
      double y = (double)(pt.yint())/((double)(1u << (maxD-1)));
      double z = (double)(pt.zint())/((double)(1u << (maxD-1)));
      double fac = ((hxOct*hxOct*hxOct)/8.0);
      unsigned int indices[8];
      da->getNodeIndices(indices); 
      unsigned char childNum = da->getChildNumber();
      unsigned char hnMask = da->getHangingNodeIndex(idx);
      unsigned char elemType = 0;
      GET_ETYPE_BLOCK(elemType,hnMask,childNum)
        for(unsigned int j = 0; j < 8; j++) {
          double integral = 0.0;
          //Quadrature Rule
          for(int m = 0; m < 3; m++) {
            for(int n = 0; n < 3; n++) {
              for(int p = 0; p < 3; p++) {
                double fnVal = 0.0;
                double xPt = ( (hxOct*(1.0 +gPts[m])*0.5) + x );
                double yPt = ( (hxOct*(1.0 + gPts[n])*0.5) + y );
                double zPt = ( (hxOct*(1.0 + gPts[p])*0.5) + z );
                EVAL_FN(xPt,yPt,zPt,fnVal) 
                  integral +=
                  (wts[m]*wts[n]*wts[p]*fnVal*
                   ShapeFnStencil[childNum][elemType][j][m][n][p]);
              }
            }
          }
          inarray[indices[j]] += (fac*integral);
        }//end for j
    }//end for i
  }//end if active

  da->vecRestoreBuffer(in,inarray,false,false,false,1);

  PetscFunctionReturn(0);
}

#undef square
#undef EVAL_FN 

PetscErrorCode ComputeSol(ot::DAMG damg,Vec expectedSol) {
  PetscFunctionBegin;	 	 
  ot::DA* da = damg->da;
  PetscScalar *solArray;
  VecZeroEntries(expectedSol);
  da->vecGetBuffer(expectedSol,solArray,false,false,false,1);

  unsigned int maxD;
  unsigned int balOctmaxD;

  if(da->iAmActive()) {
    maxD = da->getMaxDepth();
    balOctmaxD = maxD - 1;
    for(da->init<ot::DA_FLAGS::ALL>(); 
        da->curr() < da->end<ot::DA_FLAGS::ALL>();
        da->next<ot::DA_FLAGS::ALL>())  
    {
      Point pt;
      pt = da->getCurrentOffset();
      unsigned levelhere = da->getLevel(da->curr()) - 1;
      double hxOct = (double)((double)(1u << (balOctmaxD - levelhere))/(double)(1u << balOctmaxD));
      double x = (double)(pt.xint())/((double)(1u << (maxD-1)));
      double y = (double)(pt.yint())/((double)(1u << (maxD-1)));
      double z = (double)(pt.zint())/((double)(1u << (maxD-1)));

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

      for(int i = 0; i < 8; i++)
      {
        if (!(hn & (1 << i))){
          double xhere, yhere, zhere;
          xhere = x + coord[i][0]*hxOct ; yhere = y + coord[i][1]*hxOct; zhere = z + coord[i][2]*hxOct; 
          double solSum = 0.0;
          for(int freqCnt = 1; freqCnt < 10; freqCnt++) {
            double facsol = freqCnt;
            solSum  += cos(facsol*M_PI*xhere)*cos(facsol*M_PI*yhere)*cos(facsol*M_PI*zhere);				               
          }
          solArray[indices[i]] = solSum;
        }
      }
    }
  }

  da->vecRestoreBuffer(expectedSol,solArray,false,false,false,1); 

  PetscFunctionReturn(0);
}

PetscErrorCode ComputeRHS0(ot::DAMG damg,Vec in) {
  PetscFunctionBegin;
  VecZeroEntries(in);
  PetscFunctionReturn(0);
}

PetscErrorCode ComputeRandomRHS(ot::DAMG damg,Vec in) {
  PetscFunctionBegin;
  PetscRandom rctx;  
  PetscRandomCreate(PETSC_COMM_WORLD,&rctx);
  PetscRandomSetType(rctx,PETSCRAND48);
  PetscInt randomSeed = 12345;
  PetscOptionsGetInt(0,"-randomSeed",&randomSeed,0);
  int rank;
  MPI_Comm_rank(damg->comm,&rank);
  if(!rank) {
    std::cout<<"Using Random Seed: "<<randomSeed<<std::endl;
  }
  PetscRandomSetSeed(rctx,randomSeed);
  PetscRandomSeed(rctx);
  PetscRandomSetFromOptions(rctx);
  VecSetRandom(in,rctx);
  PetscRandomDestroy(rctx);

  PetscReal norm2;
  PetscReal normInf;
  VecNorm(in, NORM_INFINITY, &normInf);
  VecNorm(in, NORM_2, &norm2);
  if(!rank) {
    std::cout<<"Random Rhs norm-2: "<<norm2<<" normInf: "<<normInf<<std::endl; 
  }

  PetscFunctionReturn(0);
}

//Consistent Random RHS
PetscErrorCode ComputeRHS3(ot::DAMG damg,Vec in) {
  PetscFunctionBegin;

  int rank;
  MPI_Comm_rank(damg->comm,&rank);
  Vec tmp;
  VecDuplicate(in,&tmp);
  ComputeRandomRHS(damg,tmp);
  Mat massMat;
  CreateAndComputeMassMatrix(damg->da,&massMat);
  MatMult(massMat,tmp,in);
  MatDestroy(massMat);
  VecDestroy(tmp);

  PetscReal norm2;
  PetscReal normInf;
  VecNorm(in, NORM_INFINITY, &normInf);
  VecNorm(in, NORM_2, &norm2);
  if(!rank) {
    std::cout<<"End of RHS-3: Rhs norm-2: "<<norm2<<" normInf: "<<normInf<<std::endl; 
  }
  PetscFunctionReturn(0);
}

//Random Solution Consistent RHS
PetscErrorCode ComputeRHS4(ot::DAMG damg, Vec in) {
  PetscFunctionBegin;
  int rank;
  MPI_Comm_rank(damg->comm, &rank);
  Vec tmp;
  VecDuplicate(in, &tmp);
  ComputeRandomRHS(damg, tmp);

  MatMult(damg->J, tmp, in);
  VecDestroy(tmp);

  PetscReal norm2;
  PetscReal normInf;
  VecNorm(in, NORM_INFINITY, &normInf);
  VecNorm(in, NORM_2, &norm2);

  if(!rank) {
    std::cout<<" End of RHS-4: Rhs norm-2: "<<norm2<<" normInf: "<<normInf<<std::endl; 
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ComputeRHS1(ot::DAMG damg,Vec in) {
  PetscFunctionBegin;	 	 
  ot::DA* da = damg->da;
  Vec tmp;
  VecDuplicate(in,&tmp);
  PetscScalar *inarray;
  VecZeroEntries(tmp);
  da->vecGetBuffer(tmp,inarray,false,false,false,1);

  int rank;
  MPI_Comm_rank(damg->comm,&rank);
  unsigned int maxD;
  unsigned int balOctmaxD;

  if(da->iAmActive()) {
    maxD = da->getMaxDepth();
    balOctmaxD = maxD - 1;
    for(da->init<ot::DA_FLAGS::ALL>(); da->curr() < da->end<ot::DA_FLAGS::ALL>(); da->next<ot::DA_FLAGS::ALL>())  
    {
      Point pt;
      pt = da->getCurrentOffset();
      unsigned levelhere = da->getLevel(da->curr()) - 1;
      double hxOct = (double)((double)(1u << (balOctmaxD - levelhere))/(double)(1u << balOctmaxD));
      double x = (double)(pt.xint())/((double)(1u << (maxD-1)));
      double y = (double)(pt.yint())/((double)(1u << (maxD-1)));
      double z = (double)(pt.zint())/((double)(1u << (maxD-1)));

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

      for(int i = 0; i < 8; i++)
      {
        if (!(hn & (1 << i))){
          double xhere, yhere, zhere;
          xhere = x + coord[i][0]*hxOct ; yhere = y + coord[i][1]*hxOct; zhere = z + coord[i][2]*hxOct; 
          double rhsSum = 0.0;
          for(int freqCnt = 1; freqCnt < 10; freqCnt++) {
            double facsol = freqCnt;
            rhsSum  += (1.0 + 3*facsol*facsol*M_PI*M_PI)*cos(facsol*M_PI*xhere)*
              cos(facsol*M_PI*yhere)*cos(facsol*M_PI*zhere);				               
          }
          inarray[indices[i]] = rhsSum;
        }
      }
    }
  }

  da->vecRestoreBuffer(tmp,inarray,false,false,false,1); 

  PetscReal norm2;
  PetscReal normInf;

  VecNorm(tmp, NORM_INFINITY, &normInf);
  VecNorm(tmp, NORM_2, &norm2);
  if(!rank) {
    std::cout<<"tmpRhs norm-2: "<<norm2<<" normInf: "<<normInf<<std::endl; 
  }

  Mat massMat;
  CreateAndComputeMassMatrix(da,&massMat);
  MatMult(massMat,tmp,in);
  MatDestroy(massMat);
  VecDestroy(tmp);

  VecNorm(in, NORM_INFINITY, &normInf);
  VecNorm(in, NORM_2, &norm2);
  if(!rank) {
    std::cout<<"End of RHS-1: Rhs norm-2: "<<norm2<<" normInf: "<<normInf<<std::endl; 
  }

  PetscFunctionReturn(0);
}



/*
   RHS corresponding to the homogeneous neumann, variable coefficient, scalar, poisson problem
   with solution = cos(2*pi*x)cos(2*pi*y)cos(2*pi*z) and
   epsilon = 1 + 10^6[cos^2(2*pi*x) + cos^2(2*pi*y) + cos^2(2*pi*z)] and
   alpha  = 1
   */
PetscErrorCode ComputeRHS5(ot::DAMG damg,Vec in) {
  PetscFunctionBegin;	 	 
  ot::DA* da = damg->da;
  Vec tmp;
  VecDuplicate(in, &tmp);
  PetscScalar *inarray;
  VecZeroEntries(tmp);
  da->vecGetBuffer(tmp, inarray, false, false, false, 1);

  PetscReal lapFac = 0.0;
  PetscTruth optFound;
  PetscOptionsGetReal("lap","-MatPropFac",&lapFac,&optFound);
  int rank;
  MPI_Comm_rank(damg->comm,&rank);
  unsigned int maxD;
  unsigned int balOctmaxD;

  if(da->iAmActive()) {
    maxD = da->getMaxDepth();
    balOctmaxD = maxD - 1;
    for(da->init<ot::DA_FLAGS::ALL>(); da->curr() < da->end<ot::DA_FLAGS::ALL>(); da->next<ot::DA_FLAGS::ALL>())  
    {
      Point pt;
      pt = da->getCurrentOffset();
      unsigned levelhere = da->getLevel(da->curr()) - 1;
      double hxOct = (double)((double)(1u << (balOctmaxD - levelhere))/(double)(1u << balOctmaxD));
      double x = (double)(pt.xint())/((double)(1u << (maxD-1)));
      double y = (double)(pt.yint())/((double)(1u << (maxD-1)));
      double z = (double)(pt.zint())/((double)(1u << (maxD-1)));

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

      for(int i = 0; i < 8; i++)
      {
        if (!(hn & (1 << i))){
          double xhere, yhere, zhere;
          xhere = x + (coord[i][0]*hxOct);
          yhere = y + (coord[i][1]*hxOct);
          zhere = z + (coord[i][2]*hxOct); 
          double solVal = cos(2.0*M_PI*xhere)*cos(2.0*M_PI*yhere)*cos(2.0*M_PI*zhere); 
          double epVal = (1.0 + (lapFac*(
                  (cos(2.0*M_PI*xhere)*cos(2.0*M_PI*xhere)) +
                  (cos(2.0*M_PI*yhere)*cos(2.0*M_PI*yhere)) +
                  (cos(2.0*M_PI*zhere)*cos(2.0*M_PI*zhere)) 
                  )));
          double rhsVal = (solVal + (12.0*(M_PI*M_PI)*solVal*epVal)
              - (8.0*(M_PI*M_PI)*lapFac*solVal*(
                  (sin(2.0*M_PI*xhere)*sin(2.0*M_PI*xhere)) +
                  (sin(2.0*M_PI*yhere)*sin(2.0*M_PI*yhere)) +
                  (sin(2.0*M_PI*zhere)*sin(2.0*M_PI*zhere)) 
                  )));
          inarray[indices[i]] = rhsVal;
        }
      }
    }
  }

  da->vecRestoreBuffer(tmp, inarray, false, false, false, 1); 

  PetscReal norm2;
  PetscReal normInf;

  VecNorm(tmp, NORM_INFINITY, &normInf);
  VecNorm(tmp, NORM_2, &norm2);
  if(!rank) {
    std::cout<<"tmpRhs-5 norm-2: "<<norm2<<" normInf: "<<normInf<<std::endl; 
  }

  Mat massMat;
  CreateAndComputeMassMatrix(da,&massMat);
  MatMult(massMat,tmp,in);
  MatDestroy(massMat);
  VecDestroy(tmp);

  VecNorm(in, NORM_INFINITY, &normInf);
  VecNorm(in, NORM_2, &norm2);
  if(!rank) {
    std::cout<<"End of RHS-5: Rhs norm-2: "<<norm2<<" normInf: "<<normInf<<std::endl; 
  }

  PetscFunctionReturn(0);
}

/**
  RHS assembled by multiplying the J matrix with the solution
  */
PetscErrorCode ComputeRHS6(ot::DAMG damg,Vec in) {
  PetscFunctionBegin;	 	 
  ot::DA* da = damg->da;
  Vec tmp;
  VecDuplicate(in, &tmp);
  SetSolution5(da,tmp);
  MatMult(damg->J, tmp, in);
  VecDestroy(tmp);
  PetscFunctionReturn(0);
}

/**
  RHS corresponding to the homogeneous neumann, variable coefficient, scalar, poisson problem
  with solution = cos(2*pi*x)cos(2*pi*y)cos(2*pi*z) and
  epsilon = 1 + 10^6[cos^2(2*pi*x) + cos^2(2*pi*y) + cos^2(2*pi*z)] and
  alpha  = 1
  This is the same as RHStypes 5 and 6. It is only assembled differently.
  Rhs is assembled using 3-pt gauss-quadrature integration per dimension. This is exact to
  polynomials of degree 5 or less
  w1 =8/9, x1 =0, w2=w3 =5/9, x2 = +sqrt(3/5), x3 = -sqrt(3/5)
  int_a_to_b(f) = ((b-a)/2)*sum(wi*f((((b-a)/2)*xi) + ((b+a)/2)))
  */
PetscErrorCode ComputeRHS7(ot::DAMG damg,Vec in) {
  PetscFunctionBegin;	 	 
  ot::DA* da = damg->da;
  PetscScalar *inarray;
  VecZeroEntries(in);
  da->vecGetBuffer(in,inarray,false,false,false,1);

  PetscReal lapFac = 0.0;
  PetscTruth optFound;
  PetscOptionsGetReal("lap","-MatPropFac",&lapFac,&optFound);

  unsigned int maxD;
  unsigned int balOctmaxD;

  double wts[3] = { (8.0/9.0), (5.0/9.0), (5.0/9.0) };
  double gPts[3] = { 0.0, sqrt((3.0/5.0)), -sqrt((3.0/5.0)) };

  if(da->iAmActive()) {
    maxD = da->getMaxDepth();
    balOctmaxD = maxD - 1;
    for(da->init<ot::DA_FLAGS::ALL>();
        da->curr() < da->end<ot::DA_FLAGS::ALL>();
        da->next<ot::DA_FLAGS::ALL>())  
    {
      Point pt;
      pt = da->getCurrentOffset();
      unsigned int idx = da->curr();
      unsigned levelhere = (da->getLevel(idx) - 1);
      double hxOct = (double)((double)(1u << (balOctmaxD - levelhere))/(double)(1u << balOctmaxD));
      double x = (double)(pt.xint())/((double)(1u << (maxD-1)));
      double y = (double)(pt.yint())/((double)(1u << (maxD-1)));
      double z = (double)(pt.zint())/((double)(1u << (maxD-1)));
      double fac = ((hxOct*hxOct*hxOct)/8.0);
      unsigned int indices[8];
      da->getNodeIndices(indices); 
      unsigned char childNum = da->getChildNumber();
      unsigned char hnMask = da->getHangingNodeIndex(idx);
      unsigned char elemType = 0;
      GET_ETYPE_BLOCK(elemType,hnMask,childNum)
        for(unsigned int j = 0; j < 8; j++) {
          double integral = 0.0;
          //Quadrature Rule
          for(int m = 0; m < 3; m++) {
            for(int n = 0; n < 3; n++) {
              for(int p = 0; p < 3; p++) {
                double xPt = ( (hxOct*(1.0 +gPts[m])*0.5) + x );
                double yPt = ( (hxOct*(1.0 + gPts[n])*0.5) + y );
                double zPt = ( (hxOct*(1.0 + gPts[p])*0.5) + z );
                double solVal = cos(2.0*M_PI*xPt)*cos(2.0*M_PI*yPt)*cos(2.0*M_PI*zPt); 
                double epVal = (1.0 + (lapFac*(
                        (cos(2.0*M_PI*xPt)*cos(2.0*M_PI*xPt)) +
                        (cos(2.0*M_PI*yPt)*cos(2.0*M_PI*yPt)) +
                        (cos(2.0*M_PI*zPt)*cos(2.0*M_PI*zPt)) 
                        )));
                double rhsVal = (solVal + (12.0*(M_PI*M_PI)*solVal*epVal)
                    - (8.0*(M_PI*M_PI)*lapFac*solVal*(
                        (sin(2.0*M_PI*xPt)*sin(2.0*M_PI*xPt)) +
                        (sin(2.0*M_PI*yPt)*sin(2.0*M_PI*yPt)) +
                        (sin(2.0*M_PI*zPt)*sin(2.0*M_PI*zPt)) 
                        )));
                integral += (wts[m]*wts[n]*wts[p]*rhsVal*
                    ShapeFnStencil[childNum][elemType][j][m][n][p]);
              }
            }
          }
          inarray[indices[j]] += (fac*integral);
        }//end for j
    }//end for i
  }//end if active

  da->vecRestoreBuffer(in,inarray,false,false,false,1);

  PetscFunctionReturn(0);
}

/**
  RHS corresponding to the homogeneous neumann, variable coefficient, scalar, poisson problem
  with solution = cos(2*pi*x)cos(2*pi*y)cos(2*pi*z) and
  epsilon = 1 + 10^6[cos^2(2*pi*x) + cos^2(2*pi*y) + cos^2(2*pi*z)] and
  alpha  = 1
  This is the same as RHStypes 5, 6 and 7. It is only assembled differently.
  Rhs is assembled using arbitrary order(>=2 and <=7) gauss-quadrature integration per dimension. This is exact to
  polynomials of degree 5 or less
  w1 =8/9, x1 =0, w2=w3 =5/9, x2 = +sqrt(3/5), x3 = -sqrt(3/5)
  int_a_to_b(f) = ((b-a)/2)*sum(wi*f((((b-a)/2)*xi) + ((b+a)/2)))
  */
PetscErrorCode ComputeRHS8(ot::DAMG damg,Vec in) {
  PetscFunctionBegin;	 	 
  ot::DA* da = damg->da;
  PetscScalar *inarray;
  VecZeroEntries(in);
  da->vecGetBuffer(in,inarray,false,false,false,1);

  PetscReal lapFac = 0.0;
  PetscTruth optFound;
  PetscOptionsGetReal("lap","-MatPropFac",&lapFac,&optFound);

  PetscInt numGaussPts = 0;
  PetscOptionsGetInt(0,"-numGaussPts",&numGaussPts,0);

  assert(numGaussPts <= 7);
  assert(numGaussPts >= 2);

  unsigned int maxD;
  unsigned int balOctmaxD;

  std::vector<std::vector<double> > wts(6);
  std::vector<std::vector<double> > gPts(6);

  for(int k = 0; k < 6; k++) {
    wts[k].resize(k+2);
    gPts[k].resize(k+2);
  }

  //2-pt rule
  wts[0][0] = 1.0; wts[0][1] = 1.0;
  gPts[0][0] = sqrt(1.0/3.0); gPts[0][1] = -sqrt(1.0/3.0); 

  //3-pt rule
  wts[1][0] = 0.88888889;  wts[1][1] = 0.555555556;  wts[1][2] = 0.555555556;
  gPts[1][0] = 0.0;  gPts[1][1] = 0.77459667;  gPts[1][2] = -0.77459667;

  //4-pt rule
  wts[2][0] = 0.65214515;  wts[2][1] = 0.65214515;
  wts[2][2] = 0.34785485; wts[2][3] = 0.34785485;  
  gPts[2][0] = 0.33998104;  gPts[2][1] = -0.33998104;
  gPts[2][2] = 0.86113631; gPts[2][3] = -0.86113631;

  //5-pt rule
  wts[3][0] = 0.568888889;  wts[3][1] = 0.47862867;  wts[3][2] =  0.47862867;
  wts[3][3] = 0.23692689; wts[3][4] = 0.23692689;
  gPts[3][0] = 0.0;  gPts[3][1] = 0.53846931; gPts[3][2] = -0.53846931;
  gPts[3][3] = 0.90617985; gPts[3][4] = -0.90617985;

  //6-pt rule
  wts[4][0] = 0.46791393;  wts[4][1] = 0.46791393;  wts[4][2] = 0.36076157;
  wts[4][3] = 0.36076157; wts[4][4] = 0.17132449; wts[4][5] = 0.17132449;
  gPts[4][0] = 0.23861918; gPts[4][1] = -0.23861918; gPts[4][2] = 0.66120939;
  gPts[4][3] = -0.66120939; gPts[4][4] = 0.93246951; gPts[4][5] = -0.93246951;

  //7-pt rule
  wts[5][0] = 0.41795918;  wts[5][1] = 0.38183005;  wts[5][2] = 0.38183005;
  wts[5][3] = 0.27970539;  wts[5][4] = 0.27970539; 
  wts[5][5] = 0.12948497; wts[5][6] = 0.12948497;
  gPts[5][0] = 0.0;  gPts[5][1] = 0.40584515;  gPts[5][2] = -0.40584515;
  gPts[5][3] = 0.74153119;  gPts[5][4] = -0.74153119;
  gPts[5][5] = 0.94910791; gPts[5][6] = -0.94910791;

  if(da->iAmActive()) {
    maxD = da->getMaxDepth();
    balOctmaxD = maxD - 1;
    for(da->init<ot::DA_FLAGS::ALL>();
        da->curr() < da->end<ot::DA_FLAGS::ALL>();
        da->next<ot::DA_FLAGS::ALL>())  
    {
      Point pt;
      pt = da->getCurrentOffset();
      unsigned int idx = da->curr();
      unsigned levelhere = (da->getLevel(idx) - 1);
      double hxOct = (double)((double)(1u << (balOctmaxD - levelhere))/(double)(1u << balOctmaxD));
      double x = (double)(pt.xint())/((double)(1u << (maxD-1)));
      double y = (double)(pt.yint())/((double)(1u << (maxD-1)));
      double z = (double)(pt.zint())/((double)(1u << (maxD-1)));
      double fac = ((hxOct*hxOct*hxOct)/8.0);
      unsigned int indices[8];
      da->getNodeIndices(indices); 
      unsigned char childNum = da->getChildNumber();
      unsigned char hnMask = da->getHangingNodeIndex(idx);
      unsigned char elemType = 0;
      GET_ETYPE_BLOCK(elemType,hnMask,childNum)
        for(unsigned int j = 0; j < 8; j++) {
          double integral = 0.0;
          //Quadrature Rule
          for(int m = 0; m < numGaussPts; m++) {
            for(int n = 0; n < numGaussPts; n++) {
              for(int p = 0; p < numGaussPts; p++) {
                double xPt = ( (hxOct*(1.0 +gPts[numGaussPts-2][m])*0.5) + x );
                double yPt = ( (hxOct*(1.0 + gPts[numGaussPts-2][n])*0.5) + y );
                double zPt = ( (hxOct*(1.0 + gPts[numGaussPts-2][p])*0.5) + z );
                double solVal = cos(2.0*M_PI*xPt)*cos(2.0*M_PI*yPt)*cos(2.0*M_PI*zPt); 
                double epVal = (1.0 + (lapFac*(
                        (cos(2.0*M_PI*xPt)*cos(2.0*M_PI*xPt)) +
                        (cos(2.0*M_PI*yPt)*cos(2.0*M_PI*yPt)) +
                        (cos(2.0*M_PI*zPt)*cos(2.0*M_PI*zPt)) 
                        )));
                double rhsVal = (solVal + (12.0*(M_PI*M_PI)*solVal*epVal)
                    - (8.0*(M_PI*M_PI)*lapFac*solVal*(
                        (sin(2.0*M_PI*xPt)*sin(2.0*M_PI*xPt)) +
                        (sin(2.0*M_PI*yPt)*sin(2.0*M_PI*yPt)) +
                        (sin(2.0*M_PI*zPt)*sin(2.0*M_PI*zPt)) 
                        )));
                double ShFnVal = ( ShapeFnCoeffs[childNum][elemType][j][0] + 
                    (ShapeFnCoeffs[childNum][elemType][j][1]*gPts[numGaussPts-2][m]) +
                    (ShapeFnCoeffs[childNum][elemType][j][2]*gPts[numGaussPts-2][n]) +
                    (ShapeFnCoeffs[childNum][elemType][j][3]*gPts[numGaussPts-2][p]) +
                    (ShapeFnCoeffs[childNum][elemType][j][4]*gPts[numGaussPts-2][m]*
                     gPts[numGaussPts-2][n]) +
                    (ShapeFnCoeffs[childNum][elemType][j][5]*gPts[numGaussPts-2][n]*
                     gPts[numGaussPts-2][p]) +
                    (ShapeFnCoeffs[childNum][elemType][j][6]*gPts[numGaussPts-2][p]*
                     gPts[numGaussPts-2][m]) +
                    (ShapeFnCoeffs[childNum][elemType][j][7]*gPts[numGaussPts-2][m]*
                     gPts[numGaussPts-2][n]*gPts[numGaussPts-2][p]) );
                integral += (wts[numGaussPts-2][m]*wts[numGaussPts-2][n]
                    *wts[numGaussPts-2][p]*rhsVal*ShFnVal);
              }
            }
          }
          inarray[indices[j]] += (fac*integral);
        }//end for j
    }//end for i
  }//end if active

  da->vecRestoreBuffer(in,inarray,false,false,false,1);

  PetscFunctionReturn(0);
}



/**
  Sets cos(2*pi*x)cos(2*pi*y)cos(2*pi*z) at the node values
  */
PetscErrorCode SetSolution5(ot::DA* da , Vec tmp) {
  PetscFunctionBegin;	 	 
  PetscScalar *inarray;
  VecZeroEntries(tmp);
  da->vecGetBuffer(tmp, inarray, false, false, false, 1);

  unsigned int maxD;
  unsigned int balOctmaxD;

  if(da->iAmActive()) {
    maxD = da->getMaxDepth();
    balOctmaxD = maxD - 1;
    for(da->init<ot::DA_FLAGS::ALL>(); da->curr() < da->end<ot::DA_FLAGS::ALL>(); da->next<ot::DA_FLAGS::ALL>())  
    {
      Point pt;
      pt = da->getCurrentOffset();
      unsigned levelhere = da->getLevel(da->curr()) - 1;
      double hxOct = (double)((double)(1u << (balOctmaxD - levelhere))/(double)(1u << balOctmaxD));
      double x = (double)(pt.xint())/((double)(1u << (maxD-1)));
      double y = (double)(pt.yint())/((double)(1u << (maxD-1)));
      double z = (double)(pt.zint())/((double)(1u << (maxD-1)));

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

      for(int i = 0; i < 8; i++)
      {
        if (!(hn & (1 << i))){
          double xhere, yhere, zhere;
          xhere = x + (coord[i][0]*hxOct);
          yhere = y + (coord[i][1]*hxOct);
          zhere = z + (coord[i][2]*hxOct); 
          double solVal = cos(2.0*M_PI*xhere)*cos(2.0*M_PI*yhere)*cos(2.0*M_PI*zhere); 
          inarray[indices[i]] = solVal;
        }
      }
    }
  }

  da->vecRestoreBuffer(tmp, inarray, false, false, false, 1); 

  PetscFunctionReturn(0);
}

/*
   Measure L-2 error for the variable coefficient scalar
   elliptic homogeneous neumann problem  with RHS given by RHS5
   */
double ComputeError5(ot::DA* da, Vec in)
{

  if(da->iAmActive()) {
    double wts[3] = { (8.0/9.0), (5.0/9.0), (5.0/9.0) };
    double gPts[3] = { 0.0, sqrt((3.0/5.0)), -sqrt((3.0/5.0)) };
    unsigned int maxD = da->getMaxDepth();
    unsigned int balOctmaxD = (maxD - 1);
    double error = 0.0;
    PetscScalar *inarray;
    da->vecGetBuffer(in, inarray, false, false, true, 1);
    da->ReadFromGhostsBegin<PetscScalar>(inarray, 1);
    da->ReadFromGhostsEnd<PetscScalar>(inarray);
    for(da->init<ot::DA_FLAGS::WRITABLE>();
        da->curr() < da->end<ot::DA_FLAGS::WRITABLE>();
        da->next<ot::DA_FLAGS::WRITABLE>())  
    {
      Point pt;
      pt = da->getCurrentOffset();
      unsigned int idx = da->curr();
      unsigned int levelhere = (da->getLevel(idx) - 1);
      double hxOct = (double)((double)(1u << (balOctmaxD - levelhere))/(double)(1u << balOctmaxD));
      double x = (double)(pt.xint())/((double)(1u << (maxD-1)));
      double y = (double)(pt.yint())/((double)(1u << (maxD-1)));
      double z = (double)(pt.zint())/((double)(1u << (maxD-1)));
      double fac = ((hxOct*hxOct*hxOct)/8.0);
      unsigned int indices[8];
      da->getNodeIndices(indices); 
      unsigned char childNum = da->getChildNumber();
      unsigned char hnMask = da->getHangingNodeIndex(idx);
      unsigned char elemType = 0;
      GET_ETYPE_BLOCK(elemType,hnMask,childNum)
        double integral = 0.0;
      //Quadrature Rule
      for(int m = 0; m < 3; m++) {
        for(int n = 0; n < 3; n++) {
          for(int p = 0; p < 3; p++) {
            double xPt = ( (hxOct*(1.0 +gPts[m])*0.5) + x );
            double yPt = ( (hxOct*(1.0 + gPts[n])*0.5) + y );
            double zPt = ( (hxOct*(1.0 + gPts[p])*0.5) + z );
            double fnVal = cos(2.0*M_PI*xPt)*cos(2.0*M_PI*yPt)*cos(2.0*M_PI*zPt); 
            double uhval= 0.0; 
            for(unsigned int j = 0; j < 8; j++) {
              uhval += (inarray[indices[j]]*ShapeFnStencil[childNum][elemType][j][m][n][p]);
            }
            integral += (wts[m]*wts[n]*wts[p]*(fnVal - uhval)*(fnVal - uhval));
          }// end for p
        } // end for n
      }//end for m
      error += integral*fac;
    }//end for i
    da->vecRestoreBuffer(in, inarray, false, false, true, 1);
    double totalError = 0.0;
    par::Mpi_Reduce<double>(&error, &totalError, 1, MPI_SUM, 0, da->getCommActive());
    return(sqrt(totalError));
  } else {
    return(0.0);
  }//end if-else active
}//end fn.

double TestError5(ot::DA* da) {
  Vec tmp;
  da->createVector(tmp, false, false, 1);
  SetSolution5(da, tmp);
  double error = ComputeError5(da, tmp);
  VecDestroy(tmp);
  return error;
}

//Force is a sum of delta functions.
//The location and the magnitude of the delta
//functions is given in the file "deltaSources_<rank>_<npes>.txt"
PetscErrorCode ComputeRHS9(ot::DAMG damg, Vec in) {
  PetscFunctionBegin;	 	 

  ot::DA* da = damg->da;

  MPI_Comm commAll = da->getComm();

  int rankAll = da->getRankAll();
  int npesAll = da->getNpesAll();

  unsigned int dim = da->getDimension();
  unsigned int maxDepth = da->getMaxDepth();
  int npesActive = da->getNpesActive();

  if( npesActive < npesAll ) {
    unsigned int tmpArr[3];
    tmpArr[0] = dim;
    tmpArr[1] = maxDepth;
    tmpArr[2] = npesActive;

    par::Mpi_Bcast<unsigned int>(tmpArr, 3, 0, commAll);

    dim = tmpArr[0];
    maxDepth = tmpArr[1];
    npesActive = static_cast<int>(tmpArr[2]);
  }

  unsigned int balOctMaxD = (maxDepth - 1);

  char fname[250];
  sprintf(fname,"deltaSources_%d_%d.txt", rankAll, npesAll);

  FILE* inFile = fopen(fname,"r");
  if(!inFile) {
    std::cout<<"Unable to open "<<fname<<" for reading."<<std::endl;
  }

  unsigned int numLocalDelta; 
  fscanf(inFile,"%u",&numLocalDelta);

  std::vector<ot::NodeAndValues<double, 4> > tnAndVal(numLocalDelta);

  for(unsigned int i = 0; i < numLocalDelta; i++) {
    double x, y, z, v;
    fscanf(inFile,"%lf",&x);
    fscanf(inFile,"%lf",&y);
    fscanf(inFile,"%lf",&z);
    fscanf(inFile,"%lf",&v);

    assert(x >= 0.0);
    assert(y >= 0.0);
    assert(z >= 0.0);
    assert(x <= 1.0);//Can be equal to 1.0 for domain boundaries
    assert(y <= 1.0);//Can be equal to 1.0 for domain boundaries
    assert(z <= 1.0);//Can be equal to 1.0 for domain boundaries

    unsigned int xint = static_cast<unsigned int>(x*static_cast<double>(1u << balOctMaxD));
    unsigned int yint = static_cast<unsigned int>(y*static_cast<double>(1u << balOctMaxD));
    unsigned int zint = static_cast<unsigned int>(z*static_cast<double>(1u << balOctMaxD));

    tnAndVal[i].node = ot::TreeNode(xint, yint, zint, maxDepth, dim, maxDepth);

    tnAndVal[i].values[0] = x;
    tnAndVal[i].values[1] = y;
    tnAndVal[i].values[2] = z;
    tnAndVal[i].values[3] = v;
  }//end for i

  fclose(inFile);

  std::vector<ot::TreeNode> minBlocks = da->getMinAllBlocks();

  if(npesActive < npesAll) {
    bool* activeStates = new bool[npesAll];

    activeStates[0] = true;
    for(int i = 1; i < npesActive; i++) {
      activeStates[i] = false;
    }
    for(int i = npesActive; i < npesAll; i++) {
      activeStates[i] = true;
    }

    MPI_Comm tmpComm;

    par::splitComm2way(activeStates, &tmpComm, commAll);

    if(rankAll == 0) {
      par::Mpi_Bcast<ot::TreeNode>(&(*(minBlocks.begin())) , npesActive, 0, tmpComm);
    }
    if(rankAll >= npesActive) {
      minBlocks.resize(npesActive);
      par::Mpi_Bcast<ot::TreeNode>(&(*(minBlocks.begin())) , npesActive, 0, tmpComm);
    }

    delete [] activeStates;
  }

  unsigned int* part = new unsigned int[tnAndVal.size()];

  int *sendCnt = new int[npesAll];
  for(int i = 0; i < npesAll; i++) {
    sendCnt[i] = 0;
  }

  for(int i = 0; i < tnAndVal.size(); i++) {
    seq::maxLowerBound<ot::TreeNode>(minBlocks, tnAndVal[i].node, part[i], NULL, NULL);
    sendCnt[part[i]]++; 
  }

  int *recvCnt = new int[npesAll];

  par::Mpi_Alltoall<int>( sendCnt, recvCnt, 1, commAll);

  int *sendOffsets = new int[npesAll];
  int *recvOffsets = new int[npesAll];

  sendOffsets[0] = 0;
  recvOffsets[0] = 0;
  for(int i = 1; i < npesAll; i++) {
    sendOffsets[i] = sendOffsets[i - 1] + sendCnt[i - 1];
    recvOffsets[i] = recvOffsets[i - 1] + recvCnt[i - 1];
  }

  //We can simply communicate the doubles instead, but then we will need to
  //recreate octants from doubles to do further processing.
  std::vector<ot::NodeAndValues<double, 4> > sendList(sendOffsets[npesAll] + sendCnt[npesAll]);
  std::vector<ot::NodeAndValues<double, 4> > recvList(recvOffsets[npesAll] + recvCnt[npesAll]);

  int* tmpSendCnt = new int[npesAll];
  for(int i = 0; i < npesAll; i++) {
    tmpSendCnt[i] = 0;
  }

  for(int i = 0; i < tnAndVal.size(); i++) {
    sendList[sendOffsets[part[i]] + tmpSendCnt[part[i]]] = tnAndVal[i];
    tmpSendCnt[part[i]]++;
  }
  delete [] tmpSendCnt;
  tnAndVal.clear();

  par::Mpi_Alltoallv_sparse<ot::NodeAndValues<double, 4> >( &(*(sendList.begin())), sendCnt,
      sendOffsets, &(*(recvList.begin())), recvCnt, recvOffsets, commAll);

  sendList.clear();

  delete [] part;
  delete [] sendCnt;
  delete [] recvCnt;
  delete [] sendOffsets;
  delete [] recvOffsets;

  sort(recvList.begin(), recvList.end());

  //Now the points are sorted and aligned with the DA partition

  //In PETSc's debug mode they use an Allreduce where all processors (active
  //and inactive) participate. Hence, VecZeroEntries must be called by active
  //and inactive processors.
  VecZeroEntries(in);

  if(!(da->iAmActive())) {
    assert(recvList.empty());
  } else {
    PetscScalar *inarray;
    da->vecGetBuffer(in, inarray, false, false, false, 1);

    unsigned int ptsCtr = 0;

    for(da->init<ot::DA_FLAGS::WRITABLE>();
        da->curr() < da->end<ot::DA_FLAGS::WRITABLE>();
        da->next<ot::DA_FLAGS::WRITABLE>())  
    {
      Point pt;
      pt = da->getCurrentOffset();
      unsigned int idx = da->curr();
      unsigned levelhere = da->getLevel(idx);
      unsigned int xint = pt.xint();
      unsigned int yint = pt.yint();
      unsigned int zint = pt.zint();
      ot::TreeNode currOct(xint, yint, zint, levelhere, dim, maxDepth);
      while( (ptsCtr < recvList.size()) && 
          (recvList[ptsCtr].node < currOct) ) {
        ptsCtr++;
      }
      double hxOct = (double)((double)(1u << (maxDepth - levelhere))/(double)(1u << balOctMaxD));
      double x = static_cast<double>(xint)/((double)(1u << balOctMaxD));
      double y = static_cast<double>(yint)/((double)(1u << balOctMaxD));
      double z = static_cast<double>(zint)/((double)(1u << balOctMaxD));
      unsigned int indices[8];
      da->getNodeIndices(indices); 
      unsigned char childNum = da->getChildNumber();
      unsigned char hnMask = da->getHangingNodeIndex(idx);
      unsigned char elemType = 0;
      GET_ETYPE_BLOCK(elemType,hnMask,childNum)
        while( (ptsCtr < recvList.size()) && 
            ( currOct.isAncestor(recvList[ptsCtr].node) || 
              ( currOct == (recvList[ptsCtr].node) ) ) ) {
          for(unsigned int j = 0; j < 8; j++) {
            double xLoc = (((2.0/hxOct)*(recvList[ptsCtr].values[0] - x)) - 1.0);
            double yLoc = (((2.0/hxOct)*(recvList[ptsCtr].values[1] - y)) - 1.0);
            double zLoc = (((2.0/hxOct)*(recvList[ptsCtr].values[2] - z)) - 1.0);
            double ShFnVal = ( ShapeFnCoeffs[childNum][elemType][j][0] + 
                (ShapeFnCoeffs[childNum][elemType][j][1]*xLoc) +
                (ShapeFnCoeffs[childNum][elemType][j][2]*yLoc) +
                (ShapeFnCoeffs[childNum][elemType][j][3]*zLoc) +
                (ShapeFnCoeffs[childNum][elemType][j][4]*xLoc*yLoc) +
                (ShapeFnCoeffs[childNum][elemType][j][5]*yLoc*zLoc) +
                (ShapeFnCoeffs[childNum][elemType][j][6]*zLoc*xLoc) +
                (ShapeFnCoeffs[childNum][elemType][j][7]*xLoc*yLoc*zLoc) );
            inarray[indices[j]] += ((recvList[ptsCtr].values[3])*ShFnVal);
          }//end for j
          ptsCtr++;
        }
    }//end for i

    da->WriteToGhostsBegin<PetscScalar>(inarray, 1);
    da->WriteToGhostsEnd<PetscScalar>(inarray, 1);

    da->vecRestoreBuffer(in, inarray, false, false, false, 1);

  }//end if active

  PetscFunctionReturn(0);

}//end fn.


