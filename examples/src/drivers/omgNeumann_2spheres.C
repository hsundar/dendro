#include <mpi.h>
#include <cstdio>
#include "oda.h"
#include "omg.h"
#include "Point.h"
#include "parUtils.h"
#include "octUtils.h"
#include "TreeNode.h"
#include "handleStencils.h"
#include <cstdlib>
#include "sys.h"
#include "externVars.h"
#include <sstream>
#include "omgNeumann.h"

//Don't want time to be synchronized. Need to check load imbalance.
#ifdef MPI_WTIME_IS_GLOBAL
#undef MPI_WTIME_IS_GLOBAL
#endif

#ifdef PETSC_USE_LOG

int Jac2DiagEvent;
int Jac2MultEvent;
int Jac2FinestDiagEvent;
int Jac2FinestMultEvent;

#endif

const double w = 1;  // frequency for test purposes
const double c_x =0.5;
const double c_y =0.5;
const double c_z =0.5;
const double inneR = 0.2;
const double outeR = 0.4;
const double penalty = 1e8;

void CalcVarCoef(const std::vector<double> & pts, std::vector<double> & coef)
{
  const double inS = inneR*inneR;
  const double outS = outeR*outeR;
  
  for(int i=0,j=0; i<pts.size(); i+=3,j+=2)
  {
    double s = (pts[i]-c_x)*(pts[i]-c_x) + (pts[i+1]-c_y)*(pts[i+1]-c_y) + (pts[i+2]-c_z)*(pts[i+2]-c_z);   
    
    if ( s>=inS && s<=outS )
    {
      // we are inside the domain
      coef[j]=1;   // weight of -div grad u
      coef[j+1]=0; // weight of u
    }
    else
    {
      // we are outside the domain
      coef[j]=1;
      coef[j+1]=penalty;
    }
  }
}

void CalcVarRHS(const std::vector<double> & pts, std::vector<double> & rhs)
{
  const double inS = inneR*inneR;
  const double outS = outeR*outeR;
  const double l = outeR-inneR;
  
  for(int i=0,j=0; i<pts.size(); i+=3,j++)
  {
   
    double s = (pts[i]-c_x)*(pts[i]-c_x) + (pts[i+1]-c_y)*(pts[i+1]-c_y) + (pts[i+2]-c_z)*(pts[i+2]-c_z);    
    if ( s>=inS && s<=outS )
    {
      // we are inside the domain
      double r = sqrt(s);
      rhs[j]= w*w*M_PI*M_PI/l/l*sin(w*M_PI*(r-inneR)/l)/r; 
    }
    else
    {
      // we are outside the domain
      rhs[j]=0;
    }
  }
}

int main(int argc, char ** argv )
{
  PetscInitialize(&argc,&argv,"options","DENDRO example\n");
  ot::RegisterEvents();

  ot::DAMG_Initialize(MPI_COMM_WORLD);

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if (argc<2) 
  {
    std::cout<<"Usage: mpirun -np <numProcs> "<<argv[0]<<" <filePrefix>"<<std::endl;
    return(-1);
  }

  std::ostringstream fName;
  fName<<argv[1]<<rank<<"_"<<size<<".ot";

  std::vector<ot::TreeNode> octs;
  ot::readNodesFromFile(const_cast<char*>(fName.str().c_str()), octs);
  
  Vec sol;
  ot::DAMG * damg;
  
  // This function pointer will be called while setting up the solver on the coarsest grid if not all processes are active on the coarse grid 
  ot::getPrivateMatricesForKSP_Shell = getPrivateMatricesForKSP_Shell_Jac2;

  solve_neumann_oct(
    /* input parameters: */
     octs,
     CalcVarCoef,
     CalcVarRHS,
     100,

     /* output parameters */
     sol,
     damg
    );
 
  // compare solution with the exact one
  unsigned numMultigridLevels = damg[0]->totalLevels;
  ot::DA* da=damg[numMultigridLevels-1]->da;
  PetscScalar *solArray;
  da->vecGetBuffer(sol,solArray,false,false,true,1);
  double maxDiff = 0;
  unsigned maxDepth=da->getMaxDepth()-1;
  double l = outeR - inneR;

  if(da->iAmActive())
  { 
    double outeS = outeR*outeR;
    double inneS = inneR*inneR;
    
    // first loop over those octants which do not contain ghost nodes.
    for(da->init<ot::DA_FLAGS::INDEPENDENT>(); da->curr() < da->end<ot::DA_FLAGS::INDEPENDENT>(); da->next<ot::DA_FLAGS::INDEPENDENT>())  
    {
      Point pt = da->getCurrentOffset();
      unsigned levelhere = da->getLevel(da->curr()) - 1;
      double hxOct = ldexp(1.0,-levelhere);
      double x = ldexp(pt.x(),-maxDepth );
      double y = ldexp(pt.y(),-maxDepth );
      double z = ldexp(pt.z(),-maxDepth );

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

      double s = (x+hxOct/2-c_x)*(x+hxOct/2-c_x) + (y+hxOct/2-c_y)*(y+hxOct/2-c_y) + (z+hxOct/2-c_z)*(z+hxOct/2-c_z) ;
      if (s<inneS || s>outeS)
	continue; // octant outside the domain, don't account for the error on any of its nodes

      for(int i = 0; i < 8; i++)
	if (!(hn & (1 << i)))
	{
	  double xhere, yhere, zhere;
	  xhere = x + coord[i][0]*hxOct; 
	  yhere = y + coord[i][1]*hxOct; 
	  zhere = z + coord[i][2]*hxOct; 

	  double rhere = sqrt( (xhere-c_x)*(xhere-c_x) + (yhere-c_y)*(yhere-c_y) + (zhere-c_z)*(zhere-c_z) );
	  
	  double currDiff;
	  if (rhere<inneR || rhere>outeR)
	    continue; // solution should be zero here, don't account for error
	  else
	  {
	    // debug
	    // cout<<rhere<<" "<< solArray[indices[i]]<<endl;
	    currDiff = fabs ( solArray[indices[i]] - sin(w*M_PI*(rhere-inneR)/l)/rhere );
	  }
	  
	  if (maxDiff < currDiff)
	    maxDiff=currDiff;
	}
    }  
    
    // now loop over those octants which DO contain ghost nodes; now we need to  check if the node is ghost, and skip it if it is
    size_t myFirstNode = da->getIdxElementBegin();
    size_t postFirstNode = da->getIdxPostGhostBegin();
    for(da->init<ot::DA_FLAGS::DEPENDENT>(); da->curr() < da->end<ot::DA_FLAGS::DEPENDENT>(); da->next<ot::DA_FLAGS::DEPENDENT>())  
    {
      Point pt = da->getCurrentOffset();
      unsigned levelhere = da->getLevel(da->curr()) - 1;
      double hxOct = ldexp(1.0,-levelhere);
      double x = ldexp(pt.x(),-maxDepth );
      double y = ldexp(pt.y(),-maxDepth );
      double z = ldexp(pt.z(),-maxDepth );

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

      double s = (x+hxOct/2-c_x)*(x+hxOct/2-c_x) + (y+hxOct/2-c_y)*(y+hxOct/2-c_y) + (z+hxOct/2-c_z)*(z+hxOct/2-c_z) ;
      if (s<inneS || s>outeS)
	continue; // octant outside the domain, don't account for the error on any of its nodes
      
      for(int i = 0; i < 8; i++)
	if ( indices[i]>=myFirstNode && indices[i]<postFirstNode && !(hn & (1 << i)) )
	{
	  double xhere, yhere, zhere;
	  xhere = x + coord[i][0]*hxOct; 
	  yhere = y + coord[i][1]*hxOct; 
	  zhere = z + coord[i][2]*hxOct; 

	  double rhere = sqrt( (xhere-c_x)*(xhere-c_x) + (yhere-c_y)*(yhere-c_y) + (zhere-c_z)*(zhere-c_z) );
	  
	  double currDiff;
	  if (rhere<inneR || rhere>outeR)
	    continue; // solution should be zero here, don't account for error
	  else
	    currDiff = fabs ( solArray[indices[i]] - sin(w*M_PI*(rhere-inneR)/l)/rhere );
	  
	  if (maxDiff < currDiff)
	    maxDiff=currDiff;
	}
    }  
  }
  da->vecRestoreBuffer(sol,solArray,false,false,true,1); 

  if (!rank)
  {
    double gMaxDiff;
    par::Mpi_Reduce<double>(&maxDiff, &gMaxDiff, 1, MPI_MAX, 0, MPI_COMM_WORLD);  
    std::cout<<"Maximum difference: "<<gMaxDiff<<std::endl;
  }
  else {
    par::Mpi_Reduce<double>(&maxDiff, &maxDiff, 1, MPI_MAX, 0, MPI_COMM_WORLD);  
  }
 
  // destroy multigrid context; this destroys the solution as well 
  DAMGDestroy(damg);
  
  ot::DAMG_Finalize();
  PetscFinalize();
  return 0;
}

