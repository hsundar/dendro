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


void CalcVarCoef(const std::vector<double> & pts, std::vector<double> & coef)
{
  for(int i=0,j=0; i<pts.size(); i+=3,j+=2)
  {
    double s = pts[i]*pts[i] + pts[i+1]*pts[i+1] + pts[i+2]*pts[i+2];
    coef[j] = s+1;
    coef[j+1] = s ;
  }
}

void CalcVarRHS(const std::vector<double> & pts, std::vector<double> & rhs)
{
  for(int i=0,j=0; i<pts.size(); i+=3,j++)
  {
    double s = pts[i]*pts[i] + pts[i+1]*pts[i+1] + pts[i+2]*pts[i+2];
    rhs[j] = (s/*mass_factor*/ + (s+1)*3*w*w*M_PI*M_PI)*cos(w*M_PI*pts[i])*cos(w*M_PI*pts[i+1])*cos(w*M_PI*pts[i+2]) + 2*pts[i]*w*M_PI*sin(w*M_PI*pts[i])*cos(w*M_PI*pts[i+1])*cos(w*M_PI*pts[i+2]) + 2*pts[i+1]*w*M_PI*cos(w*M_PI*pts[i])*sin(w*M_PI*pts[i+1])*cos(w*M_PI*pts[i+2]) + 2*pts[i+2]*w*M_PI*cos(w*M_PI*pts[i])*cos(w*M_PI*pts[i+1])*sin(w*M_PI*pts[i+2]);
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
  fName<<argv[1]<<rank<<"_"<<size<<".pts";

  std::vector<double> pts;
  ot::readPtsFromFile(const_cast<char*>(fName.str().c_str()), pts);
  
  Vec sol;
  ot::DAMG * damg;
  
  // This function pointer will be called while setting up the solver on the coarsest grid if not all processes are active on the coarse grid 
  ot::getPrivateMatricesForKSP_Shell = getPrivateMatricesForKSP_Shell_Jac2;

  solve_neumann(
    /* input parameters: */
     pts,
     CalcVarCoef,
     CalcVarRHS,
     30, /* upper bound */

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

  if(da->iAmActive())
  { 
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

      for(int i = 0; i < 8; i++)
	if (!(hn & (1 << i)))
	{
	  double xhere, yhere, zhere;
	  xhere = x + coord[i][0]*hxOct; 
	  yhere = y + coord[i][1]*hxOct; 
	  zhere = z + coord[i][2]*hxOct; 

	  double currDiff = fabs ( solArray[indices[i]] - cos(w*M_PI*xhere)*cos(w*M_PI*yhere)*cos(w*M_PI*zhere) );	
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

      for(int i = 0; i < 8; i++)
	if ( indices[i]>=myFirstNode && indices[i]<postFirstNode && !(hn & (1 << i)) )
	{
	  double xhere, yhere, zhere;
	  xhere = x + coord[i][0]*hxOct; 
	  yhere = y + coord[i][1]*hxOct; 
	  zhere = z + coord[i][2]*hxOct; 

	  double currDiff = fabs ( solArray[indices[i]] - cos(w*M_PI*xhere)*cos(w*M_PI*yhere)*cos(w*M_PI*zhere) );	
	  if (maxDiff < currDiff)
	     maxDiff=currDiff;
	}
    }  
  }
  da->vecRestoreBuffer(sol,solArray,false,false,true,1); 

  double gMaxDiff;
  par::Mpi_Reduce<double>(&maxDiff, &gMaxDiff, 1, MPI_MAX, 0, MPI_COMM_WORLD);  
  if (!rank)
  {
    std::cout<<"Maximum difference: "<<gMaxDiff<<std::endl;
  }
    
  double solmax,solmin;
  VecMax(sol,NULL,&solmax);
  VecMin(sol,NULL,&solmin);

  if (!rank)
    std::cout<<"solution min: " << solmin << "; solution max: " << solmax << std::endl;

  // destroy multigrid context; this destroys the solution as well 
  DAMGDestroy(damg);
  
  ot::DAMG_Finalize();
  PetscFinalize();
  return 0;
}

