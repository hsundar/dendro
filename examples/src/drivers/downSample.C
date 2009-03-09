
#include "mpi.h"
#include "petsc.h"
#include "sys.h"
#include "TreeNode.h"
#include "parUtils.h"
#include "omg.h"
#include "oda.h"
#include <cstdlib>
#include "colors.h"
#include "externVars.h"
#include "dendro.h"

static char help[] = "Solves 3D PBE using FEM, MG, DA and Matrix-Free methods.\n\n";
//Don't want time to be synchronized. Need to check load imbalance.
#ifdef MPI_WTIME_IS_GLOBAL
#undef MPI_WTIME_IS_GLOBAL
#endif

int main(int argc, char ** argv ) {	
  int size, rank;
  unsigned int numPts;
  char pFile[50],oFile[50];
  double gSize[3];
  unsigned int ptsLen;
  unsigned int maxNumPts= 1;
  unsigned int dim=3;
  unsigned int maxDepth=30;
  DendroIntL locSz, totalSz;
  std::vector<ot::TreeNode> linOct, balOct;
  std::vector<double> pts;

  PetscInitialize(&argc,&argv,"optionsDS",help);
  ot::RegisterEvents();

  ot::DAMG_Initialize(MPI_COMM_WORLD);

#ifdef PETSC_USE_LOG
  int stages[3];
  PetscLogStageRegister(&stages[0],"P2O.");
  PetscLogStageRegister(&stages[1],"Bal");  
  PetscLogStageRegister(&stages[2],"Solve");  
#endif

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(argc < 2) {
    std::cerr << "Usage: " << argv[0] << "inpfile  maxDepth[30] \
      dim[3] maxNumPtsPerOctant[1] " << std::endl;
    return -1;
  }
  if(argc > 2) {
    maxDepth = atoi(argv[2]);
  }
  if(argc > 3) {
    dim = atoi(argv[3]);
  }
  if(argc > 4) {
    maxNumPts = atoi(argv[4]);
  }

  sprintf(pFile,"%s%d_%d.pts",argv[1],rank,size);
  sprintf(oFile,"%s_Out%d_%d.pts",argv[1],rank,size);

  //Points2Octree....
  MPI_Barrier(MPI_COMM_WORLD);	
  if(!rank){
    std::cout << " reading  "<<pFile<<std::endl; // Point size
  }
  ot::readPtsFromFile(pFile, pts);
  if(!rank){
    std::cout << " finished reading  "<<pFile<<std::endl; // Point size
  }
  MPI_Barrier(MPI_COMM_WORLD);	
  ptsLen = pts.size();
  std::vector<ot::TreeNode> tmpNodes;
  for(unsigned int i = 0; i < ptsLen; i += 3) {
    if( (pts[i] > 0.0) &&
        (pts[i+1] > 0.0)  
        && (pts[i+2] > 0.0) &&
        ( (static_cast<unsigned int>(pts[i]*(static_cast<double>(1u<<maxDepth)))) < (1u<<maxDepth))  &&
        ( (static_cast<unsigned int>(pts[i+1]*(static_cast<double>(1u<<maxDepth)))) < (1u<<maxDepth))  &&
        ( (static_cast<unsigned int>(pts[i+2]*(static_cast<double>(1u<<maxDepth)))) < (1u<<maxDepth)) ) {
      tmpNodes.push_back( ot::TreeNode(
            static_cast<unsigned int>(pts[i]*static_cast<double>(1u<<maxDepth)),
            static_cast<unsigned int>(pts[i+1]*static_cast<double>(1u<<maxDepth)),
            static_cast<unsigned int>(pts[i+2]*static_cast<double>(1u<<maxDepth)),
            maxDepth,dim,maxDepth) );
    }
  }
  pts.clear();
  par::removeDuplicates<ot::TreeNode>(tmpNodes,false,MPI_COMM_WORLD);	
  linOct = tmpNodes;
  tmpNodes.clear();
  par::partitionW<ot::TreeNode>(linOct, NULL,MPI_COMM_WORLD);
  // reduce and only print the total ...
  locSz = linOct.size();
  par::Mpi_Reduce<DendroIntL>(&locSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);	
  if(rank==0) {
    std::cout<<"# pts= " << totalSz<<std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);	

  pts.resize(3*(linOct.size()));
  ptsLen = (3*(linOct.size()));
  for(unsigned int i = 0; i < linOct.size(); i++) {
    pts[3*i] = (((double)(linOct[i].getX())) + 0.5)/((double)(1u << maxDepth));
    pts[(3*i)+1] = (((double)(linOct[i].getY())) +0.5)/((double)(1u << maxDepth));
    pts[(3*i)+2] = (((double)(linOct[i].getZ())) +0.5)/((double)(1u << maxDepth));
  }//end for i
  linOct.clear();
  gSize[0] = 1.;
  gSize[1] = 1.;
  gSize[2] = 1.;

  MPI_Barrier(MPI_COMM_WORLD);	
#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[0]);
#endif
  ot::points2Octree(pts, gSize, linOct, dim, maxDepth, maxNumPts, MPI_COMM_WORLD);
#ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif
  // reduce and only print the total ...
  locSz = linOct.size();
  par::Mpi_Reduce<DendroIntL>(&locSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  if(rank==0) {
    std::cout<<"# of Unbalanced Octants: " << totalSz<<std::endl;
  }
  pts.clear();

  //Balancing...
  MPI_Barrier(MPI_COMM_WORLD);	
#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[1]);
#endif
  ot::balanceOctree (linOct, balOct, dim, maxDepth, true, MPI_COMM_WORLD, NULL, NULL);
#ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif
  linOct.clear();
  // compute total inp size and output size
  locSz = balOct.size();
  par::Mpi_Reduce<DendroIntL>(&locSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);

  if(!rank) {
    std::cout << "# of Balanced Octants: "<< totalSz << std::endl;       
  }

  MPI_Barrier(MPI_COMM_WORLD);	
#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[2]);
#endif


  ot::DAMG        *damg;    
  int       nlevels = 1; //number of multigrid levels
  unsigned int       dof =1;// degrees of freedom per node  

  PetscInt nlevelsPetscInt = nlevels;
  PetscOptionsGetInt(0, "-nlevels", &nlevelsPetscInt, 0);
  nlevels = nlevelsPetscInt;

  // Note: The user context for all levels will be set separately later.
  MPI_Barrier(MPI_COMM_WORLD);	

  if(!rank) {
    std::cout<<"nlevels initial: "<<nlevels<<std::endl;
  }

  ot::DAMGCreateAndSetDA(PETSC_COMM_WORLD, nlevels, NULL, &damg,
      balOct, dof, 1.5, false, true);

  if(!rank) {
    std::cout<<"nlevels final: "<<nlevels<<std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);	

  if(!rank) {
    std::cout << "Created DA for all levels."<< std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  ot::PrintDAMG(damg);
  MPI_Barrier(MPI_COMM_WORLD);

  ot::DA* dac = damg[0]->da;

  std::vector<double> outputPts;

  if(dac->iAmActive()) {
    double factor = (1.0/(static_cast<double>(1u << maxDepth)));
    for(dac->init<ot::DA_FLAGS::WRITABLE>();
        dac->curr() < dac->end<ot::DA_FLAGS::WRITABLE>();
        dac->next<ot::DA_FLAGS::WRITABLE>() ) {
      Point currPt = dac->getCurrentOffset();
      outputPts.push_back(currPt.x()*factor);
      outputPts.push_back(currPt.y()*factor);
      outputPts.push_back(currPt.z()*factor);
    } 
  }

  ot::writePtsToFile(oFile, outputPts);
  outputPts.clear();

  DAMGDestroy(damg);

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

