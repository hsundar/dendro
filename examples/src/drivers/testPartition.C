#include "mpi.h"
#include "petsc.h"
#include "sys.h"

#include "parUtils.h"
#include "TreeNode.h"
#include "colors.h"
#include "oda.h"

#include <cstdlib>
#include <execinfo.h>
#include <unistd.h>
#include <cxxabi.h>

#include "externVars.h"
#include "dendro.h"
#include "rotation.h"
#include "treenode2vtk.h"


//Don't want time to be synchronized. Need to check load imbalance.
#ifdef MPI_WTIME_IS_GLOBAL
#undef MPI_WTIME_IS_GLOBAL
#endif




//@author: Milinda Fernando.
// This function is to calculate the boundary faces
// Assume that the given octree vector is sorted.
int calculateBoundaryFaces(const std::vector<ot::TreeNode> & mesh, int q) {

  ot::TreeNode first=mesh[0];
  ot::TreeNode last=mesh[mesh.size()-1];
  ot::TreeNode R (first.getDim(), first.getMaxDepth());


  int com_size=q;
  int mesh_nodes=mesh.size();
  assert(mesh_nodes>q);
  int local_mesh_size=mesh_nodes/q;
  int boundary_faces[q];



  int found_pt;
  int num_boundary_faces=0;

  ot::TreeNode top;
  ot::TreeNode bottom;
  ot::TreeNode left;
  ot::TreeNode right;
  ot::TreeNode front;
  ot::TreeNode back;
  int temp=0;
  int begin=0;
  int end=0;
  int total_boundary_faces=0;

  for (int j=0;j<q;j++){

     boundary_faces[j]=0;
     begin=j*local_mesh_size;

     temp=begin;
     end=begin + local_mesh_size;

     if(end+local_mesh_size>mesh_nodes)
       end=mesh_nodes;


     num_boundary_faces=0;

      while(temp<end)//for(int i=0;i<mesh.size();i++)
      {

      top=mesh[temp].getTop();
      bottom=mesh[temp].getBottom();
      left=mesh[temp].getLeft();
      right=mesh[temp].getRight();
      front=mesh[temp].getFront();
      back=mesh[temp].getBack();

      found_pt=(std::lower_bound(&mesh[begin], &mesh[end-1], top, std::less<ot::TreeNode>()) - &mesh[begin]);
       // std::cout<<"top:"<<found_pt<<std::endl;
      if(found_pt==0 || ((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(top))))
      {
//        if((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(top)))
//          std::cout<<"end condition"<<std::endl;

        num_boundary_faces++;
      }

      found_pt=(std::lower_bound(&mesh[begin], &mesh[end-1], bottom, std::less<ot::TreeNode>()) - &mesh[begin]);
        //std::cout<<"bottom:"<<found_pt<<std::endl;
      if(found_pt==0 || ((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(bottom))))
      {
//        if((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(top)))
//          std::cout<<"end condition"<<std::endl;


        num_boundary_faces++;
      }

      found_pt=(std::lower_bound(&mesh[begin], &mesh[end-1], left, std::less<ot::TreeNode>()) - &mesh[begin]);
        //std::cout<<"left:"<<found_pt<<std::endl;
      if(found_pt==0 || ((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(left))))
      {
//        if((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(top)))
//          std::cout<<"end condition"<<std::endl;

        num_boundary_faces++;
      }

      found_pt=(std::lower_bound(&mesh[begin], &mesh[end-1], right, std::less<ot::TreeNode>()) - &mesh[begin]);
        //std::cout<<"right:"<<found_pt<<std::endl;
      if(found_pt==0 || ((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(right))))
      {
//        if((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(top)))
//          std::cout<<"end condition"<<std::endl;

        num_boundary_faces++;
      }
      found_pt=(std::lower_bound(&mesh[begin], &mesh[end-1], front, std::less<ot::TreeNode>()) - &mesh[begin]);
        //std::cout<<"front:"<<found_pt<<std::endl;
      if(found_pt==0 || ((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(front))))
      {
//        if((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(top)))
//          std::cout<<"end condition"<<std::endl;

        num_boundary_faces++;
      }

      found_pt=(std::lower_bound(&mesh[begin], &mesh[end-1], back, std::less<ot::TreeNode>()) - &mesh[begin]);
        //std::cout<<"back:"<<found_pt<<std::endl;
      if(found_pt==0 || ((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(back)))) {

//        if((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(top)))
//          std::cout<<"end condition"<<std::endl;

        num_boundary_faces++;
      }


      temp++;
      if(temp%local_mesh_size==0)
      {
        break;
      }

    }
    boundary_faces[j]=num_boundary_faces;
    total_boundary_faces=total_boundary_faces+num_boundary_faces;
    //std::cout<<"q:"<<j<<"\t "<<"number of boundary faces:"<<boundary_faces[j]<<std::endl;

  }


  //MPI_Reduce(boundary_faces,&total_boundary_faces,q, MPI_INT,MPI_SUM,0, MPI_COMM_WORLD);
  return total_boundary_faces;


}






/** Print a demangled stack backtrace of the caller function to FILE* out. */
void handler (int sig) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // char fname[256];
  // sprintf(fname, "trace%.2d", rank);
  // FILE *out = fopen(fname, "w");
  unsigned int max_frames = 63;

  // if (!rank) {
    printf("%s---------------------------------%s\n", RED, NRM);
    printf("%sError:%s signal %d:\n", RED, NRM, sig);
    printf("%s---------------------------------%s\n", RED, NRM);
    printf("\n%s======= stack trace =======%s\n", GRN, NRM);
  // }

  // fprintf(out, "======= stack trace =======\n");

  // storage array for stack trace address data
  void *addrlist[max_frames + 1];

  // retrieve current stack addresses
  int addrlen = backtrace(addrlist, sizeof(addrlist) / sizeof(void *));

  if (addrlen == 0) {
    // if (!rank)
    fprintf(stderr, "%s  <empty, possibly corrupt>%s\n",RED, NRM);

    // fprintf(out, "    <empty, possibly corrupt>\n");
    return;
  }

  // resolve addresses into strings containing "filename(function+address)",
  // this array must be free()-ed
  char **symbollist = backtrace_symbols(addrlist, addrlen);

  // allocate string which will be filled with the demangled function name
  size_t funcnamesize = 256;
  char *funcname = (char *) malloc(funcnamesize);

  // iterate over the returned symbol lines. skip the first, it is the
  // address of this function.
  for (int i = 1; i < addrlen; i++) {
    char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

    // find parentheses and +address offset surrounding the mangled name:
    // ./module(function+0x15c) [0x8048a6d]
    for (char *p = symbollist[i]; *p; ++p) {
      if (*p == '(')
        begin_name = p;
      else if (*p == '+')
        begin_offset = p;
      else if (*p == ')' && begin_offset) {
        end_offset = p;
        break;
      }
    }

    if (begin_name && begin_offset && end_offset
        && begin_name < begin_offset) {
      *begin_name++ = '\0';
      *begin_offset++ = '\0';
      *end_offset = '\0';

      // mangled name is now in [begin_name, begin_offset) and caller
      // offset in [begin_offset, end_offset). now apply
      // __cxa_demangle():

      int status;
      char *ret = abi::__cxa_demangle(begin_name,
                                      funcname, &funcnamesize, &status);
      if (status == 0) {
        funcname = ret; // use possibly realloc()-ed string
        // if (!rank)
        printf("%s[%.2d]%s%s : %s%s%s : \n",RED,rank,YLW, symbollist[i], MAG, funcname,NRM);

        // fprintf(out, "%s : %s : ", symbollist[i], funcname);
      }
      else {
        // demangling failed. Output function name as a C function with
        // no arguments.
        // if (!rank)
        printf("%s[%.2d]%s%s : %s%s()%s : \n", RED,rank, YLW, symbollist[i], GRN,begin_name, NRM);

        // fprintf(out, "%s : %s() : ", symbollist[i], begin_name);
      }
      size_t p = 0;
      char syscom[256];
      while(symbollist[i][p] != '(' && symbollist[i][p] != ' ' && symbollist[i][p] != 0)
        ++p;

      sprintf(syscom,"addr2line %p -e %.*s", addrlist[i], p, symbollist[i]);
      //last parameter is the file name of the symbol
      system(syscom);
    }
    else {
      // couldn't parse the line? print the whole line.
      // if (!rank)
      printf("%sCouldn't Parse:%s  %s\n", RED, NRM, symbollist[i]);

      // fprintf(out, "  %s\n", symbollist[i]);
    }
  }

  free(funcname);
  free(symbollist);
  // fclose(out);

  exit(1);
}

int main(int argc, char **argv) {

  int size, rank;
  bool incCorner = 1;
  // char bFile[50];
  char ptsFileName[256];
  bool compressLut = false;

  unsigned int ptsLen;
  unsigned int maxNumPts = 1;
  unsigned int dim = 3;
  unsigned int maxDepth = 8;
  double gSize[3];
  //initializeHilbetTable(2);
  G_MAX_DEPTH = maxDepth;
  G_dim = dim;

  _InitializeHcurve();

  double localTime, totalTime;
  double startTime, endTime;
  DendroIntL localSz, totalSz;
  std::vector<double> pts;
  std::vector<ot::TreeNode> linOct, balOct;

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " inpfile " << std::endl;
    return -1;
  }

  PetscInitialize(&argc, &argv, "options.hs", NULL);
  ot::RegisterEvents();
  ot::DA_Initialize(MPI_COMM_WORLD);

  signal(SIGSEGV, handler);   // install our handler
  signal(SIGTERM, handler);   // install our handler

#ifdef PETSC_USE_LOG
  int stages[3];
  PetscLogStageRegister("P2O.", &stages[0]);
  PetscLogStageRegister("Bal", &stages[1]);
  PetscLogStageRegister("ODACreate", &stages[2]);
#endif

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  sprintf(ptsFileName, "%s%d_%d.pts", argv[1], rank, size);

  //Read pts from files
  if (!rank) {
    std::cout << RED " Reading  " << argv[1] << NRM << std::endl; // Point size
  }
  ot::readPtsFromFile(ptsFileName, pts);

  if (!rank) {
    std::cout << GRN " Finished reading  " << argv[1] << NRM << std::endl; // Point size
  }

  ptsLen = pts.size();

  std::vector<ot::TreeNode> tmpNodes;
  for (int i = 0; i < ptsLen; i += 3) {
    if ((pts[i] > 0.0) &&
        (pts[i + 1] > 0.0)
        && (pts[i + 2] > 0.0) &&
        (((unsigned int) (pts[i] * ((double) (1u << maxDepth)))) < (1u << maxDepth)) &&
        (((unsigned int) (pts[i + 1] * ((double) (1u << maxDepth)))) < (1u << maxDepth)) &&
        (((unsigned int) (pts[i + 2] * ((double) (1u << maxDepth)))) < (1u << maxDepth))) {
      tmpNodes.push_back(ot::TreeNode((unsigned int) (pts[i] * (double) (1u << maxDepth)),
                                      (unsigned int) (pts[i + 1] * (double) (1u << maxDepth)),
                                      (unsigned int) (pts[i + 2] * (double) (1u << maxDepth)),
                                      maxDepth, dim, maxDepth));
    }
  }
  pts.clear();

  // treeNodesTovtk(tmpNodes, rank, "input_points");

  // std::cout << rank << "removeDuplicates" << std::endl;
  par::removeDuplicates<ot::TreeNode>(tmpNodes, false, MPI_COMM_WORLD);

  linOct = tmpNodes;
  tmpNodes.clear();
  // std::cout << rank << "partition" << std::endl;
  par::partitionW<ot::TreeNode>(linOct, NULL, MPI_COMM_WORLD);

  // treeNodesTovtk(linOct,rank,"par_1");
  treeNodesTovtk(linOct, rank, "input_points");

  // reduce and only print the total ...
  localSz = linOct.size();
  par::Mpi_Reduce<DendroIntL>(&localSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);

  pts.resize(3 * (linOct.size()));
  ptsLen = (3 * (linOct.size()));
  for (int i = 0; i < linOct.size(); i++) {
    pts[3 * i] = (((double) (linOct[i].getX())) + 0.5) / ((double) (1u << maxDepth));
    pts[(3 * i) + 1] = (((double) (linOct[i].getY())) + 0.5) / ((double) (1u << maxDepth));
    pts[(3 * i) + 2] = (((double) (linOct[i].getZ())) + 0.5) / ((double) (1u << maxDepth));
  } //end for i
  linOct.clear();
  gSize[0] = 1.0;
  gSize[1] = 1.0;
  gSize[2] = 1.0;

  // ========== Points2Octree ===========
  MPI_Barrier(MPI_COMM_WORLD);
  if (!rank) {
    std::cout << BLU << "===============================================" << NRM << std::endl;
    std::cout << RED " Starting Points to Octree" NRM << std::endl;
    std::cout << BLU << "===============================================" << NRM << std::endl;
  }

#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[0]);
#endif

  startTime = MPI_Wtime();
  ot::points2Octree(pts, gSize, linOct, dim, maxDepth, maxNumPts, MPI_COMM_WORLD);
  endTime = MPI_Wtime();

  par::sampleSort(linOct, balOct, MPI_COMM_WORLD);
  linOct = balOct;
  balOct.clear();

#ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif

  if (!rank) {
    std::cout << BLU << "===============================================" << NRM << std::endl;
    std::cout << RED " Ended Points to Octree" NRM << std::endl;
    std::cout << BLU << "===============================================" << NRM << std::endl;
  }

  // par::partitionW<ot::TreeNode>(linOct, NULL, MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  treeNodesTovtk(linOct, rank, "bfBalancing");

  DendroIntL num_surface=calculateBoundaryFaces(linOct,10);
  //if (!rank) {
    std::cout << BLU << "===============================================" << NRM << std::endl;
    std::cout << RED " Number of Boundary Faces: Rank:" <<rank<<"\t"<<num_surface<<NRM << std::endl;
    std::cout << BLU << "===============================================" << NRM << std::endl;
  //}
  DendroIntL total_boundary_faces=0;
  par::Mpi_Reduce<DendroIntL>(&num_surface,&total_boundary_faces,1, MPI_SUM, 0, MPI_COMM_WORLD);




  localTime = endTime - startTime;
  par::Mpi_Reduce<double>(&localTime, &totalTime, 1, MPI_MAX, 0, MPI_COMM_WORLD);
  if (!rank) {
    std::cout << GRN " P2n Time: " YLW << totalTime << NRM << std::endl;
  }
  // reduce and only print the total ...
  localSz = linOct.size();
  par::Mpi_Reduce<DendroIntL>(&localSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    std::cout << GRN " # of Unbalanced Octants: " YLW << totalSz << NRM << std::endl;
  }

  if (!rank) {
    double surf_volume_r=total_boundary_faces/(double)totalSz;
    std::cout << BLU << "===============================================" << NRM << std::endl;
    std::cout << RED " Surface to Volume ratio:" <<total_boundary_faces<<"/"<<totalSz<<"="<<surf_volume_r<<NRM << std::endl;
    std::cout << BLU << "===============================================" << NRM << std::endl;
  }

  pts.clear();

  //treeNodesTovtk(linOct, rank, "p2o_output");

  // =========== Balancing ============
  MPI_Barrier(MPI_COMM_WORLD);
  if (!rank) {
    std::cout << BLU << "===============================================" << NRM << std::endl;
    std::cout << RED " Starting 2:1 Balance" NRM << std::endl;
    std::cout << BLU << "===============================================" << NRM << std::endl;
  }
#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[1]);
#endif


  // std::cout << "input[0] is " << linOct[0] << std::endl;
  startTime = MPI_Wtime();
  ot::balanceOctree(linOct, balOct, dim, maxDepth, incCorner, MPI_COMM_WORLD, NULL, NULL);
  endTime = MPI_Wtime();

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //treeNodesTovtk(balOct,rank,"afBalancing");

#ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif
  linOct.clear();
  // compute total inp size and output size
  localSz = balOct.size();
  localTime = endTime - startTime;
  par::Mpi_Reduce<DendroIntL>(&localSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  par::Mpi_Reduce<double>(&localTime, &totalTime, 1, MPI_MAX, 0, MPI_COMM_WORLD);

  if (!rank) {
    std::cout << "# of Balanced Octants: " << totalSz << std::endl;
    std::cout << "bal Time: " << totalTime << std::endl;
  }

  treeNodesTovtk(balOct, rank, "bal_output");


  //==================ODA Meshing=================================
  if (!rank) {
    std::cout << BLU << "===============================================" << NRM << std::endl;
    std::cout << RED " Starting ODA Meshing" NRM << std::endl;
    std::cout << BLU << "===============================================" << NRM << std::endl;
  }
  //ODA ...
  MPI_Barrier(MPI_COMM_WORLD);
#ifdef PETSC_USE_LOG
  PetscLogStagePush(stages[2]);
#endif
  startTime = MPI_Wtime();
  assert(!(balOct.empty()));
  ot::DA da(balOct, MPI_COMM_WORLD, MPI_COMM_WORLD, compressLut);
  endTime = MPI_Wtime();
#ifdef PETSC_USE_LOG
  PetscLogStagePop();
#endif
  balOct.clear();
  // compute total inp size and output size
  localSz = da.getNodeSize();
  localTime = endTime - startTime;
  par::Mpi_Reduce<DendroIntL>(&localSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
  par::Mpi_Reduce<double>(&localTime, &totalTime, 1, MPI_MAX, 0, MPI_COMM_WORLD);

  if (!rank) {
    std::cout << "Total # Vertices: " << totalSz << std::endl;
    std::cout << "Time to build ODA: " << totalTime << std::endl;
  }

  //! Quality of the partition ...
  DendroIntL maxNodeSize, minNodeSize,
      maxBdyNode, minBdyNode,
      maxIndepSize, minIndepSize,
      maxElementSize, minElementSize;

  localSz = da.getNodeSize();
  par::Mpi_Reduce<DendroIntL>(&localSz, &maxNodeSize, 1, MPI_MAX, 0, MPI_COMM_WORLD);
  par::Mpi_Reduce<DendroIntL>(&localSz, &minNodeSize, 1, MPI_MIN, 0, MPI_COMM_WORLD);

  localSz = da.getBoundaryNodeSize();
  par::Mpi_Reduce<DendroIntL>(&localSz, &maxBdyNode, 1, MPI_MAX, 0, MPI_COMM_WORLD);
  par::Mpi_Reduce<DendroIntL>(&localSz, &minBdyNode, 1, MPI_MIN, 0, MPI_COMM_WORLD);

  localSz = da.getElementSize();
  par::Mpi_Reduce<DendroIntL>(&localSz, &maxElementSize, 1, MPI_MAX, 0, MPI_COMM_WORLD);
  par::Mpi_Reduce<DendroIntL>(&localSz, &minElementSize, 1, MPI_MIN, 0, MPI_COMM_WORLD);

  localSz = da.getIndependentSize();
  par::Mpi_Reduce<DendroIntL>(&localSz, &maxIndepSize, 1, MPI_MAX, 0, MPI_COMM_WORLD);
  par::Mpi_Reduce<DendroIntL>(&localSz, &minIndepSize, 1, MPI_MIN, 0, MPI_COMM_WORLD);

  if (!rank) {
    std::cout << "Nodes          \t(" << minNodeSize << ", " << maxNodeSize << ")" << std::endl;
    std::cout << "Boundary Node  \t(" << minBdyNode << ", " << maxBdyNode << ")" << std::endl;
    std::cout << "Element        \t(" << minElementSize << ", " << maxElementSize << ")" << std::endl;
    std::cout << "Independent    \t(" << minIndepSize << ", " << maxIndepSize << ")" << std::endl;
  }

  //! ========================

  if (!rank) {
    std::cout << GRN << "Finalizing ..." << NRM << std::endl;
  }

  ot::DA_Finalize();
  PetscFinalize();
} //end function

