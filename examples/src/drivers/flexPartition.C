//
// Created by hari on 10/14/15.
//

#include "mpi.h"
#include "petsc.h"

#include "parUtils.h"
#include "TreeNode.h"
#include "colors.h"

#include <cstdlib>
#include <sys.h>

#include "externVars.h"
#include "dendro.h"
#include "rotation.h"
#include "treenode2vtk.h"
#include "../../include/octreeStatistics.h"


#include <execinfo.h>

#include <cxxabi.h>

#include "externVars.h"

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



void flexiblePartitionCalculation(std::vector<ot::TreeNode>& balOct,double slack,int q)
{

  int rank,size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  unsigned long localSz=balOct.size(), globalSz;
  MPI_Allreduce(&localSz, &globalSz, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

  int slackCnt = slack*globalSz/size/q;

  if(!rank) std::cout << "slack size is " << slackCnt << " octants" << std::endl;

  //----------------------------------------------------------------------
  //   FLEX
  //----------------------------------------------------------------------

  // 1. each process sends slackCnt octs to next/prev.
  std::vector<ot::TreeNode> slack_next(slackCnt), slack_prev(slackCnt);
  MPI_Status status;

  int prev = (rank)?rank-1:MPI_PROC_NULL;
  int next = (rank+1==size)?MPI_PROC_NULL:rank+1;

  ot::TreeNode* sendPrev = &(*(balOct.begin()));
  ot::TreeNode* sendNext = &(balOct[balOct.size() - slackCnt]);

  ot::TreeNode* recvPrev = &(*(slack_prev.begin()));
  ot::TreeNode* recvNext = &(*(slack_next.begin()));

  if(!rank)
    std::cout<<"TreeNode Communication started"<<std::endl;

  //send to prev

  if (rank)
    par::Mpi_Sendrecv<ot::TreeNode>(sendPrev, slackCnt, prev, 0,
                                    recvPrev, slackCnt, prev, 0,
                                    MPI_COMM_WORLD, &status);
  // send to next
  if (rank != (size-1))
    par::Mpi_Sendrecv<ot::TreeNode>(sendNext, slackCnt, next, 0,
                                    recvNext, slackCnt, next, 0,
                                    MPI_COMM_WORLD, &status);

//  for(int i=0;i<slack_prev.size();i++)
//  std::cout<<"Recive Previous:"<<i<<" ::"<<slack_prev[i]<<std::endl;


  // Have the extra octants ...
  // 2. compute partitions.
  int part_size = localSz/q;
  int* partitions   = new int[q+1];
  int* num_faces    = new int[q+1];
  int faces;

  int faces_all[size];
  if(!rank)
    std::cout<<"TreeNode Communication Completed"<<std::endl;

  // partition 0
  partitions[0] = slackCnt;
  if (!rank) {
    num_faces[0] = calculateBoundaryFaces(balOct.begin(), balOct.begin()+part_size);
    slack_prev.insert(slack_prev.end(), balOct.begin(), balOct.begin()+part_size);
    for (int j = 0; j < slackCnt; ++j) {
      assert(j<slack_prev.size());
      faces = calculateBoundaryFaces(slack_prev.begin()+j, slack_prev.end());
      if (faces < num_faces[0]) {
        num_faces[0] = faces;
        partitions[0] = j;
      }
    }
  }
  std::vector<ot::TreeNode>::const_iterator first, last;
  // local partitions
  for (int i=1; i<q; ++i) {
    first = balOct.begin()+i*part_size;
    last = (i==(q-1))?balOct.end():balOct.begin()+(i+1)*part_size;
    assert(first<last);
    partitions[i] = i*part_size + slackCnt;
    num_faces[i] = calculateBoundaryFaces(first,last);
    first -= slackCnt;
    for (int j = 0; j < slackCnt; ++j) {
      assert(j<balOct.size());
      faces = calculateBoundaryFaces(first+j, last);
      //std::cout<<"Faces:"<<faces<<std::endl;
      if (faces < num_faces[i]) {
        num_faces[i] = faces;
        partitions[i] = j;
      }

    }
  }

  if(!rank)
    std::cout<<"Flexible part calculation completed"<<std::endl;
  //std::cout<<"Rank:"<<rank<<" Faces:"<<faces<<std::endl;

  int min_faces;
  int max_faces;
  int faces_sum;

  double mean_num_faces;

  MPI_Allreduce(&faces,&min_faces,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&faces,&max_faces,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&faces,&faces_sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  mean_num_faces=(double)faces_sum/size;

  MPI_Allgather(&faces,1,MPI_INT,faces_all,1,MPI_INT,MPI_COMM_WORLD);

//  if(!rank)
//    for(int i=0;i<size;i++)
//      std::cout<<"faces_all:"<<i<<" :"<<faces_all[i]<<std::endl;


  if (!rank) {
    std::cout << BLU << "===============================================" << NRM << std::endl;
    std::cout << RED " Boundary Surfaces (min):"<<min_faces<< NRM << std::endl;
    std::cout << RED " Boundary Surfaces (max):"<<max_faces<< NRM << std::endl;
    std::cout << RED " Boundary Surfaces (mean):"<<mean_num_faces<< NRM << std::endl;
    std::cout << BLU << "===============================================" << NRM << std::endl;
  }



  // partition q
  partitions[q] = balOct.size();

  //----------------------------------------------------------------------
  if (!rank) {
    std::cout << GRN << "Finalizing ..." << NRM << std::endl;
  }

  delete [] partitions;
  delete [] num_faces;


}





int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::cout << "Usage: " << argv[0] << " balOctsPrefix slack[0.0-1.0] num_procs" << std::endl;
    return 1;
  }

  //PetscInitialize(&argc, &argv, "options.hs", NULL);
  MPI_Init(&argc,&argv);
  int rank, size;
  char ptsFileName[256];
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  signal(SIGSEGV, handler);   // install our handler
  signal(SIGTERM, handler);   // install our handler

  double slack = atof(argv[2]);
  // @hari is this total number of procs or npes = size*q ?
  int q = atoi(argv[3]);

  // load in file ...
  std::vector<ot::TreeNode> balOct_H;
  std::vector<ot::TreeNode> balOct_M;


#ifdef HILBERT_ORDERING
  sprintf(ptsFileName, "%s_H_%d_%d.oct", argv[1], rank, size);
  ot::readNodesFromFile_binary(ptsFileName, balOct_H);

  //std::cout<<"FileName:"<<ptsFileName<<std::endl;
#else
  sprintf(ptsFileName, "%s_M_%d_%d.oct", argv[1], rank, size);
  ot::readNodesFromFile_binary(ptsFileName, balOct_M);
#endif
//  std::string filename = argv[1];
//  convert_r << filename << "_" << rank;
//  filename = convert_r.str();


  if (!rank) std::cout << "finished reading in balanced octree" << std::endl;

#ifdef HILBERT_ORDERING

  if (!rank) {
    std::cout << BLU << "===============================================" << NRM << std::endl;
    std::cout << RED " Flexible Parition for Hilbert  q:"<<q<<"\t slack:"<<slack<< NRM << std::endl;
    std::cout << BLU << "===============================================" << NRM << std::endl;
  }

  //assert(par::test::isUniqueAndSorted(balOct_H,MPI_COMM_WORLD));

  flexiblePartitionCalculation(balOct_H,slack,q);

#else
  if (!rank) {
    std::cout << BLU << "===============================================" << NRM << std::endl;
    std::cout << RED " Flexible Parition for Morton q:"<<q<<"\t slack:"<<slack<< NRM<< std::endl;
    std::cout << BLU << "===============================================" << NRM << std::endl;
  }

  flexiblePartitionCalculation(balOct_M,slack,q);
#endif


  MPI_Finalize();
  //PetscFinalize();
  return 0;
}


