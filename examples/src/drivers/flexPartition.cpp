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


int calculateBoundaryFaces(const std::vector<ot::TreeNode> & mesh, int q);

int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::cout << "Usage: " << argv[0] << " balOctsPrefix slack[0.0-1.0] num_procs" << std::endl;
    return 1;
  }

  PetscInitialize(&argc, &argv, "options.hs", NULL);

  int rank, size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double slack = atof(argv[2]);
  // @hari is this total number of procs or npes = size*q ?
  int q = atoi(argv[3]);

  // load in file ...
  std::vector<ot::TreeNode> balOct, tmpOct;

  std::ostringstream convert_r;
  std::string filename = argv[1];
  convert_r << filename << "_" << rank;
  filename = convert_r.str();
  ot::readNodesFromFile((char *)filename.c_str(), tmpOct);

  if (!rank) std::cout << "finished reading in balanced octree" << std::endl;

#ifdef HILBERT_ORDERING
  par::sampleSort(tmpOct, balOct, MPI_COMM_WORLD);
  if (!rank) std::cout << "finished sorting balanced octree" << std::endl;
#else
  std::swap(tmpOct, balOct);
#endif

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

  int prev = (!rank)?rank-1:MPI_PROC_NULL;
  int next = (rank+1==size)?MPI_PROC_NULL:rank+1;

  ot::TreeNode* sendPrev = &(*(balOct.begin()));
  ot::TreeNode* sendNext = &(balOct[balOct.size() - slackCnt]);

  ot::TreeNode* recvPrev = &(*(slack_prev.begin()));
  ot::TreeNode* recvNext = &(*(slack_next.begin()));

  // send to prev
  if (!rank)
    par::Mpi_Sendrecv<ot::TreeNode>(sendPrev, slackCnt, prev, 0,
                             recvPrev, slackCnt, prev, 0,
                             MPI_COMM_WORLD, &status);
  // send to next
  if (rank != (size-1))
    par::Mpi_Sendrecv<ot::TreeNode>(sendNext, slackCnt, next, 0,
                                    recvNext, slackCnt, next, 0,
                                    MPI_COMM_WORLD, &status);

  // Have the extra octants ...
  // 2. compute partitions.
  int* partitions = new int[q+1];
  // partition 0
  partitions[0] = 0;
  if (!rank) {
    slack_prev.insert(slack_prev.end(), balOct.begin(), balOct.begin() + slackCnt);
    for (int j = 0; j < slack; ++j) {

    }
  }
  // local partitions
  for (int i=1; i<q; ++i) {

  }
  // partition q
  partitions[q] = balOct.size();

  //----------------------------------------------------------------------
  if (!rank) {
    std::cout << GRN << "Finalizing ..." << NRM << std::endl;
  }

  PetscFinalize();
  return 0;
}

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
  unsigned long total_boundary_faces=0;

  unsigned long min_faces=1<<30;
  unsigned long max_faces=0;

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
    total_boundary_faces += num_boundary_faces;
    if (min_faces > num_boundary_faces) min_faces = num_boundary_faces;
    if (max_faces < num_boundary_faces) max_faces = num_boundary_faces;

    //std::cout<<"q:"<<j<<"\t "<<"number of boundary faces:"<<boundary_faces[j]<<std::endl;

  }

  unsigned long global_sum, global_max, global_min;

  MPI_Reduce(&total_boundary_faces, &global_sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&max_faces, &global_max, 1, MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&min_faces, &global_min, 1, MPI_UNSIGNED_LONG, MPI_MIN, 0, MPI_COMM_WORLD);

  return total_boundary_faces;
}
