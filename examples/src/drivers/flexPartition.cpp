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
int calculateBoundaryFaces(const std::vector<ot::TreeNode>::const_iterator &beg, const std::vector<ot::TreeNode>::const_iterator &end);

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
  // npes = size*q ?
  int q = atoi(argv[3]);

  // load in file ...
  std::vector<ot::TreeNode> balOct, tmpOct;

  std::ostringstream convert_r;
  std::string filename = argv[1];
  convert_r << filename << "_" << rank;
  filename = convert_r.str();
  ot::readNodesFromFile((char *)filename.c_str(), tmpOct);

  if (!rank) std::cout << "finished reading in balanced octree" << std::endl;

  // @milinda this might not be needed anymore
#ifdef HILBERT_ORDERING
  par::sampleSort(tmpOct, balOct, MPI_COMM_WORLD);
  if (!rank) std::cout << "finished sorting balanced octree" << std::endl;
#else
  std::swap(tmpOct, balOct);
#endif

  unsigned long localSz=balOct.size(), globalSz;
  MPI_Allreduce(&localSz, &globalSz, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

  int slackCnt = slack*globalSz/size/q;
  int part_size = localSz/q;

  if(!rank) std::cout << "slack size is " << slackCnt << " octants" << std::endl;

  assert(slackCnt < part_size);


  //----------------------------------------------------------------------
  //   FLEX
  //----------------------------------------------------------------------

  // 1. each process sends slackCnt octs to next/prev.
  std::vector<ot::TreeNode> ghosted(localSz+2*slackCnt);
  MPI_Status status;

  ghosted.insert(ghosted.begin()+slackCnt, balOct.begin(), balOct.end());


  int prev = (!rank)?rank-1:MPI_PROC_NULL;
  int next = (rank+1==size)?MPI_PROC_NULL:rank+1;

  ot::TreeNode* sendPrev = &(*(balOct.begin()));
  ot::TreeNode* sendNext = &(balOct[balOct.size() - slackCnt]);

  ot::TreeNode* recvPrev = &(*(ghosted.begin()));
  ot::TreeNode* recvNext = &(*(ghosted.begin()+localSz+slackCnt));

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
  int* partitions   = new int[q+1];
  int* num_faces    = new int[q+1];
  int faces;

  // partition 0
  partitions[0] = slackCnt;
  if (!rank) {
    // find first
    num_faces[0] = 6*localSz; // calculateBoundaryFaces(balOct.begin(), balOct.begin()+part_size);
    for (int j = 0; j < 2*slackCnt; ++j) {
      faces = calculateBoundaryFaces(ghosted.begin()+j, ghosted.begin()+slackCnt+part_size);
      if (faces < num_faces[0]) {
        num_faces[0] = faces; // we will update this later ...
        partitions[0] = j;
      }
    }
  } else {
    num_faces[0] = 6*localSz; // calculateBoundaryFaces(balOct.begin(), balOct.begin()+part_size);
    for (int j = 0; j < slackCnt; ++j) {
      faces = calculateBoundaryFaces(ghosted.begin()+slackCnt+j, ghosted.begin()+slackCnt+part_size);
      if (faces < num_faces[0]) {
        num_faces[0] = faces; // we will update this later ...
        partitions[0] = slackCnt+j;
      }
    }
  }
  std::vector<ot::TreeNode>::const_iterator first, last;

  // local partitions
  for (int i=1; i<q; ++i) {
    first = ghosted.begin() + partitions[i-1];
    num_faces[i] = 6*localSz;
    last = ghosted.begin() + i*part_size;
    for (int j = 0; j < 2*slackCnt; ++j) {
      faces = calculateBoundaryFaces(first, last + j);
      if (faces < num_faces[i]) {
        num_faces[i] = faces;
        partitions[i] = i*part_size + j;
      }
    }
  }
  // partition q
  if (rank == size-1) {
    first = ghosted.begin() + partitions[q-1];
    num_faces[q] = 6*localSz;
    last = ghosted.begin() + q*part_size;
    for (int j = 0; j < slackCnt; ++j) {
      faces = calculateBoundaryFaces(first, last + j);
      if (faces < num_faces[q]) {
        num_faces[q] = faces;
        partitions[q] = q*part_size + j;
      }
    }
  } else {
    first = ghosted.begin() + partitions[q-1];
    num_faces[q] = 6*localSz;
    last = ghosted.begin() + q*part_size;
    for (int j = 0; j < 2*slackCnt; ++j) {
      faces = calculateBoundaryFaces(first, last + j);
      if (faces < num_faces[q]) {
        num_faces[q] = faces;
        partitions[q] = q*part_size + j;
      }
    }
  }

  // @milinda now update num_faces
  for (int i=0; i<q; ++i) {
    num_faces[i] = calculateBoundaryFaces(ghosted.begin()+partitions[i], ghosted.begin()+partitions[i+1]);
    // num_elems = partitions[i+1] - partitions[i];
  }

  //----------------------------------------------------------------------
  if (!rank) {
    std::cout << GRN << "Finalizing ..." << NRM << std::endl;
  }

  delete [] partitions;
  delete [] num_faces;

  PetscFinalize();
  return 0;
}

//@author: Milinda Fernando.
// This function is to calculate the boundary faces
// Assume that the given octree vector is sorted.
int calculateBoundaryFaces(const std::vector<ot::TreeNode>::iterator first, const std::vector<ot::TreeNode>::iterator last) {

  ot::TreeNode R (first->getDim(), first->getMaxDepth());

  int mesh_nodes = last - first;
  assert(mesh_nodes > 1);

  int found_pt;
  bool isBdyElem;
  int num_boundary_faces=0;
  int num_boundary_elems=0;

  ot::TreeNode top;
  ot::TreeNode bottom;
  ot::TreeNode left;
  ot::TreeNode right;
  ot::TreeNode front;
  ot::TreeNode back;

  auto temp = first;
  auto last_local = last-1;

  unsigned long min_faces=1<<30;
  unsigned long max_faces=0;

  while(temp < last) {
    isBdyElem = false;
    top    = temp->getTop();
    bottom = temp->getBottom();
    left   = temp->getLeft();
    right  = temp->getRight();
    front  = temp->getFront();
    back   = temp->getBack();

    found_pt=(std::lower_bound(first, last, top, std::less<ot::TreeNode>()) - first);
    if(found_pt==0 || ( (found_pt==(last-1-first)) && (! last_local->isAncestor(top)))) {
      num_boundary_faces++;
      isBdyElem=true;
    }

    found_pt=(std::lower_bound(first, last, bottom, std::less<ot::TreeNode>()) - first);
    if(found_pt==0 || ((found_pt==(last-1-first)) && (!last_local->isAncestor(bottom)))) {
      num_boundary_faces++;
      isBdyElem=true;
    }

    found_pt=(std::lower_bound(first, last, left, std::less<ot::TreeNode>()) - first);
    if(found_pt==0 || ((found_pt==(last-1-first)) && (!last_local->isAncestor(left)))) {
      num_boundary_faces++;
      isBdyElem=true;
    }

    found_pt=(std::lower_bound(first, last, right, std::less<ot::TreeNode>()) - first);
    if(found_pt==0 || ((found_pt==(last-1-first)) && (!last_local->isAncestor(right)))) {
      num_boundary_faces++;
      isBdyElem=true;
    }
    found_pt=(std::lower_bound(first, last, front, std::less<ot::TreeNode>()) - first);
    if(found_pt==0 || ((found_pt==(last-1-first)) && (!last_local->isAncestor(front)))) {
      num_boundary_faces++;
      isBdyElem=true;
    }

    found_pt=(std::lower_bound(first, last, back, std::less<ot::TreeNode>()) - first);
    if(found_pt==0 || ((found_pt==(last-1-first)) && (!last_local->isAncestor(back)))) {
      num_boundary_faces++;
      isBdyElem=true;
    }

    if (isBdyElem) num_boundary_elems++;

    temp++;
  }

  return num_boundary_faces;
}
