//@author: Milinda Fernando.
// School of Computing
// University of Utah

// Contains the functions to calculate the mesh statistics.


#ifndef OCTREE_STATISTICS
#define OCTREE_STATISTICS

#include<vector>
#include <algorithm>
#include "mpi.h"
#include "petsc.h"
#include "sys.h"
#include "TreeNode.h"





//q: Number of pseudo processors
// mesh: List of octants
// stat[0] : min
// stat[1] : max
// stat[2] : mean
void calculateBoundaryFaces(const std::vector<ot::TreeNode> & mesh, int q,double* stat);
int calculateBoundaryFaces(const std::vector<ot::TreeNode>::const_iterator &beg, const std::vector<ot::TreeNode>::const_iterator &end);



#endif