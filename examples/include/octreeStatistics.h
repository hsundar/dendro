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

//q: Number of pseudo processors
// mesh: List of octants
// slack: flexibility should be between [0,1]
#endif