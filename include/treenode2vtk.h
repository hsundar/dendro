#ifndef TREENODE22VTK
#define TREENODE22VTK

#include "TreeNode.h"
#include "Point.h"


#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#define VTK_HEXAHEDRON 12

// convert one treenode to vtk file compatibel string
//std::string treeNodeTovtk(const TreeNode& T,int mpi_rank,int hindex);

void treeNodesTovtk(std::vector<ot::TreeNode>& nodes,int mpi_rank,std::string vtk_file_name);




#endif
