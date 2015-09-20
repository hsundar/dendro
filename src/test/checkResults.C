/**
  @author Rahul Sampath, rahul.sampath@gmail.com
  */

#include "mpi.h"
#include "TreeNode.h"
#include "testUtils.h"
#include "externVars.h"
#include <cstdio>
#include <iostream>

int main(int argc, char **argv) {
  int res;
  unsigned int d;
  unsigned int mD;
  unsigned int numNodes;

  if (argc < 2) {
    std::cerr << "Usage: exe inpFile incCorn[1] maxLevDiff[1]" << std::endl;
    exit(1);
  }

  FILE *infile = fopen(argv[1], "r");
  bool incCorn = 1;
  unsigned int maxLevDiff = 1;
  if (argc > 2) {
    incCorn = (bool) atoi(argv[2]);
  }
  if (argc > 3) {
    maxLevDiff = atoi(argv[3]);
  }
  res = fscanf(infile, "%u", &d);
  res = fscanf(infile, "%u", &mD);
  res = fscanf(infile, "%u", &numNodes);
  std::vector<ot::TreeNode> nodes(numNodes);
  std::cout << " Starting to Read file: " << argv[1] << std::endl;
  for (unsigned int i = 0; i < numNodes; i++) {
    unsigned int myX, myY, myZ, myLev;
    res = fscanf(infile, "%u", &myX);
    res = fscanf(infile, "%u", &myY);
    res = fscanf(infile, "%u", &myZ);
    res = fscanf(infile, "%u", &myLev);
    ot::TreeNode tmpNode(myX, myY, myZ, myLev, d, mD);
    nodes[i] = tmpNode;
  }
  fclose(infile);
  std::cout << " Finished Reading " << nodes.size() << " octants." << std::endl;

  std::cout << " Is octree Sorted? : " << (seq::test::isSorted<ot::TreeNode>(nodes)) << std::endl;
  std::cout << " Is octree Unique & Sorted? : " <<
  (seq::test::isUniqueAndSorted<ot::TreeNode>(nodes)) << std::endl;
  std::cout << " Is octree Linear? : " << (ot::test::isLinear(nodes)) << std::endl;
  std::cout << " Is octree Complete? : " << (ot::test::isComplete(nodes)) << std::endl;
  std::cout << " Is octree Balanced? : " <<
  (ot::test::isBalanced(d, mD, "failedCheck", nodes, incCorn, maxLevDiff)) << std::endl;

}

