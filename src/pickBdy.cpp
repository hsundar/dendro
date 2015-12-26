
/**
  @file pickBdy.C
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  */

#include "TreeNode.h"
#include "parUtils.h"

#ifdef __DEBUG__
#ifndef __DEBUG_OCT__
#define __DEBUG_OCT__
#endif
#endif

namespace ot {

  //Basic idea: Any octant with all possible neighbours inside blocks are
  //considered to be stable. This works on the assumption that the blocks are
  //internally balanced and inter-block boundaries are also balanced. So if all
  //possible neighbours of an octant are inside these blocks, then its entire
  //insulation layer is stable. Although, initially developed for checking
  //unstable inter-processor boundaries only,this idea is also applicable to
  //test stable remote octants in computing sendNodes for stage 2 communication
  //for balancing.
  //Assumptions: Blocks are sorted and unique and linear, nodes is sorted and unique and linear.
  //Every element in nodes is a decendant of or equal to some block.
  //Blocks are complete i.e. there are no gaps in the morton-space covered by
  //the blocks. 
  //Some Blocks may be empty.

#define PICK_IPBB_SET_TN_VALUE res.push_back(nodes[nodeCnt]);

#define PICK_IPBB_SET_UI_VALUE res.push_back(nodeCnt);

#define PICK_INTER_PROCESSOR_BOUNDARY_BLOCK(SetValueLine) {\
  res.clear();\
  unsigned int dim = firstBlock.getDim();\
  unsigned int maxDepth = firstBlock.getMaxDepth();\
  ot::TreeNode root(dim, maxDepth);\
  for(unsigned int nodeCnt = 0; nodeCnt < nodes.size(); nodeCnt++) {\
    unsigned int myMinX = nodes[nodeCnt].minX();\
    unsigned int myMinY = nodes[nodeCnt].minY();\
    unsigned int myMinZ = nodes[nodeCnt].minZ();\
    unsigned int myMaxX = nodes[nodeCnt].maxX();\
    unsigned int myMaxY = nodes[nodeCnt].maxY();\
    unsigned int myMaxZ = nodes[nodeCnt].maxZ();\
    unsigned int myLen = (myMaxX - myMinX);\
    unsigned int myLevel = nodes[nodeCnt].getLevel();\
    unsigned int negX = ((myMinX > 0) ? (myMinX - myLen) : myMinX);\
    unsigned int negY = ((myMinY > 0) ? (myMinY - myLen) : myMinY);\
    unsigned int negZ = ((myMinZ > 0) ? (myMinZ - myLen) : myMinZ);\
    unsigned int posX = ((myMaxX < (1u << maxDepth)) ? (myMaxX + myLen -1) : (myMaxX-1));\
    unsigned int posY = ((myMaxY < (1u << maxDepth)) ? (myMaxY + myLen -1) : (myMaxY-1));\
    unsigned int posZ = ((myMaxZ < (1u << maxDepth)) ? (myMaxZ + myLen -1) : (myMaxZ-1));\
    ot::TreeNode negCorner(negX, negY, negZ, myLevel, dim, maxDepth);\
    ot::TreeNode posCorner(posX, posY, posZ, maxDepth, dim, maxDepth);\
    bool add = true;\
    if( (negCorner >= nodes[0] && negCorner<=nodes[nodes.size()-1].getDLD() ) && ( posCorner >=nodes[0] && posCorner <= nodes[nodes.size()-1].getDLD()) ) {\
      add = false;\
    }\
    if (add) {\
      SetValueLine\
    }\
  }\
  return 1;\
}

//New Implementation

int pickInterProcessorBoundaryNodes(const std::vector<ot::TreeNode> & nodes,
    std::vector<unsigned int> & res, const ot::TreeNode & firstBlock,
    const ot::TreeNode & lastBlock) {

  PROF_PICK_BND_BEGIN

  PICK_INTER_PROCESSOR_BOUNDARY_BLOCK(PICK_IPBB_SET_UI_VALUE)

  PROF_PICK_BND_END
}

int pickInterProcessorBoundaryNodes(const std::vector<ot::TreeNode> & nodes,
    std::vector<ot::TreeNode> & res, const ot::TreeNode & firstBlock,
    const ot::TreeNode & lastBlock) {
  PROF_PICK_BND_BEGIN

  res.clear();


  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);


#ifdef HILBERT_ORDERING


  std::vector<ot::TreeNode> neighbourOct;
  ot::TreeNode tmp;
  for(int i=0;i<nodes.size();i++)
  {
    tmp=nodes[i].getTop();
    if(!tmp.isRoot())
      neighbourOct.push_back(tmp);

    tmp=nodes[i].getBottom();
    if(!tmp.isRoot())
      neighbourOct.push_back(tmp);

    tmp=nodes[i].getRight();
    if(!tmp.isRoot())
      neighbourOct.push_back(tmp);

    tmp=nodes[i].getLeft();
    if(!tmp.isRoot())
      neighbourOct.push_back(tmp);

    tmp=nodes[i].getFront();
    if(!tmp.isRoot())
      neighbourOct.push_back(tmp);

    tmp=nodes[i].getBack();
    if(!tmp.isRoot())
      neighbourOct.push_back(tmp);

    bool add=false;
    for(int j=0;j<neighbourOct.size();j++)
    {
      if(!(nodes[0]<= neighbourOct[j] && neighbourOct[j]<=nodes[nodes.size()-1]))
      {
        add=true;
        break;
      }
    }

    if(add)
    {
      res.push_back(nodes[i]);
    }

    neighbourOct.clear();

  }
#else
  unsigned int dim = firstBlock.getDim();
  unsigned int maxDepth = firstBlock.getMaxDepth();
  ot::TreeNode root(dim, maxDepth);
  for(unsigned int nodeCnt = 0; nodeCnt < nodes.size(); nodeCnt++) {
    unsigned int myMinX = nodes[nodeCnt].minX();
    unsigned int myMinY = nodes[nodeCnt].minY();
    unsigned int myMinZ = nodes[nodeCnt].minZ();
    unsigned int myMaxX = nodes[nodeCnt].maxX();
    unsigned int myMaxY = nodes[nodeCnt].maxY();
    unsigned int myMaxZ = nodes[nodeCnt].maxZ();
    unsigned int myLen = (myMaxX - myMinX);
    unsigned int myLevel = nodes[nodeCnt].getLevel();
    unsigned int negX = ((myMinX > 0) ? (myMinX - myLen) : myMinX);
    unsigned int negY = ((myMinY > 0) ? (myMinY - myLen) : myMinY);
    unsigned int negZ = ((myMinZ > 0) ? (myMinZ - myLen) : myMinZ);
    unsigned int posX = ((myMaxX < (1u << maxDepth)) ? (myMaxX + myLen -1) : (myMaxX-1));
    unsigned int posY = ((myMaxY < (1u << maxDepth)) ? (myMaxY + myLen -1) : (myMaxY-1));
    unsigned int posZ = ((myMaxZ < (1u << maxDepth)) ? (myMaxZ + myLen -1) : (myMaxZ-1));
    //@milinda @hari: This is not working for Hilbert ordering. Because we can't guarantee if negCorner and PosCorner is is the current partition current considering node also in the partition.
    ot::TreeNode negCorner(negX, negY, negZ, myLevel, dim, maxDepth);
    ot::TreeNode posCorner(posX, posY, posZ, maxDepth, dim, maxDepth);
    bool add = true;
    if( (negCorner >= firstBlock) && ( posCorner <= lastBlock.getDLD()) ) {\
      add = false;
    }
    if (add) {
      res.push_back(nodes[nodeCnt]);
    }
  }
#endif

  treeNodesTovtk(res,rank,"oda_res");

  return 1;
  
  PROF_PICK_BND_END
}


/*
#define PICK_INTER_PROCESSOR_BOUNDARY_BLOCK(SetValueLine) {\
TreeNode root(dim,maxDepth);\
res.clear();\
unsigned int nodeCnt = 0;\
unsigned int blkCnt = 0;\
while ( ( blkCnt < blocks.size() ) && (  nodeCnt < nodes.size() )  ) {\
if ( blocks[blkCnt] <= nodes[nodeCnt] ) {\
*//*Check if Ancestor*//* \        
if ( blocks[blkCnt].isAncestor(nodes[nodeCnt]) ) {\
  *//*Check if this node touches the internal boundary of the block.*//* \
    unsigned int myMinX = nodes[nodeCnt].minX();\
    unsigned int myMinY = nodes[nodeCnt].minY();\
    unsigned int myMinZ = nodes[nodeCnt].minZ();\
    unsigned int myMaxX = nodes[nodeCnt].maxX();\
    unsigned int myMaxY = nodes[nodeCnt].maxY();\
    unsigned int myMaxZ = nodes[nodeCnt].maxZ();\
    unsigned int bMinX = blocks[blkCnt].minX();\
    unsigned int bMinY = blocks[blkCnt].minY();\
    unsigned int bMinZ = blocks[blkCnt].minZ();\
    unsigned int bMaxX = blocks[blkCnt].maxX();\
    unsigned int bMaxY = blocks[blkCnt].maxY();\
    unsigned int bMaxZ = blocks[blkCnt].maxZ();\
    if ( (myMinX == bMinX) || (myMinY == bMinY) || (myMinZ == bMinZ) ||\
        (myMaxX == bMaxX) || (myMaxY == bMaxY) || (myMaxZ == bMaxZ) ) {\
      *//*node touches the boundary of this block from the inside*//* \
        std::vector<TreeNode> nh = nodes[nodeCnt].getAllNeighbours();\
        par::makeVectorUnique<ot::TreeNode>(nh, false);\
        unsigned int st = 0;\
        if (nh[0] == root) {\
          st = 1;\
        }\
      bool add = true;\
        *//*An indirect and cheap way to assert existence*//* \
        *//*of an ancestor is to*//* \
        *//*compare nh with the span of the domain owned by this processor
        *//* \
        *//*Note the second condition. You must compare the DLDs on both
        sides.*//*\
        *//* Just comparing reducedNh with last block or last block's DLD is
        not enough*//*\
        *//* This is because, the last block could be a decendant of
        reducedNh*//*\
        if(st < nh.size()) {\
          if( (nh[st] >= blocks[0]) &&\
              ((nh[nh.size()-1].getDLD())\
               <= (blocks[blocks.size()-1].getDLD())) ) {\
            add = false;\
          }\
        } else {\
          add = false;\
        }\
      if (add) {\
        SetValueLine\
      }\
    }\
  *//*There could be more nodes that are decendants of the same
    block.*//* \
    nodeCnt++;\
} else if ( blocks[blkCnt] == nodes[nodeCnt] ) {\
  std::vector<TreeNode> nh = nodes[nodeCnt].getAllNeighbours();\
    par::makeVectorUnique<TreeNode>(nh,false);\
    unsigned int st = 0;\
    if (nh[0] == root) {\
      st = 1;\
    }\
  bool add = true;\
    *//*An indirect and cheap way to assert existence*//* \
    *//*of an ancestor is to*//* \
    *//*compare nh with the span of the domain owned by this processor *//* \
    if(st < nh.size()) {\
      if( (nh[st] >= blocks[0]) &&\
          ((nh[nh.size()-1].getDLD())\
           <= (blocks[blocks.size()-1].getDLD())) ) {\
        add = false;\
      }\
    } else {\
      add = false;\
    }\
  if (add) {\
    SetValueLine\
  }\
  blkCnt++;\
} else {\
  *//*Neither Ancestor nor equal.*//* \
    blkCnt++;\
}\
  } else {\
    *//*Block can't be an Ancestor of the node.*//* \
      nodeCnt++;\
  }\
  }\
  return 1;\
  }
*/

/*
//Old Implementation
int pickInterProcessorBoundaryNodes(const std::vector<ot::TreeNode> & blocks, 
const std::vector<ot::TreeNode> &nodes, 
std::vector<unsigned int >& res, 
unsigned int dim, unsigned int maxDepth) {
PICK_INTER_PROCESSOR_BOUNDARY_BLOCK(PICK_IPBB_SET_UI_VALUE)
} // end pickboundary

int pickInterProcessorBoundaryNodes(const std::vector<ot::TreeNode> & blocks, 
const std::vector<ot::TreeNode> &nodes, 
std::vector<ot::TreeNode>& res, 
unsigned int dim, unsigned int maxDepth) {
PICK_INTER_PROCESSOR_BOUNDARY_BLOCK(PICK_IPBB_SET_TN_VALUE)
} // end pickboundary
*/

#undef PICK_IPBB_SET_TN_VALUE 
#undef PICK_IPBB_SET_UI_VALUE 
#undef PICK_INTER_PROCESSOR_BOUNDARY_BLOCK


// Oldest Implementation ...  
/*
   int pickInterProcessorBoundaryNodes(const std::vector<ot::TreeNode> & blocks, const std::vector<ot::TreeNode> &nodes, std::vector<ot::TreeNode>& res, unsigned int dim, unsigned int maxDepth){
   TreeNode root(dim,maxDepth);
   unsigned int resLen=0;
   res.resize(nodes.size());
   for(unsigned int i=0;i<nodes.size();i++) {
//A node is considered stable only if it is stable in all directions.
std::vector<TreeNode> nh = nodes[i].getAllNeighbours();
bool add=false;
for(unsigned int k=0;k<nh.size();k++) {
if(nh[k]==root){
//Stable from this direction 
continue;
}
bool foundAnc = false;
for(unsigned int j=0;j<blocks.size();j++) {
if((blocks[j] == nh[k]) || (blocks[j].isAncestor(nh[k]))){
foundAnc=true;
break;
}
}//end for j
if(!foundAnc) {
//Unstable from this direction 
add=true;
break;
}
}//end for k
if(add){
res[resLen] = nodes[i];
resLen++;
}
}//end for i
res.resize(resLen);
return 1;
}//end function
*/

}//end namespace



