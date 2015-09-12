
/**
  @file TreeNode.C
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  */

//#include "TreeNode.h"
#include "parUtils.h"
#include "seqUtils.h"
#include <vector>
#include <iterator>
#include "Point.h"


#include "TreeNode.h"


#ifdef __DEBUG__
#ifndef __DEBUG_TN__
#define __DEBUG_TN__
#endif
#endif

namespace ot {

std::vector<TreeNode> TreeNode::getAllNeighbours() const {
  /*
     0 = Left;  1 =  Right;  2 =  Front;  3 = Back;   4 = LeftBack;  5 = RightBack;  6 = LeftFront;  7 = RightFront;  8 = Top;
     9 = TopRight;  10 =  TopBack;  11 =  TopRightBack;  12 =  Bottom;  13 =  BottomBack;  14 =  TopLeft;  15 =  BottomLeft; 
     16 =  BottomRight;  17 =  TopFront;  18 =  BottomFront;  19 =  TopLeftFront;  20 =  TopRightFront;  21 =  BottomLeftFront;
     22 =  BottomRightFront;  23 =  TopLeftBack;  24 = BottomLeftBack;  25 = BottomRightBack;
     */
  std::vector<TreeNode> neighList;

  if (m_uiDim == 3) {
    neighList.resize(26);
    neighList[0] = getLeft();
    neighList[1] =  getRight();
    neighList[2] =  getFront();
    neighList[3] =  getBack();
    neighList[4] =  getLeftBack();
    neighList[5] = getRightBack();
    neighList[6] =  getLeftFront();
    neighList[7] =  getRightFront();
    neighList[8] =  getTop();
    neighList[9] = getTopRight();
    neighList[10] =  getTopBack();
    neighList[11] =  getTopRightBack();
    neighList[12] =  getBottom();
    neighList[13] =  getBottomBack();
    neighList[14] =  getTopLeft();
    neighList[15] =  getBottomLeft();
    neighList[16] =  getBottomRight();
    neighList[17] =  getTopFront();
    neighList[18] =  getBottomFront();
    neighList[19] =  getTopLeftFront();
    neighList[20] =  getTopRightFront();
    neighList[21] =   getBottomLeftFront();
    neighList[22] =   getBottomRightFront();
    neighList[23] =  getTopLeftBack();
    neighList[24] = getBottomLeftBack();
    neighList[25] = getBottomRightBack();
  } else if (m_uiDim == 2) {
    neighList.resize(8);
    neighList[0] = getLeft();
    neighList[1] =  getRight();
    neighList[2] =  getFront();
    neighList[3] =  getBack();
    neighList[4] =  getLeftBack();
    neighList[5] = getRightBack();
    neighList[6] =  getLeftFront();
    neighList[7] =  getRightFront();
  } else {
    neighList.resize(2);
    neighList[0] = getLeft();
    neighList[1] =  getRight();
  }
  return neighList;
}

std::vector<std::vector<TreeNode> > TreeNode::getAllB_Neighbours() const {
  /*
     0 = Left;  1 =  Right;  2 =  Front;  3 = Back;   4 = LeftBack;  5 = RightBack;  6 = LeftFront;  7 = RightFront;  8 = Top;
     9 = TopRight;  10 =  TopBack;  11 =  TopRightBack;  12 =  Bottom;  13 =  BottomBack;  14 =  TopLeft;  15 =  BottomLeft; 
     16 =  BottomRight;  17 =  TopFront;  18 =  BottomFront;  19 =  TopLeftFront;  20 =  TopRightFront;  21 =  BottomLeftFront;
     22 =  BottomRightFront;  23 =  TopLeftBack;  24 = BottomLeftBack;  25 = BottomRightBack;
     */

  std::vector<std::vector<TreeNode> > neighList;

  if (m_uiDim == 3) {
    neighList.resize(26);
    neighList[0] = getB_Left();
    neighList[1] =  getB_Right();
    neighList[2] =  getB_Front();
    neighList[3] =  getB_Back();
    neighList[4] =  getB_LeftBack();
    neighList[5] = getB_RightBack();
    neighList[6] =  getB_LeftFront();
    neighList[7] =  getB_RightFront();
    neighList[8] =  getB_Top();
    neighList[9] = getB_TopRight();
    neighList[10] =  getB_TopBack();
    neighList[11] =  getB_TopRightBack();
    neighList[12] =  getB_Bottom();
    neighList[13] =  getB_BottomBack();
    neighList[14] =  getB_TopLeft();
    neighList[15] =  getB_BottomLeft();
    neighList[16] =  getB_BottomRight();
    neighList[17] =  getB_TopFront();
    neighList[18] =  getB_BottomFront();
    neighList[19] =  getB_TopLeftFront();
    neighList[20] =  getB_TopRightFront();
    neighList[21] =   getB_BottomLeftFront();
    neighList[22] =   getB_BottomRightFront();
    neighList[23] =  getB_TopLeftBack();
    neighList[24] = getB_BottomLeftBack();
    neighList[25] = getB_BottomRightBack();
  } else if (m_uiDim == 2) {
    neighList.resize(8);
    neighList[0] = getB_Left();
    neighList[1] =  getB_Right();
    neighList[2] =  getB_Front();
    neighList[3] =  getB_Back();
    neighList[4] =  getB_LeftBack();
    neighList[5] = getB_RightBack();
    neighList[6] =  getB_LeftFront();
    neighList[7] =  getB_RightFront();
  } else {
    neighList.resize(2);
    neighList[0] = getB_Left();
    neighList[1] =  getB_Right();
  }
  return neighList;
}


int TreeNode::addBalancingDescendants(TreeNode other, std::vector<TreeNode>& seeds, bool incCorner) const {
#ifdef __DEBUG_TN__
  assert((other.getLevel()) > (getLevel() + 1));
  assert(areComparable(*this, other));
  assert(!(this->isAncestor(other)));
#endif
  TreeNode root(m_uiDim, m_uiMaxDepth);
  std::vector<TreeNode> nhs = other.getAllNeighbours();
  /*
     0 = Left;  1 =  Right;  2 =  Front;  3 = Back;   4 = LeftBack;  5 = RightBack;  6 = LeftFront;  7 = RightFront;  8 = Top;
     9 = TopRight;  10 =  TopBack;  11 =  TopRightBack;  12 =  Bottom;  13 =  BottomBack;  14 =  TopLeft;  15 =  BottomLeft; 
     16 =  BottomRight;  17 =  TopFront;  18 =  BottomFront;  19 =  TopLeftFront;  20 =  TopRightFront;  21 =  BottomLeftFront;
     22 =  BottomRightFront;  23 =  TopLeftBack;  24 = BottomLeftBack;  25 = BottomRightBack;
     */
  std::vector<unsigned int> dirs;
  if (m_uiDim == 1) {
    //ignore incCorner
    dirs.resize(2);
    for (unsigned int i = 0; i < dirs.size(); i++) {
      dirs[i] = i;
    }
  } else if (m_uiDim == 2) {
    if (incCorner) {
      dirs.resize(8);
      for (unsigned int i = 0; i < dirs.size(); i++) {
        dirs[i] = i;
      }
    } else {
      dirs.resize(4);
      for (unsigned int i = 0; i < dirs.size(); i++) {
        dirs[i] = i;
      }
    }
  } else {
    if (incCorner) {
      dirs.resize(26);
      for (unsigned int i = 0; i < dirs.size(); i++) {
        dirs[i] = i;
      }
    } else {
      dirs.resize(18);
      for (int i = 0; i < 11; i++) {
        dirs[i] = i;
      }
      //skip 11 = TRBk
      for (int i = 12; i < 19; i++) {
        dirs[i - 1] = i;
      }
    }
  }
  unsigned int oldSz = static_cast<unsigned int>(seeds.size());
  seeds.resize(seeds.size() + dirs.size());
  for (unsigned int i = 0; i < dirs.size(); i++) {
#ifdef __DEBUG_TN__
    assert(areComparable(*this, nhs[dirs[i]]));
#endif
    if (this->isAncestor(nhs[dirs[i]])) {
      seeds[oldSz + i] = nhs[dirs[i]].getParent();
    } else {
      seeds[oldSz + i] = root;
    } //end if - else
  } //end for i
  return 1;
} //end function


bool TreeNode::isBoundaryOctant(const TreeNode& block, int type, unsigned char *flags) const {
  unsigned char _flags = 0;

  unsigned int _x = block.getX();
  unsigned int _y = block.getY();
  unsigned int _z = block.getZ();	
  unsigned int _d = block.getLevel();

  /*
  // Block has to be an ancestor of the octant or equal to the octant.
  if( (!block.isAncestor(*this)) && (block != *this) ) {
  if (flags) {
   *flags = _flags;
   }
   return false;
   }
   */

  if ((type & NEGATIVE) == NEGATIVE) {
    // test if any of the anchor values matches those of the block ...
    if (m_uiX == _x) _flags |= X_NEG_BDY;
    if (m_uiY == _y) _flags |= Y_NEG_BDY;
    if (m_uiZ == _z) _flags |= Z_NEG_BDY;
  }

  if ((type & POSITIVE) == POSITIVE) {
    unsigned int len  = (unsigned int)(1u << (m_uiMaxDepth - getLevel()));
    unsigned int blen = ((unsigned int)(1u << (m_uiMaxDepth - _d))) - len;

    if (m_uiX == (_x + blen))  _flags |= X_POS_BDY;
    if (m_uiY == (_y + blen))  _flags |= Y_POS_BDY;
    if (m_uiZ == (_z + blen))  _flags |= Z_POS_BDY;
  }

  if (flags) {
    *flags = _flags;
  }
  if (_flags) {
    return true;
  }
  return false;
} //end function

bool TreeNode::isBoundaryOctant(int type, unsigned char *flags) const {
  unsigned char _flags = 0;
  if ((type & NEGATIVE) == NEGATIVE) {
    // test if any of the anchor values is zero ...  (sufficient ??? )
    if (!m_uiX) _flags |= X_NEG_BDY;
    if (!m_uiY) _flags |=  Y_NEG_BDY;
    if (!m_uiZ) _flags |=   Z_NEG_BDY;
  }

  if ((type & POSITIVE) == POSITIVE) {
    unsigned int len  = (unsigned int)(1u << (m_uiMaxDepth - getLevel()));
    unsigned int blen = ((unsigned int)(1u << m_uiMaxDepth)) - len;

    if (m_uiX == blen)  _flags |= X_POS_BDY;
    if (m_uiY == blen)  _flags |= Y_POS_BDY;
    if (m_uiZ == blen)  _flags |= Z_POS_BDY;
  }

  if (flags) *flags = _flags;
  if (_flags) return true;

  return false;
} //end function

void TreeNode::printTreeNode()
{
  
  std::cout<<"Coordinates: "<<this->m_uiX<<","<<this->m_uiY<<","<<this->m_uiZ<<std::endl;
  std::cout<<":Level: "<<this->m_uiLevel<<"\t Max Depth: "<<this->m_uiMaxDepth<<std::endl;

}

int TreeNode  ::addChildren(std::vector<ot::TreeNode>& children) const {
  unsigned int dim = m_uiDim;
  unsigned int maxDepth = m_uiMaxDepth;
  unsigned int childrenSz = children.size();
  children.resize(childrenSz + (1 << dim));
  
  //#define MORTON_ORDERING
  
  if ((m_uiLevel & ot::TreeNode::MAX_LEVEL) == maxDepth) {
    for (int i = 0; i < (1 << dim); i++) {
      children[childrenSz + i] = *this;
    }
    return 1;
  }
  //The check that lev < maxD is taken care of in the constructor.

  //Order: X first, Y next and Z last

  unsigned int len = (unsigned int)(1u << (maxDepth - ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1)));

  TreeNode   tmpNode0(1, m_uiX, m_uiY, m_uiZ, ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
  children[childrenSz + 0] = tmpNode0;

  TreeNode   tmpNode1(1, (m_uiX + len), m_uiY, m_uiZ, ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
  children[childrenSz + 1] = tmpNode1;

  if (dim >= 2) {
    TreeNode   tmpNode2(1, m_uiX, (m_uiY + len), m_uiZ, ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
    children[childrenSz + 2] = tmpNode2;

    TreeNode   tmpNode3(1, (m_uiX + len), (m_uiY + len), m_uiZ, ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
    children[childrenSz + 3] = tmpNode3;
  }

  if (dim == 3) {
    TreeNode   tmpNode4(1, m_uiX, m_uiY, (m_uiZ + len), ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
    children[childrenSz + 4] = tmpNode4;

    TreeNode   tmpNode5(1, (m_uiX + len), m_uiY, (m_uiZ + len), ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
    children[childrenSz + 5] = tmpNode5;

    TreeNode   tmpNode6(1, m_uiX, (m_uiY + len), (m_uiZ + len), ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
    children[childrenSz + 6] = tmpNode6;

    TreeNode   tmpNode7(1, (m_uiX + len), (m_uiY + len), (m_uiZ + len), ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
    children[childrenSz + 7] = tmpNode7;
  } //end if

#ifdef HILBERT_ORDERING
#pragma message("===FIX ME===")
  std::sort(children.begin(), children.end());
#endif
  return 1;
} //end function

int TreeNode  ::addBrothers(std::vector<TreeNode>& bros) const {
  unsigned int dim = m_uiDim;
  bros.resize(((1 << dim) - 1));
  if ((this->m_uiLevel & ot::TreeNode::MAX_LEVEL) == 0) {
    for (int i = 0; i < ((1 << dim) - 1); i++) {
      bros[i] = *this;
    } //end for
    return 1;
  } //end if

  TreeNode parent = this->getParent();
  std::vector<TreeNode> childrenOfParent;
  parent.addChildren(childrenOfParent);
  int k = 0;
  for (int i = 0; i < (1 << dim); i++) {
    if (childrenOfParent[i] != (*this)) {
      bros[k] = childrenOfParent[i];
      k++;
    }
  } //end for
  childrenOfParent.clear();
  return 1;
} //end fn.

std::vector<TreeNode> TreeNode::getSearchKeys(bool incCorners) {

  std::cout<<"Get Search Keys Called"<<std::endl;
  
  unsigned int myK = this->getChildNumber();
  //Morton Order: X,Y,Z
  bool zdir = (myK >= 4);
  bool ydir = ((myK - ((zdir) ? 4 : 0)) >= 2);
  bool xdir = (myK % 2);

  unsigned int xCor = this->minX();
  unsigned int yCor = this->minY();
  unsigned int zCor = this->minZ();
  if (xdir) {
    xCor = ((this->maxX()) - 1);
  }
  if (ydir) {
    yCor = ((this->maxY()) - 1);
  }
  if (zdir) {
    zCor = ((this->maxZ()) - 1);
  }

  TreeNode searchCorner(1, xCor, yCor, zCor, m_uiMaxDepth, m_uiDim, m_uiMaxDepth);
  std::vector<TreeNode> tList;
  if (xdir) {
    tList.push_back(searchCorner.getRight());
  } else {
    tList.push_back(searchCorner.getLeft());
  }
  if (m_uiDim > 1) {
    if (ydir) {
      tList.push_back(searchCorner.getBack());
    } else {
      tList.push_back(searchCorner.getFront());
    }
    if (incCorners || m_uiDim > 2) {
      if (ydir) {
        if (xdir) {
          tList.push_back(searchCorner.getRightBack());
        } else {
          tList.push_back(searchCorner.getLeftBack());
        }
      } else {
        if (xdir) {
          tList.push_back(searchCorner.getRightFront());
        } else {
          tList.push_back(searchCorner.getLeftFront());
        }
      }
    }
  }
  if (m_uiDim > 2) {
    if (zdir) {
      tList.push_back(searchCorner.getTop());
      if (xdir) {
        tList.push_back(searchCorner.getTopRight());
      } else {
        tList.push_back(searchCorner.getTopLeft());
      }
      if (ydir) {
        tList.push_back(searchCorner.getTopBack());
      } else {
        tList.push_back(searchCorner.getTopFront());
      }
    } else {
      tList.push_back(searchCorner.getBottom());
      if (xdir) {
        tList.push_back(searchCorner.getBottomRight());
      } else {
        tList.push_back(searchCorner.getBottomLeft());
      }
      if (ydir) {
        tList.push_back(searchCorner.getBottomBack());
      } else {
        tList.push_back(searchCorner.getBottomFront());
      }
    }
    if (incCorners) {
      if (zdir) {
        if (xdir) {
          if (ydir) {
            tList.push_back(searchCorner.getTopRightBack());
          } else {
            tList.push_back(searchCorner.getTopRightFront());
          }
        } else {
          if (ydir) {
            tList.push_back(searchCorner.getTopLeftBack());
          } else {
            tList.push_back(searchCorner.getTopLeftFront());
          }
        }
      } else {
        if (xdir) {
          if (ydir) {
            tList.push_back(searchCorner.getBottomRightBack());
          } else {
            tList.push_back(searchCorner.getBottomRightFront());
          }
        } else {
          if (ydir) {
            tList.push_back(searchCorner.getBottomLeftBack());
          } else {
            tList.push_back(searchCorner.getBottomLeftFront());
          }
        }
      }
    }
  }

  return tList;
} //end function

std::vector<bool> ot::TreeNode  ::getMorton() const {
  unsigned int dim = getDim();
  unsigned int maxDepth = getMaxDepth();
  unsigned int Ln = 1;
  if (maxDepth > 0) {
    Ln = binOp::binLength(maxDepth);
  }
  unsigned int const N = (dim * maxDepth) + Ln;
  std::vector<bool> numBin(N);
  std::vector<bool> xBin(maxDepth);
  std::vector<bool> yBin(maxDepth);
  std::vector<bool> zBin(maxDepth);
  std::vector<bool> levBin(Ln);

  //create Binary Representations
  binOp::toBin(getX(), maxDepth, xBin);
  if (dim > 1) {binOp::toBin(getY(), maxDepth, yBin); }
  if (dim > 2) {binOp::toBin(getZ(), maxDepth, zBin); }
  binOp::toBin(getLevel(), Ln, levBin);

  //Interleave bits
  if (dim > 2) {
    for (unsigned int j = 0; j < maxDepth; j++) {
      numBin[(j * dim)] = zBin[j];
      numBin[((j * dim) + 1)] = yBin[j];
      numBin[((j * dim) + 2)] = xBin[j];
    } //end for
  } else if (dim > 1) {
    for (unsigned int j = 0; j < maxDepth; j++) {
      numBin[(j * dim)] = yBin[j];
      numBin[((j * dim) + 1)] = xBin[j];
    } //end for
  } else {
    for (unsigned int j = 0; j < maxDepth; j++) {
      numBin[j] = xBin[j];
    } //end for
  } //end if-else

  //Append level
  for (unsigned int j = 0; j < Ln; j++) {
    numBin[((maxDepth * dim) + j)] = levBin[j];
  }

  return numBin;
}

//Constructors...
TreeNode :: TreeNode() {
  m_uiX = m_uiY = m_uiZ = m_uiLevel =
     m_uiWeight = m_uiDim = m_uiMaxDepth = 0;
  
  calculateTreeNodeRotation();  
     
}

TreeNode  :: TreeNode(const int dummy, const unsigned int x, const unsigned int y,
                      const unsigned int z, const unsigned int lev, const unsigned int dim, const unsigned int maxDepth) {
  m_uiDim = dim;
  m_uiMaxDepth = maxDepth;
  m_uiX = x;
  if (dim > 1) {m_uiY = y; } else {m_uiY = 0; }
  if (dim > 2) {m_uiZ = z; } else {m_uiZ = 0; }

  m_uiLevel = lev;
  m_uiWeight = 1;
  
  calculateTreeNodeRotation();
  
} //end function

TreeNode  :: TreeNode(const unsigned int dim, const unsigned int maxDepth) {
  m_uiX = 0;
  m_uiY = 0;
  m_uiZ = 0;
  m_uiLevel = 0;
  m_uiWeight = 1;
  m_uiDim = dim;
  m_uiMaxDepth = maxDepth;
#ifdef __DEBUG_TN__
  if ((dim != 1) && (dim != 2) && (dim != 3)) {
    std::cout << "Wrong Value for dim: " << dim << std::endl;
  }
#endif
  assert((dim == 1) || (dim == 2) || (dim == 3));
  
  calculateTreeNodeRotation();
} //end function

TreeNode  :: TreeNode(const unsigned int x, const unsigned int y,
                      const unsigned int z, const unsigned int lev, const unsigned int dim, const unsigned int maxDepth) {
  m_uiDim = dim;
  m_uiMaxDepth = maxDepth;
  m_uiX = x;
  if (dim > 1) {m_uiY = y; } else {m_uiY = 0; }
  if (dim > 2) {m_uiZ = z; } else {m_uiZ = 0; }

  m_uiLevel = lev;
  m_uiWeight = 1;

#ifdef __DEBUG_TN__
  if ((dim != 1) && (dim != 2) && (dim != 3)) {
    std::cout << "Wrong Value for dim: " << dim << std::endl;
  }
  assert(m_uiX < ((unsigned int)(1u << maxDepth)));
  assert((m_uiX % ((unsigned int)(1u << (maxDepth - lev)))) == 0);
  assert(m_uiY < ((unsigned int)(1u << maxDepth)));
  assert((m_uiY % ((unsigned int)(1u << (maxDepth - lev)))) == 0);
  assert(m_uiZ < ((unsigned int)(1u << maxDepth)));
  assert((m_uiZ % ((unsigned int)(1u << (maxDepth - lev)))) == 0);
  assert((dim == 1) || (dim == 2) || (dim == 3));
#endif
    
  calculateTreeNodeRotation();

} //end function

//copy constructor
TreeNode  :: TreeNode(const TreeNode& other) {
  
  m_uiX = other.m_uiX;
  m_uiY = other.m_uiY;
  m_uiZ = other.m_uiZ;
  m_uiLevel = other.m_uiLevel;
  m_uiWeight = other.m_uiWeight;
  m_uiDim = other.m_uiDim;
  m_uiMaxDepth = other.m_uiMaxDepth;
  
//   if(m_uiDim==2)
//   {
//     rotation=new int[4];
//     rot_index=new int[4];
//     for(int i=0;i<4;i++){
//       rotation[i]=other.rotation[i];
//       rot_index[i]=other.rot_index[i];
//     }
//   }else if(m_uiDim==3)
//   {
//     rotation=new int[8];
//     rot_index=new int[8];
//     
//     for(int i=0;i<8;i++){
//       rotation[i]=other.rotation[i];
//       rot_index[i]=other.rot_index[i];
//     }
//     
//   }
  
  calculateTreeNodeRotation();
  
} //end function

/*
 *Operator overiding implementations. 
 *Comparison operators are based on Hilbert ordering of the points. 
 */


std::ostream& operator<<(std::ostream& os, TreeNode const& other) {
  return (os << other.getX() << " " << other.getY() << " " << other.getZ() << " " << other.getLevel());
} //end fn.

/*
 * Date: 30/05/2015
 * ===========Hilbert Ordering Implementations==========
 * ====================START============================
 */

void TreeNode::calculateTreeNodeRotation()
{
  
  unsigned int xl = 0;
  unsigned int yl = 0;
  unsigned int zl = 0;

  unsigned int len = 1 << this->m_uiMaxDepth;
  int count = 0;
  unsigned int index1 = 0;
  unsigned int index2 = 0;
  
  unsigned int ncaX,ncaY,ncaZ,ncaLev; // considering the current node as the NCA. 
  ncaX=this->m_uiX;
  ncaY=this->m_uiY;
  ncaZ=this->m_uiZ;
  ncaLev=this->m_uiLevel;
  int index_temp;

  if (m_uiDim == 2) {
     //rotation=new int[4]{0,1,2,3};
     //rot_index=new int[4]{0,1,2,3};
     int current_rot=0;
          
    while ((xl != ncaX || yl != ncaY || zl != ncaZ || count != ncaLev)) {
      
      len >>=1;
      index1 = 0;
      if (ncaX >= (len + xl)) {
        index1 += 1;
        xl += len;
        if (ncaY < (len + yl)) index1 += 2;
      }
      if (ncaY >= (len + yl)) {index1 += 1; yl += len; }

      //rotate(index1, rotation, rot_index, 2);
      index_temp=rotations_2d[current_rot].rot_index[index1];
      current_rot=HILBERT_2D_TABLE[current_rot][index_temp];
      
      count++;

    }

    rotation_2d=rotations_2d[current_rot];
    

  } else if (m_uiDim == 3) {
//     rotation = new int[8]{ 0, 1, 2, 3, 4, 5, 6, 7 }; // Initial rotation
//     rot_index = new int[8]{ 0, 1, 2, 3, 4, 5, 6, 7 }; // Initial rotation indices
    
    int current_rot=0;
    
    while ((xl != ncaX || yl != ncaY || zl != ncaZ || count != ncaLev)/*&& len >0*/) {

      len >>= 1;

      index1 = 0;
      if (ncaZ < (len + zl)) {
        if (ncaX >= (len + xl)) {
          index1 += 1;
          xl += len;
          if (ncaY < (len + yl)) index1 += 2;
        }
        if (ncaY >= (len + yl)) {
          index1 += 1;
          yl += len;
        }
      } else {
        index1 = 4;
        zl += len;
        if (ncaX < (len + xl)) {
          index1 += 1;
          if (ncaY < (len + yl)) index1 += 2;
        } else {
          xl += len;
        }
        if (ncaY >= (len + yl)) {
          index1 += 1;
          yl += len;
        }
      }

      index_temp=rotations_3d[current_rot].rot_index[index1];
      current_rot=HILBERT_3D_TABLE[current_rot][index_temp];
      count++;

    }
    rotation_3d=rotations_3d[current_rot];
  
}
//std::cout<<"Tree Node Rotation Calculation Completed"<<std::endl;

}


/*
 * ===========Hilbert Ordering Implementations==========
 * ====================END============================
 */

//Assignment operator
//No checks for dim or maxD. It's ok to change dim and maxD using the
//assignment operator.
TreeNode& TreeNode  :: operator = (TreeNode   const& other) {
  if (this == (&other)) {return *this;}
  m_uiX = other.m_uiX;
  m_uiY = other.m_uiY;
  m_uiZ = other.m_uiZ;
  m_uiLevel = other.m_uiLevel;
  m_uiWeight = other.m_uiWeight;
  m_uiDim = other.m_uiDim;
  m_uiMaxDepth = other.m_uiMaxDepth;
  calculateTreeNodeRotation();
  
//   if(m_uiDim==2)
//   {
//     rotation=new int[4];
//     rot_index=new int[4];
//     for(int i=0;i<4;i++){
//       rotation[i]=other.rotation[i];
//       rot_index[i]=other.rot_index[i];
//     }
//   }else if(m_uiDim==3)
//   {
//     rotation=new int[8];
//     rot_index=new int[8];
//     
//     for(int i=0;i<8;i++){
//       rotation[i]=other.rotation[i];
//       rot_index[i]=other.rot_index[i];
//     }
//     
//   }
  
  return *this;
} //end fn.

int TreeNode::pickInternalBoundaryCells(std::vector<TreeNode>& allInternal,
                                        std::vector<TreeNode>& boundary)  const {
  unsigned int dim = this->getDim();

  //Pre-alloc
  boundary.resize(allInternal.size());
  unsigned int boundarySz = 0;

  int myMinX = this->minX();
  int myMinY = this->minY();
  int myMinZ = this->minZ();
  int myMaxX = this->maxX();
  int myMaxY = this->maxY();
  int myMaxZ = this->maxZ();

  for (unsigned int i = 0; i < allInternal.size(); i++) {
#ifdef __DEBUG_TN__
    assert(areComparable(*this, allInternal[i]));
    if (!(this->isAncestor(allInternal[i]))) {
      std::cout << (*this) << RED << " is not an Ancestor of " << NRM << allInternal[i] << std::endl;
    }
    assert(this->isAncestor(allInternal[i]));
#endif
    int testMinX = allInternal[i].minX();
    int testMinY = allInternal[i].minY();
    int testMinZ = allInternal[i].minZ();
    int testMaxX = allInternal[i].maxX();
    int testMaxY = allInternal[i].maxY();
    int testMaxZ = allInternal[i].maxZ();
    if ((myMinX == testMinX) || (myMaxX == testMaxX)) {
      boundary[boundarySz] = (allInternal[i]);
      boundarySz++;
    } else if (dim > 1) {
      if ((myMinY == testMinY) || (myMaxY == testMaxY)) {
        boundary[boundarySz] = (allInternal[i]);
        boundarySz++;
      } else if (dim > 2) {
        if ((myMinZ == testMinZ) || (myMaxZ == testMaxZ)) {
          boundary[boundarySz] = (allInternal[i]);
          boundarySz++;
        } //end if
      } //end if-else-if
    } //end if-else-if on boundary
  } //end for

  boundary.resize(boundarySz);
  return 1;
} //end function


int TreeNode::splitInternalAndBoundaryCells(std::vector<TreeNode>& allInternal,
                                            std::vector<TreeNode>& onlyInternal, std::vector<TreeNode>& boundary) const {
  unsigned int dim = this->getDim();

  //Pre-alloc
  boundary.resize(allInternal.size());
  onlyInternal.resize(allInternal.size());
  unsigned int boundarySz = 0;
  unsigned int onlyInternalSz = 0;

  int myMinX = this->minX();
  int myMinY = this->minY();
  int myMinZ = this->minZ();
  int myMaxX = this->maxX();
  int myMaxY = this->maxY();
  int myMaxZ = this->maxZ();

  for (unsigned int i = 0; i < allInternal.size(); i++) {
#ifdef __DEBUG_TN__
    assert(areComparable(*this, allInternal[i]));
    if (!(this->isAncestor(allInternal[i]))) {
      std::cout << (*this) << RED << " is not an Ancestor of " << NRM << allInternal[i] << std::endl;
    }
    assert(this->isAncestor(allInternal[i]));
#endif
    int testMinX = allInternal[i].minX();
    int testMinY = allInternal[i].minY();
    int testMinZ = allInternal[i].minZ();
    int testMaxX = allInternal[i].maxX();
    int testMaxY = allInternal[i].maxY();
    int testMaxZ = allInternal[i].maxZ();
    if ((myMinX == testMinX) || (myMaxX == testMaxX)) {
      boundary[boundarySz] = (allInternal[i]);
      boundarySz++;
    } else if (dim > 1) {
      if ((myMinY == testMinY) || (myMaxY == testMaxY)) {
        boundary[boundarySz] = (allInternal[i]);
        boundarySz++;
      } else if (dim > 2) {
        if ((myMinZ == testMinZ) || (myMaxZ == testMaxZ)) {
          boundary[boundarySz] = (allInternal[i]);
          boundarySz++;
        } else {
          onlyInternal[onlyInternalSz] = (allInternal[i]);
          onlyInternalSz++;
        } //end if
      } else {
        onlyInternal[onlyInternalSz] = (allInternal[i]);
        onlyInternalSz++;
      } //end if-else-if
    } else {
      onlyInternal[onlyInternalSz] = (allInternal[i]);
      onlyInternalSz++;
    } //end if-else-if on boundary
  } //end for

  boundary.resize(boundarySz);
  onlyInternal.resize(onlyInternalSz);
  return 1;
} //end function

int TreeNode ::balanceSubtree(std::vector<TreeNode>& inp,
                              std::vector<TreeNode>& out, bool incCorner,
                              bool isSortedCompleteLinear) const {

  PROF_BAL_SUBTREE_BEGIN

     unsigned int dim = getDim();
  unsigned int maxDepth = getMaxDepth();
  unsigned int maxCornerX = this->maxX();
  unsigned int maxCornerY = this->maxY();
  unsigned int maxCornerZ = this->maxZ();
  unsigned int minLevInp = maxDepth;
  unsigned int maxLevInp = 0;

  out.clear();
  if ((maxDepth == 0) || (getLevel() == maxDepth) || (inp.empty())) {
    //No decendants.
    PROF_BAL_SUBTREE_END
  }

  std::vector<TreeNode> * workVec = NULL;
  unsigned int *workVecSz = NULL;

  if (maxDepth - getLevel()) {
    workVec = new std::vector<TreeNode>[maxDepth - getLevel()];
    workVecSz = new unsigned int[maxDepth - getLevel()];
  }

  for (unsigned int i = 0; i < (maxDepth - getLevel()); i++) {
    workVecSz[i] = 0;
  }

  for (unsigned int i = 0; i < inp.size(); i++) {
#ifdef __DEBUG_TN__
    assert(areComparable(*this, inp[i]));
#endif
    if (this->isAncestor(inp[i])) {
      unsigned int currOctLev = inp[i].getLevel();
      workVecSz[maxDepth - currOctLev]++;
      if (currOctLev < minLevInp) {
        minLevInp = currOctLev;
      }
      if (currOctLev > maxLevInp) {
        maxLevInp = currOctLev;
      }
    } //end if decendant of block
  } //end for

  if (isSortedCompleteLinear && ((maxLevInp - minLevInp) < 2)) {
    //Already balanced
    out = inp;

    if (workVecSz) {
      delete[] workVecSz;
      workVecSz = NULL;
    }

    if (workVec) {
      delete[] workVec;
      workVec = NULL;
    }

    PROF_BAL_SUBTREE_END
  }

  for (unsigned int i = 0; i < (maxDepth - getLevel()); i++) {
    workVec[i].resize(workVecSz[i]);
    workVecSz[i] = 0;
  }

  for (unsigned int i = 0; i < inp.size(); i++) {
#ifdef __DEBUG_TN__
    assert(areComparable(*this, inp[i]));
#endif
    if (this->isAncestor(inp[i])) {
      workVec[maxDepth - (inp[i].getLevel())][workVecSz[maxDepth - (inp[i].getLevel())]] = inp[i];
      workVecSz[maxDepth - (inp[i].getLevel())]++;
    }
  } //end for

  if (workVecSz) {
    delete[] workVecSz;
    workVecSz = NULL;
  }

  TreeNode tmpRoot(getDim(), getMaxDepth());
  for (unsigned int i = 0; i < maxDepth - getLevel() - 1; i++) {
    //This also sorts.
    seq::makeVectorUnique(workVec[i], false);
    //Remove Brothers from a sorted list.
    std::vector<TreeNode> tmpList(workVec[i].size());
    unsigned int tmpListSz = 0;
    if (!workVec[i].empty()) {
      tmpList[tmpListSz] = workVec[i][0];
      tmpListSz++;
    }
    for (unsigned int j = 1; j < workVec[i].size(); j++) {
      if (workVec[i][j - 1].getParent() != workVec[i][j].getParent()) {
        tmpList[tmpListSz] = workVec[i][j];
        tmpListSz++;
      }
    } //end for j
    workVec[i].clear();
    tmpList.resize(tmpListSz);
    workVec[i] = tmpList;
    tmpList.clear();
    unsigned int oldOutSz = out.size();
    unsigned int newOut = 0;
    out.resize(oldOutSz + (workVec[i].size() * (1 << dim)));
    unsigned int oldNextWsz = workVec[i + 1].size();
    unsigned int newNextW = 0;
    unsigned int wLen = workVec[i].size();
    if (wLen > 0) {
      workVec[i + 1].resize(oldNextWsz + (26 * (workVec[i].size())));
    }
    for (unsigned int j = 0; j < wLen; j++) {
      out[oldOutSz + newOut] = workVec[i][j];
      newOut++;
      std::vector<TreeNode> bros;
      workVec[i][j].addBrothers(bros);
      for (int k = 0; k < ((1 << dim) - 1); k++) {
        out[oldOutSz + newOut] = bros[k];
        newOut++;
      } //end for k
      bros.clear();

      TreeNode  parent = workVec[i][j].getParent();
      TreeNode  tmpNode(getDim(), getMaxDepth());
      //Add all the neighbours of parent to workVec[i+1]

      tmpNode = parent.getLeft();
      if (tmpNode > tmpRoot) {
        if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
            (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
            (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
          workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
          newNextW++;
        }
      }

      tmpNode = parent.getRight();
      if (tmpNode > tmpRoot) {
        if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
            (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
            (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
          workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
          newNextW++;
        }
      }

      tmpNode = parent.getTop();
      if (tmpNode > tmpRoot) {
        if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
            (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
            (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
          workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
          newNextW++;
        }
      }

      tmpNode = parent.getBottom();
      if (tmpNode > tmpRoot) {
        if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
            (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
            (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
          workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
          newNextW++;
        }
      }

      tmpNode = parent.getFront();
      if (tmpNode > tmpRoot) {
        if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
            (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
            (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
          workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
          newNextW++;
        }
      }

      tmpNode = parent.getBack();
      if (tmpNode > tmpRoot) {
        if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
            (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
            (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
          workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
          newNextW++;
        }
      }
      if (dim == 3 || incCorner) {
        tmpNode = parent.getTopLeft();
        if (tmpNode > tmpRoot) {
          if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
              (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
              (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
            workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
            newNextW++;
          }
        }

        tmpNode = parent.getTopRight();
        if (tmpNode > tmpRoot) {
          if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
              (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
              (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
            workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
            newNextW++;
          }
        }

        tmpNode = parent.getBottomLeft();
        if (tmpNode > tmpRoot) {
          if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
              (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
              (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
            workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
            newNextW++;
          }
        }

        tmpNode = parent.getBottomRight();
        if (tmpNode > tmpRoot) {
          if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
              (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
              (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
            workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
            newNextW++;
          }
        }

        tmpNode = parent.getLeftFront();
        if (tmpNode > tmpRoot) {
          if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
              (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
              (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
            workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
            newNextW++;
          }
        }

        tmpNode = parent.getRightFront();
        if (tmpNode > tmpRoot) {
          if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
              (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
              (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
            workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
            newNextW++;
          }
        }

        tmpNode = parent.getTopFront();
        if (tmpNode > tmpRoot) {
          if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
              (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
              (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
            workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
            newNextW++;
          }
        }

        tmpNode = parent.getBottomFront();
        if (tmpNode > tmpRoot) {
          if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
              (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
              (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
            workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
            newNextW++;
          }
        }
        tmpNode = parent.getLeftBack();
        if (tmpNode > tmpRoot) {
          if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
              (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
              (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
            workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
            newNextW++;
          }
        }

        tmpNode = parent.getRightBack();
        if (tmpNode > tmpRoot) {
          if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
              (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
              (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
            workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
            newNextW++;
          }
        }

        tmpNode = parent.getTopBack();
        if (tmpNode > tmpRoot) {
          if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
              (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
              (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
            workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
            newNextW++;
          }
        }

        tmpNode = parent.getBottomBack();
        if (tmpNode > tmpRoot) {
          if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
              (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
              (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
            workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
            newNextW++;
          }
        }


      } //end if-incCorner or dim=3
      if (incCorner) {
        tmpNode = parent.getTopLeftFront();
        if (tmpNode > tmpRoot) {
          if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
              (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
              (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
            workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
            newNextW++;
          }
        }

        tmpNode = parent.getTopRightFront();
        if (tmpNode > tmpRoot) {
          if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
              (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
              (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
            workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
            newNextW++;
          }
        }

        tmpNode = parent.getBottomLeftFront();
        if (tmpNode > tmpRoot) {
          if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
              (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
              (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
            workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
            newNextW++;
          }
        }

        tmpNode = parent.getBottomRightFront();
        if (tmpNode > tmpRoot) {
          if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
              (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
              (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
            workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
            newNextW++;
          }
        }

        tmpNode = parent.getTopLeftBack();
        if (tmpNode > tmpRoot) {
          if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
              (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
              (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
            workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
            newNextW++;
          }
        }

        tmpNode = parent.getTopRightBack();
        if (tmpNode > tmpRoot) {
          if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
              (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
              (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
            workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
            newNextW++;
          }
        }

        tmpNode = parent.getBottomLeftBack();
        if (tmpNode > tmpRoot) {
          if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
              (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
              (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {
            workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
            newNextW++;
          }
        }

        tmpNode = parent.getBottomRightBack();
        if (tmpNode > tmpRoot) {
          if ((tmpNode.m_uiX >= this->m_uiX) && (tmpNode.m_uiY >= this->m_uiY) &&
              (tmpNode.m_uiZ >= this->m_uiZ) && (tmpNode.m_uiX < maxCornerX) &&
              (tmpNode.m_uiY < maxCornerY) && (tmpNode.m_uiZ < maxCornerZ)) {

            workVec[i + 1][oldNextWsz + newNextW] = tmpNode;
            newNextW++;
          }
        }
      } //end if-incCorner
    } //end for j
    workVec[i].clear();
    if (wLen > 0) {
      workVec[i + 1].resize(oldNextWsz + newNextW);
    }
  } //end for i

  if (workVec) {
    if (!(workVec[maxDepth - 1 - getLevel()].empty())) {
      unsigned int oldOutSz = out.size();
      unsigned int newOut = 0;
      out.resize(oldOutSz + (1 << dim));
      out[oldOutSz + newOut] = workVec[maxDepth - 1 - getLevel()][0];
      newOut++;
      std::vector<TreeNode> bros;
      workVec[maxDepth - 1 - getLevel()][0].addBrothers(bros);
      workVec[maxDepth - 1 - getLevel()].clear();
      for (int k = 0; k < ((1 << dim) - 1); k++) {
        out[oldOutSz + newOut] = bros[k];
        newOut++;
      } //end for k

      bros.clear();
    } //end if
    delete[] workVec;
    workVec = NULL;
  }

  //Sort and make unique.
  seq::makeVectorUnique<TreeNode>(out, false);
  //Remove overlaps.
  std::vector<unsigned int> toRem(out.size());
  unsigned int remCtr = 0;
  //compare i+1 with i
  for (unsigned int i = 1; i < out.size(); i++) {
#ifdef __DEBUG_TN__
    assert(areComparable(out[i - 1], out[i]));
#endif
    if (out[i - 1].isAncestor(out[i])) {
      toRem[remCtr] = (i - 1);
      remCtr++;
    } //end if
  } //end for i
  toRem.resize(remCtr);

  if (!toRem.empty()) {
    //Note toRem is sorted.
    std::vector<TreeNode> tmpRem(out.size() - toRem.size());
    unsigned int pos = 0; unsigned int kk = 0;
    for (unsigned int i = 0; i < toRem.size(); i++) {
      for (unsigned int j = pos; j < toRem[i]; j++) {
        tmpRem[kk] = out[j];
        kk++;
      } //end for j
      pos = toRem[i] + 1;
    } //end for i
    for (unsigned int j = pos; j < out.size(); j++) {
      tmpRem[kk] = out[j];
      kk++;
    } //end for j
    out = tmpRem;
    tmpRem.clear();
  } //end if
  toRem.clear();

  PROF_BAL_SUBTREE_END
} //end function

int TreeNode ::completeSubtree(std::vector<TreeNode>& inp,
                               std::vector<TreeNode>& out) const {

  PROF_COMPLETE_SUBTREE_BEGIN

     unsigned int dim = m_uiDim;
  unsigned int maxDepth = m_uiMaxDepth;

  out.clear();
  if ((maxDepth == 0) || (getLevel() == maxDepth) || (inp.empty())) {
    //No decendants.
    PROF_COMPLETE_SUBTREE_END
  }

  std::vector<TreeNode> * workVec = NULL;
  unsigned int *workVecSz = NULL;

  if (maxDepth - getLevel()) {
    workVec = new std::vector<TreeNode>[maxDepth - getLevel()];
    workVecSz = new unsigned int[maxDepth - getLevel()];
  }

  for (unsigned int i = 0; i < (maxDepth - getLevel()); i++) {
    workVecSz[i] = 0;
  }

  for (unsigned int i = 0; i < inp.size(); i++) {
#ifdef __DEBUG_TN__
    assert(areComparable(*this, inp[i]));
#endif
    if (this->isAncestor(inp[i])) {
      workVecSz[maxDepth - (inp[i].getLevel())]++;
    } else {
      std::cout << inp[i] << " is not a decendant of " << (*this) << std::endl;
      assert(false);
    }
  } //end for

  for (unsigned int i = 0; i < (maxDepth - getLevel()); i++) {
    workVec[i].resize(workVecSz[i]);
    workVecSz[i] = 0;
  }

  for (unsigned int i = 0; i < inp.size(); i++) {
#ifdef __DEBUG_TN__
    assert(areComparable(*this, inp[i]));
#endif
    if (this->isAncestor(inp[i])) {
      workVec[maxDepth - (inp[i].getLevel())][workVecSz[maxDepth - (inp[i].getLevel())]] = inp[i];
      workVecSz[maxDepth - (inp[i].getLevel())]++;
    }
  } //end for

  if (workVecSz) {
    delete[] workVecSz;
    workVecSz = NULL;
  }

  TreeNode tmpRoot(m_uiDim, m_uiMaxDepth);
  for (unsigned int i = 0; i < (maxDepth - getLevel() - 1); i++) {
    //This also sorts.
    seq::makeVectorUnique(workVec[i], false);
    //Remove Brothers from a sorted list.
    std::vector<TreeNode> tmpList(workVec[i].size());
    unsigned int tmpListSz = 0;
    if (!workVec[i].empty()) {
      tmpList[tmpListSz] = workVec[i][0];
      tmpListSz++;
    }
    for (unsigned int j = 1; j < workVec[i].size(); j++) {
      if (workVec[i][j - 1].getParent() != workVec[i][j].getParent()) {
        tmpList[tmpListSz] = workVec[i][j];
        tmpListSz++;
      }
    } //end for j
    workVec[i].clear();
    tmpList.resize(tmpListSz);
    workVec[i] = tmpList;
    tmpList.clear();
    unsigned int oldOutSz = out.size();
    unsigned int newOut = 0;
    out.resize(oldOutSz + (workVec[i].size() * (1 << dim)));
    unsigned int oldNextWsz = workVec[i + 1].size();
    unsigned int newNextW = 0;
    unsigned int wLen = workVec[i].size();
    if (wLen > 0) {
      workVec[i + 1].resize(oldNextWsz + (7 * (workVec[i].size())));
    }
    for (unsigned int j = 0; j < wLen; j++) {
      out[oldOutSz + newOut] = workVec[i][j];
      newOut++;
      std::vector<TreeNode> bros;
      workVec[i][j].addBrothers(bros);
      for (int k = 0; k < ((1 << dim) - 1); k++) {
        out[oldOutSz + newOut] = bros[k];
        newOut++;
      } //end for k
      bros.clear();

      TreeNode  parent = workVec[i][j].getParent();
      parent.addBrothers(bros);
      for (int k = 0; k < ((1 << dim) - 1); k++) {
        workVec[i + 1][oldNextWsz + newNextW] = bros[k];
        newNextW++;
      } //end for k
      bros.clear();
      //Add all the brothers of parent to workVec[i+1]
    } //end for j
    workVec[i].clear();
    if (wLen > 0) {
      workVec[i + 1].resize(oldNextWsz + newNextW);
    }
  } //end for i

  if (workVec) {
    if (!(workVec[maxDepth - 1 - getLevel()].empty())) {
      unsigned int oldOutSz = out.size();
      unsigned int newOut = 0;
      out.resize(oldOutSz + (1 << dim));
      out[oldOutSz + newOut] = workVec[maxDepth - 1 - getLevel()][0];
      newOut++;
      std::vector<TreeNode> bros;
      workVec[maxDepth - 1 - getLevel()][0].addBrothers(bros);
      workVec[maxDepth - 1 - getLevel()].clear();
      for (int k = 0; k < ((1 << dim) - 1); k++) {
        out[oldOutSz + newOut] = bros[k];
        newOut++;
      } //end for k

      bros.clear();
    } //end if
    delete[] workVec;
    workVec = NULL;
  }

  //sort and removde dups.
  seq::makeVectorUnique<TreeNode>(out, false);
  //Remove overlaps.
  std::vector<unsigned int> toRem(out.size());
  unsigned int remCtr = 0;
  //compare i+1 with i
  for (unsigned int i = 1; i < out.size(); i++) {
#ifdef __DEBUG_TN__
    assert(areComparable(out[i - 1], out[i]));
#endif
    if (out[i - 1].isAncestor(out[i])) {
      toRem[remCtr] = (i - 1);
      remCtr++;
    } //end if
  } //end for i
  toRem.resize(remCtr);

  if (!toRem.empty()) {
    //Note toRem is sorted.
    std::vector<TreeNode> tmpRem(out.size() - toRem.size());
    unsigned int pos = 0;
    unsigned int kk = 0;
    for (unsigned int i = 0; i < toRem.size(); i++) {
      for (unsigned int j = pos; j < toRem[i]; j++) {
        tmpRem[kk] = out[j];
        kk++;
      } //end for j
      pos = toRem[i] + 1;
    } //end for i
    for (unsigned int j = pos; j < out.size(); j++) {
      tmpRem[kk] = out[j];
      kk++;
    } //end for j
    out = tmpRem;
    tmpRem.clear();
  } //end if
  toRem.clear();

  PROF_COMPLETE_SUBTREE_END

} //end function

} //end namespace
