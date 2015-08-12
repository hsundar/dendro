
/**
  @file TreeNode.txx
  @author Rahul S. Sampath, rahul.sampath@gmail.com
 */

#ifdef __DEBUG__
#ifndef __DEBUG_TN__
#define __DEBUG_TN__
#endif
#endif

#include "colors.h"
#include <cassert>

#define SWAP(a, b) ((&(a) == &(b)) || \
                    (((a) -= (b)), ((b) += (a)), ((a) = (b) - (a))))

namespace ot {
// Inline Functions ...
inline unsigned char TreeNode::getChildNumber() const {
  unsigned int len = (1u << (m_uiMaxDepth - getLevel()));
  unsigned int len_par = (1u << (m_uiMaxDepth - getLevel() + 1u));

  unsigned int i = (m_uiX % len_par);
  unsigned int j = (m_uiY % len_par);
  unsigned int k = (m_uiZ % len_par);
  i /= len;
  j /= len;
  k /= len;

  unsigned char childNum = static_cast<unsigned char>(4 * k + 2 * j + i);
  return childNum;

} //end function

inline int TreeNode::getAnchor(unsigned int& x, unsigned int& y, unsigned int& z) const {
  x = m_uiX;
  y = m_uiY;
  z = m_uiZ;
  return 1;
}


//========== Overloaded Operators =========== //

inline bool TreeNode  :: operator ==(TreeNode   const& other)  const {
#ifdef __DEBUG_TN__
  if (((this->m_uiDim) != (other.m_uiDim)) || ((this->m_uiMaxDepth) != (other.m_uiMaxDepth))) {
    std::cout << "Me: " << (*this) << " Other: " << other << std::endl;
    std::cout << "My Dim: " << m_uiDim << " OthDim: " << other.m_uiDim << " My MaxD: " << m_uiMaxDepth << " othMD: " << other.m_uiMaxDepth << std::endl;
    assert(false);
  }
#endif
  if ((m_uiX == other.m_uiX) && (m_uiY == other.m_uiY) && (m_uiZ == other.m_uiZ) &&
      ((m_uiLevel & ot::TreeNode::MAX_LEVEL) == (other.m_uiLevel & ot::TreeNode::MAX_LEVEL))) {
    return true;
  } else {
    return false;
  }
} //end fn.

inline bool TreeNode  :: operator !=(TreeNode  const& other)  const {
#ifdef __DEBUG_TN__
  if (((this->m_uiDim) != (other.m_uiDim)) || ((this->m_uiMaxDepth) != (other.m_uiMaxDepth))) {
    std::cout << "Me: " << (*this) << " Other: " << other << std::endl;
    std::cout << "My Dim: " << m_uiDim << " OthDim: " << other.m_uiDim << " My MaxD: " << m_uiMaxDepth << " othMD: " << other.m_uiMaxDepth << std::endl;
    assert(false);
  }
#endif
  return (!((*this) == other));
} //end fn.


// The MAIN Comparison Function ...
inline bool TreeNode  :: operator  <(TreeNode   const& other)  const {
#ifdef __DEBUG_TN__
  if (((this->m_uiDim) != (other.m_uiDim)) || ((this->m_uiMaxDepth) != (other.m_uiMaxDepth))) {
    std::cout << "Me: " << (*this) << " Other: " << other << std::endl;
    std::cout << "My Dim: " << m_uiDim << " OthDim: " << other.m_uiDim << " My MaxD: " << m_uiMaxDepth << " othMD: " << other.m_uiMaxDepth << std::endl;
    assert(false);
  }
#endif

  /* -- original Morton
    // first compare the x, y, and z to determine which one dominates ...
    //Ancestor is smaller.
    if ((this->m_uiX == other.m_uiX) && (this->m_uiY == other.m_uiY) && (this->m_uiZ == other.m_uiZ)) {
      return ((this->m_uiLevel & ot::TreeNode::MAX_LEVEL) < (other.m_uiLevel & ot::TreeNode::MAX_LEVEL));
    } //end if

    unsigned int x = (m_uiX ^ other.m_uiX);
    unsigned int y = (m_uiY ^ other.m_uiY);
    unsigned int z = (m_uiZ ^ other.m_uiZ);

    //Default pref: z > y > x.
    unsigned int maxC = z;
    unsigned int yOrx = y;
    if (yOrx < x) {if ((x ^ yOrx) >= yOrx) {yOrx = x;}
    }
    if (maxC < yOrx) {if ((maxC ^ yOrx) >= maxC) {maxC = yOrx;}
    }

    if (maxC == z) {return (m_uiZ < other.m_uiZ); } else if (maxC == y) {return (m_uiY < other.m_uiY); } else {return (m_uiX < other.m_uiX); }
  -- original Morton */

  if ((this->m_uiX == other.m_uiX) && (this->m_uiY == other.m_uiY) && (this->m_uiZ == other.m_uiZ)) {
      return ((this->m_uiLevel & ot::TreeNode::MAX_LEVEL) < (other.m_uiLevel & ot::TreeNode::MAX_LEVEL));
  } //end if

  const Point p1=Point(this->m_uiX,this->m_uiY,this->m_uiZ);
  const Point p2=Point(other.m_uiX,other.m_uiY,other.m_uiZ);
   
#ifdef HILBERT_ORDERING
  #ifdef USE_NCA_PROPERTY
    #pragma message GRN"Hilbert NCA \e[0m"
    return hilbert_order_NCA(p1,p2);
  #else
    #pragma message GRN"Hilbert \e[0m"
    return hilbert_order(p1,p2);
  #endif
#else 
  #ifdef USE_NCA_PROPERTY
    #pragma message(RED"Morton NCA\e[0m")
    return morton_order_NCA(p1,p2);
  #else 
    #pragma message(RED"Morton \e[0m")
    return morton_order(p1,p2);
  #endif
#endif

} //end function

inline bool TreeNode  :: operator  <=(TreeNode  const& other)  const {
#ifdef __DEBUG_TN__
  if (((this->m_uiDim) != (other.m_uiDim)) || ((this->m_uiMaxDepth) != (other.m_uiMaxDepth))) {
    std::cout << "Me: " << (*this) << " Other: " << other << std::endl;
    std::cout << "My Dim: " << m_uiDim << " OthDim: " << other.m_uiDim << " My MaxD: " << m_uiMaxDepth << " othMD: " << other.m_uiMaxDepth << std::endl;
    assert(false);
  }
#endif
  return (((*this) < other) || ((*this) == other));
} //end fn.

inline bool TreeNode  :: operator  >(TreeNode  const& other)  const {
#ifdef __DEBUG_TN__
  if (((this->m_uiDim) != (other.m_uiDim)) || ((this->m_uiMaxDepth) != (other.m_uiMaxDepth))) {
    std::cout << "Me: " << (*this) << " Other: " << other << std::endl;
    std::cout << "My Dim: " << m_uiDim << " OthDim: " << other.m_uiDim << " My MaxD: " << m_uiMaxDepth << " othMD: " << other.m_uiMaxDepth << std::endl;
    assert(false);
  }
#endif
  return ((!((*this) < other)) && (!((*this) == other)));
} //end fn.

inline bool TreeNode  :: operator  >=(TreeNode  const& other)  const {
#ifdef __DEBUG_TN__
  if (((this->m_uiDim) != (other.m_uiDim)) || ((this->m_uiMaxDepth) != (other.m_uiMaxDepth))) {
    std::cout << "Me: " << (*this) << " Other: " << other << std::endl;
    std::cout << "My Dim: " << m_uiDim << " OthDim: " << other.m_uiDim << " My MaxD: " << m_uiMaxDepth << " othMD: " << other.m_uiMaxDepth << std::endl;
    assert(false);
  }
#endif
  return (!((*this) < other));
} //end fn.

//========== Overloaded Operators =========== //

inline TreeNode TreeNode  ::getParent() const {
  //For any node at level l, the last (maxD-l) bits are 0.
  //By convention, root's parent is also root.
  unsigned int parX, parY, parZ;
  unsigned int parLev = (((m_uiLevel & ot::TreeNode::MAX_LEVEL) > 0)
                         ? ((m_uiLevel & ot::TreeNode::MAX_LEVEL) - 1) : 0);
  parX = ((m_uiX >> (m_uiMaxDepth - parLev)) << (m_uiMaxDepth - parLev));
  parY = ((m_uiY >> (m_uiMaxDepth - parLev)) << (m_uiMaxDepth - parLev));
  parZ = ((m_uiZ >> (m_uiMaxDepth - parLev)) << (m_uiMaxDepth - parLev));
  return TreeNode(1, parX, parY, parZ, parLev, m_uiDim, m_uiMaxDepth);
} //end function

inline TreeNode TreeNode  ::getAncestor(unsigned int ancLev) const {
  //For any node at level l, the last (maxD-l) bits are 0.
  unsigned int ancX, ancY, ancZ;
  ancX = ((m_uiX >> (m_uiMaxDepth - ancLev)) << (m_uiMaxDepth - ancLev));
  ancY = ((m_uiY >> (m_uiMaxDepth - ancLev)) << (m_uiMaxDepth - ancLev));
  ancZ = ((m_uiZ >> (m_uiMaxDepth - ancLev)) << (m_uiMaxDepth - ancLev));
  return TreeNode(1, ancX, ancY, ancZ, ancLev, m_uiDim, m_uiMaxDepth);
} //end function

inline int TreeNode::setWeight(unsigned int w) {
  m_uiWeight = w;
  return 1;
}

inline int TreeNode::addWeight(unsigned int w) {
  m_uiWeight += w;
  return 1;
}

inline unsigned int TreeNode::getDim() const {
  return m_uiDim;
}

inline unsigned int TreeNode::getMaxDepth() const {
  return m_uiMaxDepth;
}

inline unsigned int TreeNode::getLevel() const {
  return (m_uiLevel & ot::TreeNode::MAX_LEVEL);
}

inline unsigned int TreeNode::getFlag() const {
  return m_uiLevel;
}

inline int TreeNode::orFlag(unsigned int w) {
  m_uiLevel = (m_uiLevel | w);
  return 1;
}

inline int TreeNode::setFlag(unsigned int w) {
  m_uiLevel = w;
  return 1;
}

inline unsigned int TreeNode::getWeight() const {
  return m_uiWeight;
}

inline unsigned int TreeNode::getX() const {
  return m_uiX;
}

inline unsigned int TreeNode::getY() const {
  return m_uiY;
}

inline unsigned int TreeNode::getZ() const {
  return m_uiZ;
}

inline unsigned int TreeNode::getParentX() const {
  return getParent().getX();
}

inline unsigned int TreeNode::getParentY() const {
  return getParent().getY();
}

inline unsigned int TreeNode::getParentZ() const {
  return getParent().getZ();
}

inline TreeNode TreeNode  ::getDFD() const {
  TreeNode dfd(1, m_uiX, m_uiY, m_uiZ, m_uiMaxDepth, m_uiDim, m_uiMaxDepth);
  return dfd;
} //end function

inline TreeNode TreeNode  ::getDLD() const {
#ifdef HILBERT_ORDERINGs
  TreeNode dld(1, minX() , minY(), maxZ() - 1, m_uiMaxDepth, m_uiDim, m_uiMaxDepth);
#else
  TreeNode dld(1, maxX() - 1, maxY() - 1, maxZ() - 1, m_uiMaxDepth, m_uiDim, m_uiMaxDepth);
#endif
  return dld;
} //end function

inline TreeNode TreeNode::getNext() const {

  TreeNode m = *this;
  unsigned int mask = (1u << (m_uiMaxDepth - getLevel()));

  int i;
  for (i = m.m_uiLevel; i >= 0; i--) {
    m.m_uiX = (m.m_uiX ^ mask);
    if ((m.m_uiX & mask)) break;
    m.m_uiY = (m.m_uiY ^ mask);
    if ((m.m_uiY & mask)) break;
    m.m_uiZ = (m.m_uiZ ^ mask);
    if ((m.m_uiZ & mask)) break;
    mask = (mask << 1);
  }
  m.m_uiLevel = i;
  return m;

}

inline TreeNode TreeNode::getFirstChild() const {
  TreeNode m = *this;
  m.m_uiLevel++;
  return m;
}

inline bool TreeNode :: isAncestor(TreeNode const& other) const {
#ifdef __DEBUG_TN__
  if (((this->m_uiDim) != (other.m_uiDim)) || ((this->m_uiMaxDepth) != (other.m_uiMaxDepth)))  {
    std::cout << "Me: " << (*this) << " Other: " << other << std::endl;
    std::cout << "My Dim: " << m_uiDim << " OthDim: " << other.m_uiDim << " My MaxD: " << m_uiMaxDepth << " othMD: " << other.m_uiMaxDepth << std::endl;
    assert(false);
  }
#endif
/*
#ifdef HILBERT_ORDERING
  unsigned int min1[3], min2[3], max1[3], max2[3];
 
  min1[0] = this->minX(); min1[1] = this->minY(); min1[2] = this->minZ();
  min2[0] = other.minX(); min2[1] = other.minY(); min2[2] = other.minZ();
   
  max1[0] = this->maxX(); max1[1] = this->maxY(); max1[2] = this->maxZ();
  max2[0] = other.maxX(); max2[1] = other.maxY(); max2[2] = other.maxZ();

  return ( (this->m_uiLevel < other.m_uiLevel) && (min2[0] >= min1[0]) && (min2[1] >= min1[1]) && (min2[2] >= min1[2]) && (max2[0] <= max1[0]) && (max2[1] <= max1[1]) && (max2[2] <= max1[2]) );
#else */
return ((other > (*this)) && (other <= (this->getDLD())));
//#endif
} //end function

inline unsigned int TreeNode   :: minX() const {
  return getX();
} //end fn.

inline  unsigned int TreeNode   ::  minY()  const {
  if (m_uiDim < 2) {return 0;}
  return getY();
} //end fn.

inline  unsigned int TreeNode   ::  minZ() const {
  if (m_uiDim < 3) {return 0;}
  return getZ();
} //end fn.

inline  unsigned int TreeNode   ::  maxX() const {
  unsigned int len = (1u << (m_uiMaxDepth - getLevel()));
  return (minX() + len);
} //end fn.

inline  unsigned int TreeNode   ::  maxY() const {
  if (m_uiDim < 2) {return 1;}
  unsigned int len = (1u << (m_uiMaxDepth - getLevel()));
  return (minY() + len);
} //end fn.

inline  unsigned int TreeNode   ::  maxZ() const {
  if (m_uiDim < 3) {return 1;}
  unsigned int len = (1u << (m_uiMaxDepth - getLevel()));
  return (minZ() + len);
} //end fn.

inline  TreeNode   TreeNode   ::  getLeft() const {
  //-ve X
  if (minX() > 0) {
    unsigned int len =  (1u << (m_uiMaxDepth - getLevel()));
    unsigned int xres = (minX() - len);
    unsigned int yres = minY();
    unsigned int zres = minZ();
    TreeNode res(1, xres, yres, zres, getLevel(), m_uiDim, m_uiMaxDepth);
    return res;
  } else {
    TreeNode res(m_uiDim, m_uiMaxDepth);
    return res;
  }
} //end fn.

inline  TreeNode   TreeNode   ::  getRight() const  {
  //+ve X
  if (maxX() < (1u << m_uiMaxDepth)) {
    unsigned int xres = maxX();
    unsigned int yres = minY();
    unsigned int zres = minZ();
    TreeNode  res(1, xres, yres, zres, getLevel(), m_uiDim, m_uiMaxDepth);
    return res;
  } else {
    TreeNode res(m_uiDim, m_uiMaxDepth);
    return res;
  }
} //end fn.

inline  TreeNode   TreeNode   ::  getBack() const {
  //+ve Y
  if ((m_uiDim > 1) && (maxY() < (1u << m_uiMaxDepth))) {
    unsigned int xres = minX();
    unsigned int yres = maxY();
    unsigned int zres = minZ();
    TreeNode res(1, xres, yres, zres, getLevel(), m_uiDim, m_uiMaxDepth);
    return res;
  } else {
    TreeNode res(m_uiDim, m_uiMaxDepth);
    return res;
  }
} //end fn.

inline TreeNode   TreeNode   ::  getFront() const {
  //-ve Y
  if (minY() > 0) {
    unsigned int len =  (1u << (m_uiMaxDepth - getLevel()));
    unsigned int xres = minX();
    unsigned int yres = minY() - len;
    unsigned int zres = minZ();
    TreeNode res(1, xres, yres, zres, getLevel(), m_uiDim, m_uiMaxDepth);
    return res;
  } else {
    TreeNode res(m_uiDim, m_uiMaxDepth);
    return res;
  }
} //end fn.

inline TreeNode   TreeNode   ::  getBottom() const {
  //-ve Z
  if (minZ() > 0) {
    unsigned int len = (1u << (m_uiMaxDepth - getLevel()));
    unsigned int xres = minX();
    unsigned int yres = minY();
    unsigned int zres = minZ() - len;
    TreeNode res(1, xres, yres, zres, getLevel(), m_uiDim, m_uiMaxDepth);
    return res;
  } else {
    TreeNode res(m_uiDim, m_uiMaxDepth);
    return res;
  }
} //end fn.

inline  TreeNode   TreeNode   ::  getTop() const {
  //+ve Z
  if ((m_uiDim == 3) && (maxZ() < (1u << m_uiMaxDepth))) {
    unsigned int xres = minX();
    unsigned int yres = minY();
    unsigned int zres = maxZ();
    TreeNode res(1, xres, yres, zres, getLevel(), m_uiDim, m_uiMaxDepth);
    return res;
  } else {
    TreeNode res(m_uiDim, m_uiMaxDepth);
    return res;
  }
} //end fn.

inline  TreeNode   TreeNode   ::  getLeftBack() const {
  return (getLeft().getBack());
} //end fn.

inline  TreeNode    TreeNode   :: getRightBack()  const {
  return (getRight().getBack());
} //end fn.

inline  TreeNode   TreeNode   ::  getLeftFront() const {
  return (getLeft().getFront());
} //end fn.

inline TreeNode   TreeNode   ::  getRightFront() const {
  return (getRight().getFront());
} //end fn.

inline TreeNode   TreeNode   ::  getBottomLeft() const {
  return (getBottom().getLeft());
} //end fn.

inline  TreeNode   TreeNode   ::  getBottomRight() const {
  return (getBottom().getRight());
} //end fn.

inline  TreeNode   TreeNode   ::  getBottomBack() const {
  return (getBottom().getBack());
} //end fn.

inline TreeNode   TreeNode   ::  getBottomFront() const {
  return (getBottom().getFront());
} //end fn.

inline TreeNode   TreeNode   ::  getBottomLeftBack() const {
  return (getBottomLeft().getBack());
} //end fn.

inline  TreeNode   TreeNode   ::  getBottomRightBack() const {
  return (getBottomRight().getBack());
} //end fn.

inline  TreeNode   TreeNode   ::  getBottomLeftFront() const {
  return (getBottomLeft().getFront());
} //end fn.

inline TreeNode   TreeNode   ::  getBottomRightFront()  const {
  return (getBottomRight().getFront());
} //end fn.

inline  TreeNode   TreeNode   ::  getTopLeft() const {
  return (getTop().getLeft());
} //end fn.

inline   TreeNode    TreeNode   :: getTopRight() const {
  return (getTop().getRight());
} //end fn.

inline TreeNode    TreeNode   :: getTopBack() const {
  return (getTop().getBack());
} //end fn.

inline  TreeNode   TreeNode   ::  getTopFront() const {
  return (getTop().getFront());
} //end fn.

inline  TreeNode   TreeNode   ::  getTopLeftBack() const {
  return (getTopLeft().getBack());
} //end fn.

inline TreeNode    TreeNode   :: getTopRightBack() const {
  return (getTopRight().getBack());
} //end fn.

inline  TreeNode   TreeNode   ::  getTopLeftFront() const {
  return (getTopLeft().getFront());
} //end fn.

inline TreeNode   TreeNode   ::  getTopRightFront() const {
  return (getTopRight().getFront());
} //end fn.

///////////////////////////////////////////////////////////////////////////////////////
//Returns a Vector of Candidates for Balanced Neigbors in a particular direction.
//The results are returned in the following order: 1) Same level, 2) Parent's
//level and 3) Children's level.
//If there is no candidate at the parent's or child's level, then the candidate at the
//same level is returned instead for that position
//Note: For the neighbours at the parent's level, the dimension is not
//considered explicitly. It is implicit in the childNum.

inline std::vector<TreeNode>  TreeNode   :: getB_Left() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "L" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(6);
  } else if (dim == 2) {
    it.resize(4);
  } else {
    it.resize(3);
  }

  it[0] = getLeft();
  unsigned int childNum = getChildNumber();

  //0,2,4,6
  if ((childNum == 0) || (childNum == 2) || (childNum == 4) || (childNum == 6)) {
    it[1] = getParent().getLeft();
  } else {
    it[1] = it[0];
  }

  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[0].getLeft();
  if (dim > 1) {
    it[3] = children[2].getLeft();
  }
  if (dim > 2) {
    it[4] = children[4].getLeft();
    it[5] = children[6].getLeft();
  }
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_Right() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "R" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(6);
  } else if (dim == 2) {
    it.resize(4);
  } else {
    it.resize(3);
  }
  //1,3,5,7
  it[0] = getRight();
  unsigned int childNum = getChildNumber();
  if ((childNum == 1) || (childNum == 3) || (childNum == 5) || (childNum == 7)) {
    it[1] = getParent().getRight();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[1].getRight();
  if (dim > 1) {
    it[3] = children[3].getRight();
  }
  if (dim > 2) {
    it[4] = children[5].getRight();
    it[5] = children[7].getRight();
  }
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_Top() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "T" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(6);
  } else if (dim == 2) {
    it.resize(0);
    return it;
  } else {
    it.resize(0);
    return it;
  }
  //4,5,6,7
  it[0] = getTop();
  unsigned int childNum = getChildNumber();
  if ((childNum  == 4) || (childNum == 5) || (childNum  == 6) || (childNum == 7)) {
    it[1] = getParent().getTop();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[4].getTop();
  it[3] = children[5].getTop();
  it[4] = children[6].getTop();
  it[5] = children[7].getTop();
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_Bottom() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "Bo" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(6);
  } else if (dim == 2) {
    it.resize(0);
    return it;
  } else {
    it.resize(0);
    return it;
  }
  //0,1,2,3
  it[0] = getBottom();
  unsigned int childNum = getChildNumber();
  if ((childNum == 0) || (childNum == 1) || (childNum == 2) || (childNum == 3)) {
    it[1] = getParent().getBottom();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[0].getBottom();
  it[3] = children[1].getBottom();
  it[4] = children[2].getBottom();
  it[5] = children[3].getBottom();
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_Front() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "F" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(6);
  } else if (dim == 2) {
    it.resize(4);
  } else {
    it.resize(0);
    return it;
  }
  //0,1,4,5
  it[0] = getFront();
  unsigned int childNum = getChildNumber();
  if ((childNum == 0) || (childNum == 1) || (childNum == 4) || (childNum == 5)) {
    it[1] = getParent().getFront();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[0].getFront();
  it[3] = children[1].getFront();
  if (dim > 2) {
    it[4] = children[4].getFront();
    it[5] = children[5].getFront();
  }
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_Back() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "Bk" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(6);
  } else if (dim == 2) {
    it.resize(4);
  } else {
    it.resize(0);
    return it;
  }
  //2,3,6,7
  it[0] = getBack();
  unsigned int childNum = getChildNumber();
  if ((childNum == 2) || (childNum == 3) || (childNum == 6) || (childNum == 7)) {
    it[1] = getParent().getBack();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[2].getBack();
  it[3] = children[3].getBack();
  if (dim == 3) {
    it[4] = children[6].getBack();
    it[5] = children[7].getBack();
  }
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_TopLeft() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "TL" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(4);
  } else if (dim == 2) {
    it.resize(0);
    return it;
  } else {
    it.resize(0);
    return it;
  }
  //4,6 -> TL
  //5,7->T but not L
  //0,2->L but not T
  it[0] = getTopLeft();
  unsigned int childNum = getChildNumber();
  if ((childNum == 4) || (childNum == 6)) {
    it[1] = getParent().getTopLeft();
  } else if ((childNum == 5) || (childNum == 7)) {
    it[1] = getParent().getTop();
  } else if ((childNum == 0) || (childNum == 2)) {
    it[1] = getParent().getLeft();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[4].getTopLeft();
  it[3] = children[6].getTopLeft();
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_TopRight() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "TR" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(4);
  } else if (dim == 2) {
    it.resize(0);
    return it;
  } else {
    it.resize(0);
    return it;
  }
  //5,7 TR
  //4,6 T but not R
  //1,3 R but not T
  it[0] = getTopRight();
  unsigned int childNum = getChildNumber();
  if ((childNum == 5) || (childNum == 7)) {
    it[1] = getParent().getTopRight();
  } else if ((childNum == 4) || (childNum == 6)) {
    it[1] = getParent().getTop();
  } else if ((childNum == 1) || (childNum == 3)) {
    it[1] = getParent().getRight();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[5].getTopRight();
  it[3] = children[7].getTopRight();
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_BottomLeft() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "BoL" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(4);
  } else if (dim == 2) {
    it.resize(0);
    return it;
  } else {
    it.resize(0);
    return it;
  }
  //0,2 BL
  //4,6 L but not B
  //1,3 B but not L
  it[0] = getBottomLeft();
  unsigned int childNum = getChildNumber();
  if ((childNum == 0) || (childNum == 2)) {
    it[1] = getParent().getBottomLeft();
  } else if ((childNum == 1) || (childNum == 3)) {
    it[1] = getParent().getBottom();
  } else if ((childNum == 4) || (childNum == 6)) {
    it[1] = getParent().getLeft();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[0].getBottomLeft();
  it[3] = children[2].getBottomLeft();
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_BottomRight() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "BoR" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(4);
  } else if (dim == 2) {
    it.resize(0);
    return it;
  } else {
    it.resize(0);
    return it;
  }
  //1,3 BR
  //0,2 B but not R
  //5,7 R but not B
  it[0] = getBottomRight();
  unsigned int childNum = getChildNumber();
  if ((childNum == 1) || (childNum == 3)) {
    it[1] = getParent().getBottomRight();
  } else if ((childNum == 0) || (childNum == 2)) {
    it[1] = getParent().getBottom();
  } else if ((childNum == 5) || (childNum == 7)) {
    it[1] = getParent().getRight();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[1].getBottomRight();
  it[3] = children[3].getBottomRight();
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_LeftFront() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "LF" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(4);
  } else if (dim == 2) {
    it.resize(3);
  } else {
    it.resize(0);
    return it;
  }
  //0,4 LF
  //1,5 F but not L
  //2,6 L but not F
  it[0] = getLeftFront();
  unsigned int childNum = getChildNumber();

  if ((childNum == 0) || (childNum == 4)) {
    it[1] = getParent().getLeftFront();
  } else if ((childNum == 1) || (childNum == 5)) {
    it[1] = getParent().getFront();
  } else if ((childNum == 2) || (childNum == 6)) {
    it[1] = getParent().getLeft();
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[0].getLeftFront();
  if (dim > 2) {
    it[3] = children[4].getLeftFront();
  }
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_RightFront() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "RF" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(4);
  } else if (dim == 2) {
    it.resize(3);
  } else {
    it.resize(0);
    return it;
  }
  //1,5 RF
  //0,4 F but not R
  //3,7 R but not F
  it[0] = getRightFront();
  unsigned int childNum = getChildNumber();
  if ((childNum == 1) || (childNum == 5)) {
    it[1] = getParent().getRightFront();
  } else if ((childNum == 0) || (childNum == 4)) {
    it[1] = getParent().getFront();
  } else if ((childNum == 3) || (childNum == 7)) {
    it[1] = getParent().getRight();
  }

  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[1].getRightFront();
  if (dim > 2) {
    it[3] = children[5].getRightFront();
  }
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_TopFront() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "TF" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(4);
  } else if (dim == 2) {
    it.resize(0);
    return it;
  } else {
    it.resize(0);
    return it;
  }
  //4,5 TF
  //6,7 T but not F
  //0,1 F but not T
  it[0] = getTopFront();
  unsigned int childNum = getChildNumber();
  if ((childNum == 4) || (childNum == 5)) {
    it[1] = getParent().getTopFront();
  } else if ((childNum == 6) || (childNum == 7)) {
    it[1] = getParent().getTop();
  } else if ((childNum == 0) || (childNum == 1)) {
    it[1] = getParent().getFront();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[4].getTopFront();
  it[3] = children[5].getTopFront();
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_BottomFront() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "BoF" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(4);
  } else if (dim == 2) {
    it.resize(0);
    return it;
  } else {
    it.resize(0);
    return it;
  }
  //0,1 BoF
  //2,3 Bo but not F
  //4,5 F but not Bo
  it[0] = getBottomFront();
  unsigned int childNum = getChildNumber();
  if ((childNum == 0) || (childNum == 1)) {
    it[1] = getParent().getBottomFront();
  } else if ((childNum == 2) || (childNum == 3)) {
    it[1] = getParent().getBottom();
  } else if ((childNum == 4) || (childNum == 5)) {
    it[1] = getParent().getFront();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[0].getBottomFront();
  it[3] = children[1].getBottomFront();
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_TopLeftFront() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "TLF" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(3);
  } else if (dim == 2) {
    it.resize(0);
    return it;
  } else {
    it.resize(0);
    return it;
  }
  //4 TLF
  //6 TL but not F
  //0 LF but not T
  //5 TF but not L
  //7 T but not LF
  //2 L but not TF
  //1 F but not TL
  it[0] = getTopLeftFront();
  unsigned int childNum = getChildNumber();
  if ((childNum == 4)) {
    it[1] = getParent().getTopLeftFront();
  } else if ((childNum == 6)) {
    it[1] = getParent().getTopLeft();
  } else if ((childNum == 0)) {
    it[1] = getParent().getLeftFront();
  } else if ((childNum == 5)) {
    it[1] = getParent().getTopFront();
  } else if ((childNum == 7)) {
    it[1] = getParent().getTop();
  } else if ((childNum == 2)) {
    it[1] = getParent().getLeft();
  } else if ((childNum == 1)) {
    it[1] = getParent().getFront();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[4].getTopLeftFront();
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_TopRightFront() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "TRF" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(3);
  } else if (dim == 2) {
    it.resize(0);
    return it;
  } else {
    it.resize(0);
    return it;
  }
  //5 TRF
  //7 TR but not F
  //1 RF but not T
  //4 TF but not R
  //6 T but not RF
  //3 R but not TF
  //0 F but not TR
  it[0] = getTopRightFront();
  unsigned int childNum = getChildNumber();
  if ((childNum == 5)) {
    it[1] = getParent().getTopRightFront();
  } else if ((childNum == 7)) {
    it[1] = getParent().getTopRight();
  } else if ((childNum == 1)) {
    it[1] = getParent().getRightFront();
  } else if ((childNum == 4)) {
    it[1] = getParent().getTopFront();
  } else if ((childNum == 6)) {
    it[1] = getParent().getTop();
  } else if ((childNum == 3)) {
    it[1] = getParent().getRight();
  } else if ((childNum == 0)) {
    it[1] = getParent().getFront();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[5].getTopRightFront();
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_BottomLeftFront() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "BoLF" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(3);
  } else if (dim == 2) {
    it.resize(0);
    return it;
  } else {
    it.resize(0);
    return it;
  }
  //0 BoLF
  //2 BoL but not F
  //4 LF but not Bo
  //1 BoF but not L
  //3 Bo but not LF
  //6 L but not BoF
  //5 F but not BoL
  it[0] = getBottomLeftFront();
  unsigned int childNum = getChildNumber();
  if ((childNum == 0)) {
    it[1] = getParent().getBottomLeftFront();
  } else if ((childNum == 2)) {
    it[1] = getParent().getBottomLeft();
  } else if ((childNum == 4)) {
    it[1] = getParent().getLeftFront();
  } else if ((childNum == 1)) {
    it[1] = getParent().getBottomFront();
  } else if ((childNum == 3)) {
    it[1] = getParent().getBottom();
  } else if ((childNum == 6)) {
    it[1] = getParent().getLeft();
  } else if ((childNum == 5)) {
    it[1] = getParent().getFront();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[0].getBottomLeftFront();
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_BottomRightFront() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "BoRF" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(3);
  } else if (dim == 2) {
    it.resize(0);
    return it;
  } else {
    it.resize(0);
    return it;
  }
  //1 BoRF
  //3 BoR
  //5 RF
  //0 BoF
  //2 Bo
  //7 R
  //4 F
  it[0] = getBottomRightFront();
  unsigned int childNum = getChildNumber();
  if ((childNum == 1)) {
    it[1] = getParent().getBottomRightFront();
  } else if ((childNum == 3)) {
    it[1] = getParent().getBottomRight();
  } else if ((childNum == 5)) {
    it[1] = getParent().getRightFront();
  } else if ((childNum == 0)) {
    it[1] = getParent().getBottomFront();
  } else if ((childNum == 2)) {
    it[1] = getParent().getBottom();
  } else if ((childNum == 7)) {
    it[1] = getParent().getRight();
  } else if ((childNum == 4)) {
    it[1] = getParent().getFront();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[1].getBottomRightFront();
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_LeftBack() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "LBk" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(4);
  } else if (dim == 2) {
    it.resize(3);
  } else {
    it.resize(0);
    return it;
  }
  //2,6 LBk
  //0,4 L not Bk
  //3,7 Bk not L
  it[0] = getLeftBack();
  unsigned int childNum = getChildNumber();
  if ((childNum == 2) || (childNum == 6)) {
    it[1] = getParent().getLeftBack();
  } else if ((childNum == 0) || (childNum == 4)) {
    it[1] = getParent().getLeft();
  } else if ((childNum == 3) || (childNum == 7)) {
    it[1] = getParent().getBack();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[2].getLeftBack();
  if (dim == 3) {
    it[3] = children[6].getLeftBack();
  }
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_RightBack() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "RBk" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(4);
  } else if (dim == 2) {
    it.resize(3);
  } else {
    it.resize(0);
    return it;
  }
  //3,7 RBk
  //1,5 R
  //2,6 Bk
  it[0] = getRightBack();
  unsigned int childNum = getChildNumber();
  if ((childNum == 3) || (childNum == 7)) {
    it[1] = getParent().getRightBack();
  } else if ((childNum == 1) || (childNum == 5)) {
    it[1] = getParent().getRight();
  } else if ((childNum == 2) || (childNum == 6)) {
    it[1] = getParent().getBack();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[3].getRightBack();
  if (dim == 3) {
    it[3] = children[7].getRightBack();
  }
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_TopBack() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "TBk" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(4);
  } else if (dim == 2) {
    it.resize(0);
    return it;
  } else {
    it.resize(0);
    return it;
  }
  //6,7 TBk
  //4,5 T
  //2,3 Bk
  it[0] = getTopBack();
  unsigned int childNum = getChildNumber();
  if ((childNum == 6) || (childNum == 7)) {
    it[1] = getParent().getTopBack();
  } else if ((childNum == 4) || (childNum == 5)) {
    it[1] = getParent().getTop();
  } else if ((childNum == 2) || (childNum == 3)) {
    it[1] = getParent().getBack();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[6].getTopBack();
  it[3] = children[7].getTopBack();
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_BottomBack() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "BoBk" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(4);
  } else if (dim == 2) {
    it.resize(0);
    return it;
  } else {
    it.resize(0);
    return it;
  }
  //2,3 BoBk
  //0,1 Bo
  //6,7 Bk
  it[0] = getBottomBack();
  unsigned int childNum = getChildNumber();
  if ((childNum == 2) || (childNum == 3)) {
    it[1] = getParent().getBottomBack();
  } else if ((childNum == 0) || (childNum == 1)) {
    it[1] = getParent().getBottom();
  } else if ((childNum == 6) || (childNum == 7)) {
    it[1] = getParent().getBack();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[2].getBottomBack();
  it[3] = children[3].getBottomBack();
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_TopLeftBack() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "TLBk" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(3);
  } else if (dim == 2) {
    it.resize(0);
    return it;
  } else {
    it.resize(0);
    return it;
  }
  //6 TLBk
  //4 TL
  //2 LBk
  //7 TBk
  //5 T
  //0 L
  //3 Bk
  it[0] = getTopLeftBack();
  unsigned int childNum = getChildNumber();
  if ((childNum == 6)) {
    it[1] = getParent().getTopLeftBack();
  } else if ((childNum == 4)) {
    it[1] = getParent().getTopLeft();
  } else if ((childNum == 2)) {
    it[1] = getParent().getLeftBack();
  } else if ((childNum == 7)) {
    it[1] = getParent().getTopBack();
  } else if ((childNum == 5)) {
    it[1] = getParent().getTop();
  } else if ((childNum == 0)) {
    it[1] = getParent().getLeft();
  } else if ((childNum == 3)) {
    it[1] = getParent().getBack();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[6].getTopLeftBack();
  children.clear();
  return it;
} //end fn.

inline std::vector<TreeNode>  TreeNode   :: getB_TopRightBack() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "TRBk" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(3);
  } else if (dim == 2) {
    it.resize(0);
    return it;
  } else {
    it.resize(0);
    return it;
  }
  //7 TRBk
  //5 TR
  //3 RBk
  //6 TBk
  //4 T
  //1 R
  //2 Bk
  it[0] = getTopRightBack();
  unsigned int childNum = getChildNumber();
  if ((childNum == 7)) {
    it[1] = getParent().getTopRightBack();
  } else if ((childNum == 5)) {
    it[1] = getParent().getTopRight();
  } else if ((childNum == 3)) {
    it[1] = getParent().getRightBack();
  } else if ((childNum == 6)) {
    it[1] = getParent().getTopBack();
  } else if ((childNum == 4)) {
    it[1] = getParent().getTop();
  } else if ((childNum == 1)) {
    it[1] = getParent().getRight();
  } else if ((childNum == 2)) {
    it[1] = getParent().getBack();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[7].getTopRightBack();
  children.clear();
  return it;
} //end fn.


inline std::vector<TreeNode>  TreeNode   :: getB_BottomLeftBack() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "BoLBk" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(3);
  } else if (dim == 2) {
    it.resize(0);
    return it;
  } else {
    it.resize(0);
    return it;
  }
  //2 BoLBk
  //0 BoL
  //6 LBk
  //3 BoBk
  //1 Bo
  //4 L
  //7 Bk
  it[0] = getBottomLeftBack();
  unsigned int childNum = getChildNumber();
  if ((childNum == 2)) {
    it[1] = getParent().getBottomLeftBack();
  } else if ((childNum == 0)) {
    it[1] = getParent().getBottomLeft();
  } else if ((childNum == 6)) {
    it[1] = getParent().getLeftBack();
  } else if ((childNum == 3)) {
    it[1] = getParent().getBottomBack();
  } else if ((childNum == 1)) {
    it[1] = getParent().getBottom();
  } else if ((childNum == 4)) {
    it[1] = getParent().getLeft();
  } else if ((childNum == 7)) {
    it[1] = getParent().getBack();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[2].getBottomLeftBack();
  children.clear();
  return it;
} //end fn.

inline std::vector<TreeNode>  TreeNode   :: getB_BottomRightBack() const {
  unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
  std::cout << RED << "BoRBk" << NRM << std::endl;
#endif
  std::vector<TreeNode> it;
  if (dim == 3) {
    it.resize(3);
  } else if (dim == 2) {
    it.resize(0);
    return it;
  } else {
    it.resize(0);
    return it;
  }
  //3 BoRBk
  //1 BoR
  //7 RBk
  //2 BoBk
  //0 Bo
  //5 R
  //6 Bk
  it[0] = getBottomRightBack();
  unsigned int childNum = getChildNumber();
  if ((childNum == 3)) {
    it[1] = getParent().getBottomRightBack();
  } else if ((childNum == 1)) {
    it[1] = getParent().getBottomRight();
  } else if ((childNum == 7)) {
    it[1] = getParent().getRightBack();
  } else if ((childNum == 2)) {
    it[1] = getParent().getBottomBack();
  } else if ((childNum == 0)) {
    it[1] = getParent().getBottom();
  } else if ((childNum == 5)) {
    it[1] = getParent().getRight();
  } else if ((childNum == 6)) {
    it[1] = getParent().getBack();
  } else {
    it[1] = it[0];
  }
  std::vector<TreeNode> children;
  addChildren(children);
  it[2] = children[3].getBottomRightBack();
  children.clear();
  return it;
} //end fn.

//{  Hilbert related changes ===================

inline void TreeNode::rotate(int index, int *current, int *rot_index, int dim) const {

  if (dim == 2) {
    index = rot_index[index];
    if (index == 0) {
      rot_index[current[1]] = 3;
      rot_index[current[3]] = 1;
      SWAP(current[1], current[3]); // RIGHT Rotate and flip orientation


    } else if (index == 3) {
      rot_index[current[0]] = 2;
      rot_index[current[2]] = 0;
      SWAP(current[0], current[2]); //LEFT Rotate and flip orientation:

    }

  } else if (dim == 3) {

    index = rot_index[index];
    if (index == 0) {
      rot_index[current[1]] = 7;
      rot_index[current[7]] = 1;
      SWAP(current[1], current[7]);
      rot_index[current[4]] = 2;
      rot_index[current[2]] = 4;
      SWAP(current[2], current[4]);

    } else if (index == 1) {
      rot_index[current[3]] = 7;
      rot_index[current[7]] = 3;
      SWAP(current[3], current[7]);
      rot_index[current[2]] = 6;
      rot_index[current[6]] = 2;
      SWAP(current[2], current[6]);

    } else if (index == 3) {
      rot_index[current[3]] = 5;
      rot_index[current[5]] = 3;
      SWAP(current[3], current[5]);

      rot_index[current[3]] = 7;
      rot_index[current[7]] = 3;
      SWAP(current[3], current[7]);

      rot_index[current[2]] = 6;
      rot_index[current[6]] = 2;
      SWAP(current[2], current[6]);

      rot_index[current[0]] = 2;
      rot_index[current[2]] = 0;
      SWAP(current[0], current[2]);
    } else if (index == 4) {
      rot_index[current[1]] = 7;
      rot_index[current[7]] = 1;
      SWAP(current[1], current[7]);

      rot_index[current[1]] = 5;
      rot_index[current[5]] = 1;
      SWAP(current[1], current[5]);

      rot_index[current[0]] = 4;
      rot_index[current[4]] = 0;
      SWAP(current[0], current[4]);

      rot_index[current[0]] = 2;
      rot_index[current[2]] = 0;
      SWAP(current[0], current[2]);

    } else if (index == 6) {
      rot_index[current[1]] = 5;
      rot_index[current[5]] = 1;
      SWAP(current[1], current[5]);

      rot_index[current[0]] = 4;
      rot_index[current[4]] = 0;
      SWAP(current[0], current[4]);

    } else if (index == 7) {

      rot_index[current[0]] = 6;
      rot_index[current[6]] = 0;
      SWAP(current[0], current[6]);

      rot_index[current[3]] = 5;
      rot_index[current[5]] = 3;
      SWAP(current[3], current[5]);
    }

  }

}

inline int findIndex(Point *pt, int x, int y, int z, int len) {
  for (int i = 0; i < len; i++) {
    if (pt[i].xint() == x && pt[i].yint() == y && pt[i].zint() == z) {
      return i;
    }
  }

}


inline bool TreeNode::hilbert_order(const Point& p1, const Point& p2) const {

  int MAX_DEAPTH = this->m_uiMaxDepth;
  int MAX_LIMIT = (1 << MAX_DEAPTH) - 1;

  int g_dim = this->m_uiDim;

  int x1 = p1.xint();
  int y1 = p1.yint();
  int z1 = p1.zint();

  int x2 = p2.xint();
  int y2 = p2.yint();
  int z2 = p2.zint();

  if (x1 == x2 && y1 == y2 && z1 == z2) {
    return false;
  }

  int index1 = 0;
  int index2 = 0;
  int min_x, min_y, min_z, max_x, max_y, max_z;

  int len = MAX_LIMIT + 1;
  int deapth = 0;
  min_x = 0;
  min_y = 0;
  min_z = 0;

  max_x = len;
  max_y = len;
  max_z = len;



  if (g_dim == 2) {
    Point pt_hilbert[4];
    Point pt_hilbert_new[4];

    pt_hilbert[0] = Point((int)min_x, (int)min_y, (int)0);
    pt_hilbert[1] = Point((int)min_x, (int)max_y, (int)0);
    pt_hilbert[2] = Point((int)max_x, (int)max_y, (int)0);
    pt_hilbert[3] = Point((int)max_x, (int)min_y, (int)0);



    while (len > 1 && deapth < MAX_DEAPTH) {

      int xl = pt_hilbert[0].xint();
      int yl = pt_hilbert[0].yint();
      int nca_index = 0;
      index1 = 0;
      index2 = 0;
      for (int i = 1; i < 4; i++) {

        if (xl > pt_hilbert[i].xint()) {
          xl = pt_hilbert[i].xint();
          yl = pt_hilbert[i].yint();
          nca_index = i;
        } else if (xl == pt_hilbert[i].xint()) {
          if (yl > pt_hilbert[i].yint()) {
            yl = pt_hilbert[i].yint();
            nca_index = i;
          }
        }


      }
      int len_nca = len / 2;
      // Checking the membership cell.
      if ((x1 - xl) < len_nca && (y1 - yl) < len_nca) {
        index1 = findIndex(pt_hilbert, xl, yl, 0, 4);
      } else if ((x1 - xl) < len_nca && (y1 - yl) >= len_nca) {
        index1 = findIndex(pt_hilbert, xl, yl + len, 0, 4);;
      } else if ((x1 - xl) >= len_nca && (y1 - yl) >= len_nca) {
        index1 = findIndex(pt_hilbert, xl + len, yl + len, 0, 4);;
      } else if ((x1 - xl) >= len_nca && (y1 - yl) < len_nca) {
        index1 = findIndex(pt_hilbert, xl + len, yl, 0, 4);;
      }

      if ((x2 - xl) < len_nca && (y2 - yl) < len_nca) {
        index2 = findIndex(pt_hilbert, xl, yl, 0, 4);
      } else if ((x2 - xl) < len_nca && (y2 - yl) >= len_nca) {
        index2 = findIndex(pt_hilbert, xl, yl + len, 0, 4);
      } else if ((x2 - xl) >= len_nca && (y2 - yl) >= len_nca) {
        index2 = findIndex(pt_hilbert, xl + len, yl + len, 0, 4);
      } else if ((x2 - xl) >= len_nca && (y2 - yl) < len_nca) {
        index2 = findIndex(pt_hilbert, xl + len, yl, 0, 4);
      }

      if (index1 < index2) {
        return true;
      } else if (index1 > index2) {
        return false;
      }

      switch (index1) {
      case 0:
        pt_hilbert_new[0] = pt_hilbert[0];
        pt_hilbert_new[1] = (pt_hilbert[0] + pt_hilbert[3]) / 2;
        pt_hilbert_new[2] = (pt_hilbert[0] + pt_hilbert[2]) / 2;
        pt_hilbert_new[3] = (pt_hilbert[0] + pt_hilbert[1]) / 2;
        break;
      case 1:
        pt_hilbert_new[0] = (pt_hilbert[0] + pt_hilbert[1]) / 2;
        pt_hilbert_new[1] = pt_hilbert[1];
        pt_hilbert_new[2] = (pt_hilbert[1] + pt_hilbert[2]) / 2;
        pt_hilbert_new[3] = (pt_hilbert[1] + pt_hilbert[3]) / 2;
        break;
      case 2:
        pt_hilbert_new[0] = (pt_hilbert[0] + pt_hilbert[2]) / 2;
        pt_hilbert_new[1] = (pt_hilbert[1] + pt_hilbert[2]) / 2;
        pt_hilbert_new[2] = pt_hilbert[2];
        pt_hilbert_new[3] = (pt_hilbert[2] + pt_hilbert[3]) / 2;
        break;
      case 3:
        pt_hilbert_new[0] = (pt_hilbert[3] + pt_hilbert[2]) / 2;
        pt_hilbert_new[1] = (pt_hilbert[1] + pt_hilbert[3]) / 2;
        pt_hilbert_new[2] = (pt_hilbert[3] + pt_hilbert[0]) / 2;
        pt_hilbert_new[3] = pt_hilbert[3];
        break;
      default:
        std::cout << "Hilbert ordering error:Invalid nearest cell" << std::endl;
        break;

      }


      pt_hilbert[0] = pt_hilbert_new[0];
      pt_hilbert[1] = pt_hilbert_new[1];
      pt_hilbert[2] = pt_hilbert_new[2];
      pt_hilbert[3] = pt_hilbert_new[3];
      len = len / 2;
      deapth += 1;

    }
    return false;
  } else if (g_dim == 3) {
    Point pt_hilbert[8];
    Point pt_hilbert_new[8];

    pt_hilbert[0] = Point((int)min_x, (int)min_y, (int)min_z);
    pt_hilbert[1] = Point((int)min_x, (int)max_y, (int)min_z);
    pt_hilbert[2] = Point((int)max_x, (int)max_y, (int)min_z);
    pt_hilbert[3] = Point((int)max_x, (int)min_y, (int)min_z);

    pt_hilbert[4] = Point((int)max_x, (int)min_y, (int)max_z);
    pt_hilbert[5] = Point((int)max_x, (int)max_y, (int)max_z);
    pt_hilbert[6] = Point((int)min_x, (int)max_y, (int)max_z);
    pt_hilbert[7] = Point((int)min_x, (int)min_y, (int)max_z);


    while (len > 1 && deapth < MAX_DEAPTH) {


      int xl = pt_hilbert[0].xint();
      int yl = pt_hilbert[0].yint();
      int zl = pt_hilbert[0].zint();
      int nca_index = 0;
      index1 = 0;
      index2 = 0;
      for (int i = 1; i < 8; i++) {

        if (xl > pt_hilbert[i].xint()) {
          xl = pt_hilbert[i].xint();
          yl = pt_hilbert[i].yint();
          nca_index = i;
        } else if (xl == pt_hilbert[i].xint()) {
          if (yl > pt_hilbert[i].yint()) {
            yl = pt_hilbert[i].yint();
            nca_index = i;
          } else if (yl == pt_hilbert[i].yint()) {
            if (zl > pt_hilbert[i].zint()) {
              zl = pt_hilbert[i].zint();
              nca_index = i;
            }
          }

        }


      }
      int len_nca = len / 2;

      if ((x1 - xl) < len_nca && (y1 - yl) < len_nca && (z1 - zl) < len_nca) {
        index1 = findIndex(pt_hilbert, xl, yl, zl, 8);

      } else if ((x1 - xl) < len_nca && (y1 - yl) >= len_nca && (z1 - zl) < len_nca) {
        index1 = findIndex(pt_hilbert, xl, yl + len, zl, 8);
      } else if ((x1 - xl) >= len_nca && (y1 - yl) >= len_nca && (z1 - zl) < len_nca) {
        index1 = findIndex(pt_hilbert, xl + len, yl + len, zl, 8);
      } else if ((x1 - xl) >= len_nca && (y1 - yl) < len_nca && (z1 - zl) < len_nca) {
        index1 = findIndex(pt_hilbert, xl + len, yl, zl, 8);
      } else if ((x1 - xl) >= len_nca && (y1 - yl) < len_nca && (z1 - zl) >= len_nca) {
        index1 = findIndex(pt_hilbert, xl + len, yl, zl + len, 8);
      } else if ((x1 - xl) >= len_nca && (y1 - yl) >= len_nca && (z1 - zl) >= len_nca) {
        index1 = findIndex(pt_hilbert, xl + len, yl + len, zl + len, 8);
      } else if ((x1 - xl) < len_nca && (y1 - yl) >= len_nca && (z1 - zl) >= len_nca) {
        index1 = findIndex(pt_hilbert, xl, yl + len, zl + len, 8);
      } else if ((x1 - xl) < len_nca && (y1 - yl) < len_nca && (z1 - zl) >= len_nca) {
        index1 = findIndex(pt_hilbert, xl, yl, zl + len, 8);
      }

      if ((x2 - xl) < len_nca && (y2 - yl) < len_nca && (z2 - zl) < len_nca) {
        index2 = findIndex(pt_hilbert, xl, yl, zl, 8);

      } else if ((x2 - xl) < len_nca && (y2 - yl) >= len_nca && (z2 - zl) < len_nca) {
        index2 = findIndex(pt_hilbert, xl, yl + len, zl, 8);
      } else if ((x2 - xl) >= len_nca && (y2 - yl) >= len_nca && (z2 - zl) < len_nca) {
        index2 = findIndex(pt_hilbert, xl + len, yl + len, zl, 8);
      } else if ((x2 - xl) >= len_nca && (y2 - yl) < len_nca && (z2 - zl) < len_nca) {
        index2 = findIndex(pt_hilbert, xl + len, yl, zl, 8);
      } else if ((x2 - xl) >= len_nca && (y2 - yl) < len_nca && (z2 - zl) >= len_nca) {
        index2 = findIndex(pt_hilbert, xl + len, yl, zl + len, 8);
      } else if ((x2 - xl) >= len_nca && (y2 - yl) >= len_nca && (z2 - zl) >= len_nca) {
        index2 = findIndex(pt_hilbert, xl + len, yl + len, zl + len, 8);
      } else if ((x2 - xl) < len_nca && (y2 - yl) >= len_nca && (z2 - zl) >= len_nca) {
        index2 = findIndex(pt_hilbert, xl, yl + len, zl + len, 8);
      } else if ((x2 - xl) < len_nca && (y2 - yl) < len_nca && (z2 - zl) >= len_nca) {
        index2 = findIndex(pt_hilbert, xl, yl, zl + len, 8);
      }


      if (index1 < index2) {
        return true;
      } else if (index1 > index2) {
        return false;
      }
      // means that index1 ==index2
      switch (index1) {
      case 0:
        pt_hilbert_new[0] = pt_hilbert[0];
        pt_hilbert_new[1] = (pt_hilbert[0] + pt_hilbert[7]) / 2;
        pt_hilbert_new[2] = (pt_hilbert[0] + pt_hilbert[4]) / 2;
        pt_hilbert_new[3] = (pt_hilbert[0] + pt_hilbert[3]) / 2;

        pt_hilbert_new[4] = (pt_hilbert[0] + pt_hilbert[2]) / 2;
        pt_hilbert_new[5] = (pt_hilbert[0] + pt_hilbert[5]) / 2;
        pt_hilbert_new[6] = (pt_hilbert[0] + pt_hilbert[6]) / 2;
        pt_hilbert_new[7] = (pt_hilbert[0] + pt_hilbert[1]) / 2;
        break;
      case 1:
        pt_hilbert_new[0] = (pt_hilbert[0] + pt_hilbert[1]) / 2;
        pt_hilbert_new[1] = pt_hilbert[1];
        pt_hilbert_new[2] = (pt_hilbert[1] + pt_hilbert[6]) / 2;
        pt_hilbert_new[3] = (pt_hilbert[1] + pt_hilbert[7]) / 2;

        pt_hilbert_new[4] = (pt_hilbert[1] + pt_hilbert[4]) / 2;
        pt_hilbert_new[5] = (pt_hilbert[1] + pt_hilbert[5]) / 2;
        pt_hilbert_new[6] = (pt_hilbert[1] + pt_hilbert[2]) / 2;
        pt_hilbert_new[7] = (pt_hilbert[1] + pt_hilbert[3]) / 2;
        break;
      case 2:
        pt_hilbert_new[0] = (pt_hilbert[2] + pt_hilbert[0]) / 2;
        pt_hilbert_new[1] = (pt_hilbert[2] + pt_hilbert[1]) / 2;
        pt_hilbert_new[2] = pt_hilbert[2];
        pt_hilbert_new[3] = (pt_hilbert[2] + pt_hilbert[3]) / 2;

        pt_hilbert_new[4] = (pt_hilbert[2] + pt_hilbert[4]) / 2;
        pt_hilbert_new[5] = (pt_hilbert[2] + pt_hilbert[5]) / 2;
        pt_hilbert_new[6] = (pt_hilbert[2] + pt_hilbert[6]) / 2;
        pt_hilbert_new[7] = (pt_hilbert[2] + pt_hilbert[7]) / 2;

        break;
      case 3:
        pt_hilbert_new[0] = (pt_hilbert[3] + pt_hilbert[6]) / 2;
        pt_hilbert_new[1] = (pt_hilbert[3] + pt_hilbert[1]) / 2;
        pt_hilbert_new[2] = (pt_hilbert[3] + pt_hilbert[0]) / 2;
        pt_hilbert_new[3] = (pt_hilbert[3] + pt_hilbert[7]) / 2;

        pt_hilbert_new[4] = (pt_hilbert[3] + pt_hilbert[4]) / 2;
        pt_hilbert_new[5] = pt_hilbert[3];
        pt_hilbert_new[6] = (pt_hilbert[3] + pt_hilbert[2]) / 2;
        pt_hilbert_new[7] = (pt_hilbert[3] + pt_hilbert[5]) / 2;
        break;

      case 4:
        pt_hilbert_new[0] = (pt_hilbert[4] + pt_hilbert[2]) / 2;
        pt_hilbert_new[1] = (pt_hilbert[4] + pt_hilbert[5]) / 2;
        pt_hilbert_new[2] = pt_hilbert[4];
        pt_hilbert_new[3] = (pt_hilbert[4] + pt_hilbert[3]) / 2;

        pt_hilbert_new[4] = (pt_hilbert[4] + pt_hilbert[0]) / 2;
        pt_hilbert_new[5] = (pt_hilbert[4] + pt_hilbert[7]) / 2;
        pt_hilbert_new[6] = (pt_hilbert[4] + pt_hilbert[6]) / 2;
        pt_hilbert_new[7] = (pt_hilbert[4] + pt_hilbert[1]) / 2;
        break;
      case 5:
        pt_hilbert_new[0] = (pt_hilbert[5] + pt_hilbert[0]) / 2;
        pt_hilbert_new[1] = (pt_hilbert[5] + pt_hilbert[1]) / 2;
        pt_hilbert_new[2] = (pt_hilbert[5] + pt_hilbert[2]) / 2;
        pt_hilbert_new[3] = (pt_hilbert[5] + pt_hilbert[3]) / 2;

        pt_hilbert_new[4] = (pt_hilbert[5] + pt_hilbert[4]) / 2;
        pt_hilbert_new[5] = pt_hilbert[5];
        pt_hilbert_new[6] = (pt_hilbert[5] + pt_hilbert[6]) / 2;
        pt_hilbert_new[7] = (pt_hilbert[5] + pt_hilbert[7]) / 2;
        break;
      case 6:
        pt_hilbert_new[0] = (pt_hilbert[6] + pt_hilbert[4]) / 2;
        pt_hilbert_new[1] = (pt_hilbert[6] + pt_hilbert[5]) / 2;
        pt_hilbert_new[2] = (pt_hilbert[6] + pt_hilbert[2]) / 2;
        pt_hilbert_new[3] = (pt_hilbert[6] + pt_hilbert[3]) / 2;

        pt_hilbert_new[4] = (pt_hilbert[6] + pt_hilbert[0]) / 2;
        pt_hilbert_new[5] = (pt_hilbert[6] + pt_hilbert[1]) / 2;
        pt_hilbert_new[6] = pt_hilbert[6];
        pt_hilbert_new[7] = (pt_hilbert[6] + pt_hilbert[7]) / 2;
        break;
      case 7:
        pt_hilbert_new[0] = (pt_hilbert[7] + pt_hilbert[6]) / 2;
        pt_hilbert_new[1] = (pt_hilbert[7] + pt_hilbert[1]) / 2;
        pt_hilbert_new[2] = (pt_hilbert[7] + pt_hilbert[2]) / 2;
        pt_hilbert_new[3] = (pt_hilbert[7] + pt_hilbert[5]) / 2;

        pt_hilbert_new[4] = (pt_hilbert[7] + pt_hilbert[4]) / 2;
        pt_hilbert_new[5] = (pt_hilbert[7] + pt_hilbert[3]) / 2;
        pt_hilbert_new[6] = (pt_hilbert[7] + pt_hilbert[0]) / 2;
        pt_hilbert_new[7] = pt_hilbert[7];
        break;
      default:
        std::cout << "Hilbert ordering error:Invalid nearest cell" << std::endl;
        break;


      }


      pt_hilbert[0] = pt_hilbert_new[0];
      pt_hilbert[1] = pt_hilbert_new[1];
      pt_hilbert[2] = pt_hilbert_new[2];
      pt_hilbert[3] = pt_hilbert_new[3];

      pt_hilbert[4] = pt_hilbert_new[4];
      pt_hilbert[5] = pt_hilbert_new[5];
      pt_hilbert[6] = pt_hilbert_new[6];
      pt_hilbert[7] = pt_hilbert_new[7];

      len = len / 2;
      deapth += 1;

    }
    return false;



  }

}


inline bool TreeNode::hilbert_order_NCA(const Point& p1, const Point& p2) const {


  int g_dim = this->m_uiDim;
  unsigned int x1 = p1.xint();
  unsigned int x2 = p2.xint();

  unsigned int y1 = p1.yint();
  unsigned int y2 = p2.yint();

  unsigned int z1 = p1.zint();
  unsigned int z2 = p2.zint();

  if (x1 == x2 && y1 == y2 && z1 == z2) {
    return false;
  }

  unsigned int maxDepth = this->m_uiMaxDepth;
  unsigned int maxDiff = (unsigned int)(std::max((std::max((x1 ^ x2), (y1 ^ y2))), (z1 ^ z2)));
  int dim = g_dim;
  unsigned int maxDiffBinLen = binOp::binLength(maxDiff);
  //Eliminate the last maxDiffBinLen bits.
  unsigned int ncaX = ((x1 >> maxDiffBinLen) << maxDiffBinLen);
  unsigned int ncaY = ((y1 >> maxDiffBinLen) << maxDiffBinLen);
  unsigned int ncaZ = ((z1 >> maxDiffBinLen) << maxDiffBinLen);
  unsigned int ncaLev = (maxDepth - maxDiffBinLen);



  unsigned int xl = 0;
  unsigned int yl = 0;
  unsigned int zl = 0;

  unsigned int len = 1 << maxDepth;
  int count = 0;
  unsigned int index1 = 0;
  unsigned int index2 = 0;


  if (g_dim == 2) {
    int rotation[4] = { 0, 1, 2, 3 };
    int rot_index[4] = { 0, 1, 2, 3 };
    while ((xl != ncaX || yl != ncaY || zl != ncaZ || count != ncaLev)) {
      len = len / 2;

      index1 = 0;
      if (ncaX >= (len + xl)) {
        index1 += 1;
        xl += len;
        if (ncaY < (len + yl)) index1 += 2;
      }
      if (ncaY >= (len + yl)) {index1 += 1; yl += len; }

      rotate(index1, rotation, rot_index, dim);


      count++;

    }

    len = len / 2;
    if ((x1 - ncaX) < len && (y1 - ncaY) < len) { // index 0
      index1 = rot_index[0];

    } else if ((x1 - ncaX) < len && (y1 - ncaY) >= len) { // index 1
      index1 = rot_index[1];
    } else if ((x1 - ncaX) >= len && (y1 - ncaY) >= len) { // index 2
      index1 = rot_index[2];
    } else if ((x1 - ncaX) >= len && (y1 - ncaY) < len) { // index 3
      index1 = rot_index[3];

    }


    if ((x2 - ncaX) < len && (y2 - ncaY) < len) { // index 0
      index2 = rot_index[0];

    } else if ((x2 - ncaX) < len && (y2 - ncaY) >= len) { // index 1
      index2 = rot_index[1];

    } else if ((x2 - ncaX) >= len && (y2 - ncaY) >= len) { // index 2
      index2 = rot_index[2];
    } else if ((x2 - ncaX) >= len && (y2 - ncaY) < len) { // index 3
      index2 = rot_index[3];

    }


  } else if (g_dim == 3) {
    int rotation[8] = { 0, 1, 2, 3, 4, 5, 6, 7 }; // Initial rotation
    int rot_index[8] = { 0, 1, 2, 3, 4, 5, 6, 7 }; // Initial rotation indices
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

      rotate(index1, rotation, rot_index, dim);
      count++;

    }

    len >>= 1;

    if ((x1 - ncaX) < len && (y1 - ncaY) < len && (z1 - ncaZ) < len) {
      index1 = rot_index[0];
    } else if ((x1 - ncaX) < len && (y1 - ncaY) >= len && (z1 - ncaZ) < len) {
      index1 = rot_index[1];
    } else if ((x1 - ncaX) >= len && (y1 - ncaY) >= len && (z1 - ncaZ) < len) {
      index1 = rot_index[2];
    } else if ((x1 - ncaX) >= len && (y1 - ncaY) < len && (z1 - ncaZ) < len) {
      index1 = rot_index[3];
    } else if ((x1 - ncaX) >= len && (y1 - ncaY) < len && (z1 - ncaZ) >= len) {
      index1 = rot_index[4];
    } else if ((x1 - ncaX) >= len && (y1 - ncaY) >= len && (z1 - ncaZ) >= len) {
      index1 = rot_index[5];
    } else if ((x1 - ncaX) < len && (y1 - ncaY) >= len && (z1 - ncaZ) >= len) {
      index1 = rot_index[6];
    } else if ((x1 - ncaX) < len && (y1 - ncaY) < len && (z1 - ncaZ) >= len) {
      index1 = rot_index[7];
    }


    if ((x2 - ncaX) < len && (y2 - ncaY) < len && (z2 - ncaZ) < len) {
      index2 = rot_index[0];
    } else if ((x2 - ncaX) < len && (y2 - ncaY) >= len && (z2 - ncaZ) < len) {
      index2 = rot_index[1];
    } else if ((x2 - ncaX) >= len && (y2 - ncaY) >= len && (z2 - ncaZ) < len) {
      index2 = rot_index[2];
    } else if ((x2 - ncaX) >= len && (y2 - ncaY) < len && (z2 - ncaZ) < len) {
      index2 = rot_index[3];
    } else if ((x2 - ncaX) >= len && (y2 - ncaY) < len && (z2 - ncaZ) >= len) {
      index2 = rot_index[4];
    } else if ((x2 - ncaX) >= len && (y2 - ncaY) >= len && (z2 - ncaZ) >= len) {
      index2 = rot_index[5];
    } else if ((x2 - ncaX) < len && (y2 - ncaY) >= len && (z2 - ncaZ) >= len) {
      index2 = rot_index[6];
    } else if ((x2 - ncaX) < len && (y2 - ncaY) < len && (z2 - ncaZ) >= len) {
      index2 = rot_index[7];
    }



  }

  return index1 < index2;

}


inline bool TreeNode::morton_order(const Point& p1, const Point& p2) const {


  unsigned int x1 = p1.xint();
  unsigned int x2 = p2.xint();

  unsigned int y1 = p1.yint();
  unsigned int y2 = p2.yint();

  unsigned int z1 = p1.zint();
  unsigned int z2 = p2.zint();

  if (x1 == x2 && y1 == y2 && z1 == z2) {
    return false;
  }

  unsigned int x = (x1 ^ x2);
  unsigned int y = (y1 ^ y2);
  unsigned int z = (z1 ^ z2);

  //Default pref: z > y > x.
  unsigned int maxC = z;
  unsigned int yOrx = y;
  if (yOrx < x) {if ((x ^ yOrx) >= yOrx) {yOrx = x;}
  }
  if (maxC < yOrx) {if ((maxC ^ yOrx) >= maxC) {maxC = yOrx;}
  }

  if (maxC == z) {
    if (z1 < z2) return true;
    else return false;

  } else if (maxC == y) {
    if (y1 < y2) return true;
    else return false;

  } else {
    if (x1 < x2) return true;
    else return false;

  }



}

inline bool TreeNode::morton_order_NCA(const Point& p1, const Point& p2) const {

  int g_dim = this->m_uiDim;

  unsigned int x1 = p1.xint();
  unsigned int x2 = p2.xint();

  unsigned int y1 = p1.yint();
  unsigned int y2 = p2.yint();

  unsigned int z1 = p1.zint();
  unsigned int z2 = p2.zint();

  if (x1 == x2 && y1 == y2 && z1 == z2) {
    return false;
  }

  unsigned int maxDepth = this->m_uiMaxDepth;
  unsigned int maxDiff = (unsigned int)(std::max((std::max((x1 ^ x2), (y1 ^ y2))), (z1 ^ z2)));

  unsigned int maxDiffBinLen = binOp::binLength(maxDiff);
  //Eliminate the last maxDiffBinLen bits.
  unsigned int ncaX = ((x1 >> maxDiffBinLen) << maxDiffBinLen);
  unsigned int ncaY = ((y1 >> maxDiffBinLen) << maxDiffBinLen);
  unsigned int ncaZ = ((z1 >> maxDiffBinLen) << maxDiffBinLen);
  unsigned int ncaLev = (maxDepth - maxDiffBinLen);

  unsigned int xl = 0;
  unsigned int yl = 0;
  unsigned int zl = 0;

  unsigned int len = 1 << maxDepth;

  len = len / (1 << (ncaLev + 1));
  unsigned int index1 = 0;
  unsigned int index2 = 0;

  if (g_dim == 2) {


   index1 = 0;
    if ( x1>=(len+ncaX) ) {
	index1 += 1;
      if (y1>=(len+ncaY)){ 
	  index1 += 2;
      }
    }else if ( y1>=(len+ncaY) ) { 
      index1 += 2;
      
    }

    index2 = 0;
    if ( x2>=(len+ncaX) ) {
	index2 += 1;
      if (y2>=(len+ncaY)){
	index2 += 2;
      }
    }else if ( y2>=(len+ncaY) ) {
      index2 += 2;
      
    }
   
   
  } else if (g_dim == 3) {
    
      index1=0;
      if ( z1 < (len + ncaZ) ) {
	 if ( x1 >= (len + ncaX) ) {
	    index1 += 1;
	    if (y1 >= (len + ncaY)) 
		index1 += 2;
	 }else if ( y1>= (len + ncaY) ) { 
	    index1 += 2; 
	 }
      } else {
	    index1 = 4;
	 if ( x1 >= (len + ncaX) ) {
	    index1 += 1;
	    if (y1 >= (len + ncaY)) 
	      index1 += 2;
	 }else if ( y1>= (len + ncaY) ) { 
	    index1 += 2; 
	 }
      }
      
      index2=0;
      if ( z2 < (len + ncaZ) ) {
	 if ( x2 >= (len + ncaX) ) {
	    index2 += 1;
	    if (y2 >= (len + ncaY)) 
		index2 += 2;
	 }else if ( y2>= (len + ncaY) ) { 
	    index2 += 2; 
	 }
      } else {
	    index2 = 4;
	 if ( x2 >= (len + ncaX) ) {
	    index2 += 1;
	    if (y2 >= (len + ncaY)) 
	      index2 += 2;
	 }else if ( y2>= (len + ncaY) ) { 
	    index2 += 2; 
	 }
      }

  }
  return index1 < index2;


}


//}  Hilbert  related changes ===================

} //end namespace ot

