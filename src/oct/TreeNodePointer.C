
/**
  @file TreeNodePointer.C
  @brief A collection of simple functions for manipulating pointer based octrees.  
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  */

#include "TreeNode.h"
#include "TreeNodePointer.h"

namespace ot {

  void appendOctantsAtLevel(const ot::TreeNodePointer & ptrOctree, 
      std::vector<ot::TreeNode> & wList, unsigned int lev) {
    if(ptrOctree.m_tnMe.getLevel() == lev) {
      wList.push_back(ptrOctree.m_tnMe);
    } else if(ptrOctree.m_tnMe.getLevel() < lev) {
      if(ptrOctree.m_tnpMyChildren) {
        for(int i = 0; i < 8; i++) {
          appendOctantsAtLevel(ptrOctree.m_tnpMyChildren[i], wList, lev);
        }
      }
    }
  }

  void convertLinearToPointer(const std::vector<ot::TreeNode> & linOct,
      ot::TreeNodePointer & ptrOct) {

    assert(!(linOct.empty()));

    unsigned int dim = linOct[0].getDim();
    unsigned int maxDepth = linOct[0].getMaxDepth();

    //Add root first
    ptrOct.m_tnMe = ot::TreeNode(dim, maxDepth);
    ptrOct.m_tnpMyChildren = NULL;

    for(unsigned int i = 0; i < linOct.size(); i++) {
      addOctantToTreeNodePointer(ptrOct, linOct[i]);
    }
  }

  //octant must be a decendant of ptrOct.m_tnMe 
  void addOctantToTreeNodePointer(ot::TreeNodePointer & ptrOct, const ot::TreeNode & octant) {
    if(ptrOct.m_tnpMyChildren == NULL) {
      ptrOct.m_tnpMyChildren = new ot::TreeNodePointer[8];
      std::vector<ot::TreeNode> tmpOcts;
      ptrOct.m_tnMe.addChildren(tmpOcts);
      for(unsigned int i = 0; i < 8; i++) {
        ptrOct.m_tnpMyChildren[i].m_tnMe = tmpOcts[i];
        ptrOct.m_tnpMyChildren[i].m_tnpMyChildren = NULL;
      }
    }
    for(unsigned int i = 0; i < 8; i++) {
      if(ptrOct.m_tnpMyChildren[i].m_tnMe.isAncestor(octant)) {
        addOctantToTreeNodePointer(ptrOct.m_tnpMyChildren[i], octant);
        break;
      }
    }
  }

  void findOctantOrFinestAncestor(ot::TreeNodePointer & octree,
      const ot::TreeNode & key, ot::TreeNodePointer* & result) {
    if(octree.m_tnMe.isAncestor(key)) {
      if(octree.m_tnpMyChildren) {
        for(unsigned int i = 0; i < 8; i++) {
          if(octree.m_tnpMyChildren[i].m_tnMe.isAncestor(key)) {
            findOctantOrFinestAncestor(octree.m_tnpMyChildren[i], key, result);
            break;
          } else if (octree.m_tnpMyChildren[i].m_tnMe == key) {
            result = &(octree.m_tnpMyChildren[i]);
            break;
          }
        }
      } else {
        result = &octree;
      }
    } else if(octree.m_tnMe == key) {
      result = &octree;
    } else {
      assert(false);
    }
  }

  //Pre-order traversal. So the resulting linear octree will be Morton sorted.
  void convertPointerToLinear(std::vector<ot::TreeNode> & linOct,
      const ot::TreeNodePointer & ptrOct) {
    if(ptrOct.m_tnpMyChildren) {
      for(unsigned int i = 0; i < 8; i++) {
        convertPointerToLinear(linOct, ptrOct.m_tnpMyChildren[i]);
      }
    } else {
      linOct.push_back(ptrOct.m_tnMe);
    }
  }

  void deleteTreeNodePointer(ot::TreeNodePointer & ptrOct) {
    if(ptrOct.m_tnpMyChildren) {
      for(unsigned int i = 0; i < 8; i++) {
        deleteTreeNodePointer((ptrOct.m_tnpMyChildren)[i]);
      }
      delete [] (ptrOct.m_tnpMyChildren);
      ptrOct.m_tnpMyChildren = NULL;
    }
  }

}//end namespace


