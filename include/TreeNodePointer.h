
/**
 * @file TreeNodePointer.h
 * @author	Rahul S. Sampath, rahul.sampath@gmail.com
 * 
**/ 

#ifndef _TREENODE_POINTER_H_
#define _TREENODE_POINTER_H_

namespace ot {

  class TreeNode;
  
  /**
    @brief A wrapper for linked list based octree representation 
    @author Rahul Sampath
    */
  struct TreeNodePointer {
      TreeNode m_tnMe;
      TreeNodePointer* m_tnpMyChildren;
  };
 
}//end namespace 


#endif /*TREENODE_POINTER_H_*/

