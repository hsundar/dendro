
/**
 * @file TreeNode.h
 * @author	Rahul S. Sampath, rahul.sampath@gmail.com
 @author Hari Sundar, hsundar@gmail.com
 * 
**/ 

#ifndef _TREENODE_H_
#define _TREENODE_H_

#include "octUtils.h"
#include "Point.h"
#include <iostream>

#ifdef __DEBUG__
#ifndef __DEBUG_TN__
#define __DEBUG_TN__
#endif
#endif

namespace ot {

  /**
    @brief A class to manage octants.
    @author Rahul Sampath
    */
  class TreeNode {
    protected:
      //Level is also used as a flag.
      unsigned int m_uiX, m_uiY, m_uiZ, m_uiLevel, m_uiWeight, m_uiDim, m_uiMaxDepth;      
      //constructor without checks: only for faster construction.
      TreeNode (const int dummy, const unsigned int x,const unsigned int y,
          const unsigned int z,const unsigned int level, const unsigned int dim,
          const unsigned int maxDepth);

    public:

      /**
        @brief The type of boundary
        */
      enum BoundaryType1 { NEGATIVE= 2, POSITIVE= 4};

      /**
        @brief The type of boundary
        */
      enum BoundaryType2 {
        X_NEG_BDY=1, Y_NEG_BDY=2, Z_NEG_BDY=4, NEG_POS_DEMARCATION=8, 
        EXTERNAL_BDY=16, X_POS_BDY=32, Y_POS_BDY=64, Z_POS_BDY=128
      };

      /**
        @brief  The type of boundary
        */
      enum BoundaryType3 {
        FACE_BDY=1, EDGE_BDY=2, CORNER_BDY=3
      };

      enum OctantFlagType {
       MAX_LEVEL=31, BOUNDARY=64, NODE=128 
      };

      /**
        @author Rahul Sampath
        @brief returns the list of decendants of this octant at the
        same level as other's parent and also touch 'other'. 
        @param other an octant at a level > this octant's level + 1
        @param seeds the output
        @param incCorners 'true' for balancing across corners as well and 'false' otherwise.  
        @return error flag

        other should not be a decendant of this octant.
        */
      int addBalancingDescendants(TreeNode other, std::vector<TreeNode> &seeds,
          bool incCorners) const;

      /**
        @author Rahul Sampath
        @brief balances the subtree formed with this octant as its root.
        @param in input vector containing the decendants of this octant
        @param out output vector  
        @param incCorners 'true' for balancing across corners as well and 'false' otherwise.  
        @param isSortedCompleteLinear 'true' if the input is sorted, complete and linear 
        */
      int balanceSubtree(std::vector<TreeNode > & in, std::vector<TreeNode > & out,
          bool incCorners, bool isSortedCompleteLinear) const;      

      /**
        @author Rahul Sampath
        @brief generates the coarsest possible complete, linear
        subtree containing the octants in 'in'
        @param in a set of decendants of this octant
        @param out the complete subtree of this octant
        */
      int completeSubtree(std::vector<TreeNode > & in, std::vector<TreeNode > & out) const;

      /**
        @author Rahul Sampath
        @brief returns the subset of the input that touch this
        octant's boundary and are this octant's decendants.
        @param in a set of decendants of this octant
        @param out the output vector
        */
      int pickInternalBoundaryCells(std::vector<TreeNode >& in, std::vector<TreeNode > & out) const;      

      /**
        @author Rahul Sampath
        @brief separates the input into two lists: 
        one containing this octant's decendants that touch this octant's boundary and
        another containing the rest of the octants in the input
        @param allInternal a set of this octant's decendants
        @param onlyInternal decendants of this octant that lie in the interior of this octant
        @param boundary decendants of this octant that lie on its boundary
        */
      int splitInternalAndBoundaryCells(std::vector<TreeNode> & allInternal,
          std::vector<TreeNode> & onlyInternal,std::vector<TreeNode> & boundary) const;

      /** @name Constructors */
      //@{

      /**
        @author Rahul Sampath
        @brief Constructs a root octant
        @param dim the dimension of the tree
        @param maxDepth The maximum depth of the tree. Must be <= MAX_LEVEL
        @see MAX_LEVEL
        */
      TreeNode (const unsigned int dim, const unsigned int maxDepth);

      /**
        @author Rahul Sampath
        @brief Constructs an octant
        @param x,y,z The coordinates of the anchor of the octant
        @param level The level of the octant
        @param dim the dimension of the tree
        @param maxDepth The maximum depth of the tree. Must be <= MAX_LEVEL
        @see MAX_LEVEL
        */
      TreeNode (const unsigned int x,const unsigned int y,const unsigned int z,
          const unsigned int level, const unsigned int dim, const unsigned int maxDepth);

      /**
        @author Rahul Sampath
        @brief Copy constructor
        */
      TreeNode (const TreeNode & other);

      /**
        @author Rahul Sampath
        @brief Default constructor. Constructs an invalid octant
        */
      TreeNode ();
      //@}

      /** @name Overloaded Operators */
      //@{

      /**
        @author Rahul Sampath
        @brief Assignment operator. No checks for dim or maxD are performed.
        It's ok to change dim and maxD of the current octant using the assignment operator.
        */
      TreeNode & operator = (TreeNode const  & other);

      /**
        @author Rahul Sampath
        @brief Two octants are equal if their respective anchors are equal and their levels are
        equal.  
        */
      bool  operator == ( TreeNode const  &other) const;

      /**
        @author Rahul Sampath
        @brief Two octants are equal if their respective anchors are equal and their levels are
        equal.  
        */
      bool  operator != (TreeNode const  &other) const;

      /**
        @author Rahul Sampath
        @author Hari Sundar
        @brief The comparisons are based on the Morton ordering of the octants 
        */
      bool  operator < ( TreeNode const  &other) const;

      /**
        @author Rahul Sampath
        @brief The comparisons are based on the Morton ordering of the octants 
        */
      bool  operator > ( TreeNode const  &other) const;

      /**
        @author Rahul Sampath
        @brief The comparisons are based on the Morton ordering of the octants 
        */
      bool  operator <= ( TreeNode const  &other) const;

      /**
        @author Rahul Sampath
        @brief The comparisons are based on the Morton ordering of the octants 
        */
      bool  operator >= ( TreeNode const  &other) const;

      /**
        @author Rahul Sampath
        @brief Appends the anchor and level of the octant to the stream
        */
      friend std::ostream & operator << (std::ostream & os,TreeNode const & node) ;  
      //@}

      /**
        @author Rahul Sampath
        @name Getters and Setters
        */
      //@{
      unsigned int getDim() const;
      unsigned int getMaxDepth() const;
      unsigned int getWeight() const;
      unsigned int getLevel() const;
      unsigned int getFlag() const;
      unsigned int getX() const;
      unsigned int getY() const;
      unsigned int getZ() const;
      int getAnchor(unsigned int &x, unsigned int&y, unsigned int&z) const;
      Point getAnchor() const { return Point(m_uiX, m_uiY, m_uiZ); };
      unsigned int getParentX() const;
      unsigned int getParentY() const;
      unsigned int getParentZ() const;    
      unsigned char getChildNumber() const;
      int setWeight(unsigned int w);
      int addWeight(unsigned int w);
      int setFlag(unsigned int w);
      int orFlag(unsigned int w);
      //@}

      /**
        @author Rahul Sampath
        @brief Checks if this octant is an ancestor of the input
        @return 'true' if this is an ancestor of 'other'
        */
      bool isAncestor(const TreeNode & other) const;

      /**
        @author Hari Sundar
        @brief flags is a datastructure which will store which boundaries were
        touched. highest 3 bits are for +Z,+y, and +x axes ... and and
        smallest 3 are for -z,-y and -x axes.
        */
      bool isBoundaryOctant(int type=POSITIVE, unsigned char *flags=NULL) const;

      /**
        @author Hari Sundar
        @brief flags is a datastructure which will store which boundaries were
        touched. highest 3 bits are for +Z,+y, and +x axes ... and and
        smallest 3 are for -z,-y and -x axes.
        */
      bool isBoundaryOctant(const TreeNode &block, int type=POSITIVE, unsigned char *flags=NULL) const;

      /**
        @author Rahul Sampath
        @brief appends the children (sorted) of this octant to 'list'
        */
      int addChildren(std::vector<TreeNode > &list) const;

      /**
        @author Rahul Sampath
        @brief appends the siblings (sorted) of this octant to 'brothers'
        */
      int addBrothers(std::vector<TreeNode>&brothers) const;

      /**
        @author Rahul Sampath
        @return the parent of this octant
        */
      TreeNode  getParent() const;

      /**
        @author Rahul Sampath
        @param The level of the ancestor
        @return the ancestor of this octant at level 'ancLev'        
        */
      TreeNode	getAncestor(unsigned int ancLev) const;

      /**
        @author Rahul Sampath
        @return the Deepest first decendant of this octant
        */
      TreeNode getDFD() const;

      /**
        @author Rahul Sampath
        @return the deepest last decendant of this octant
        */
      TreeNode getDLD() const;

      /**
       *@author Dhairya Malhotra
       *@return the next octant in morton order which has the least depth.
       */
      TreeNode getNext() const;

      /**
       *@author Dhairya Malhotra
       *@return the smallest (in morton id) of this octant
       */
      TreeNode getFirstChild() const;

      /**
        @author Rahul Sampath
        @return the Morton encoding for this octant
        */
      std::vector<bool> getMorton() const;

      unsigned int minX() const;
      unsigned int minY() const;
      unsigned int minZ() const;
      unsigned int maxX() const;
      unsigned int maxY() const;
      unsigned int maxZ() const;

      std::vector<TreeNode> getSearchKeys(bool incCorners);

      /**
        @author Rahul Sampath
        @name Get Neighbours at the same level as the current octant
        */
      //@{
      TreeNode  getLeft() const;
      TreeNode  getRight() const;
      TreeNode  getTop() const;
      TreeNode  getBottom() const;
      TreeNode  getFront() const;
      TreeNode  getBack() const; 
      TreeNode  getTopLeft() const;
      TreeNode  getTopRight() const;
      TreeNode  getBottomLeft() const; 
      TreeNode  getBottomRight() const;
      TreeNode  getLeftFront() const;
      TreeNode  getRightFront() const;
      TreeNode  getTopFront() const;
      TreeNode  getBottomFront() const;
      TreeNode  getTopLeftFront() const;
      TreeNode  getTopRightFront() const;
      TreeNode  getBottomLeftFront() const;
      TreeNode  getBottomRightFront() const;
      TreeNode  getLeftBack() const;
      TreeNode  getRightBack() const;
      TreeNode  getTopBack() const;
      TreeNode  getBottomBack() const;
      TreeNode  getTopLeftBack() const;
      TreeNode  getTopRightBack() const;
      TreeNode  getBottomLeftBack() const;
      TreeNode  getBottomRightBack() const;
      std::vector<TreeNode> getAllNeighbours() const;
      //@}

      /**
        @author Rahul Sampath
        @name Get the list of potential neighbours that will satisfy the 2:1 balance constraint
        */
      //@{
      std::vector<TreeNode > getB_Left() const;
      std::vector<TreeNode > getB_Right() const;
      std::vector<TreeNode > getB_Top() const;
      std::vector<TreeNode > getB_Bottom() const;
      std::vector<TreeNode > getB_Front() const;
      std::vector<TreeNode > getB_Back() const; 
      std::vector<TreeNode > getB_TopLeft() const;
      std::vector<TreeNode > getB_TopRight() const;
      std::vector<TreeNode > getB_BottomLeft() const; 
      std::vector<TreeNode > getB_BottomRight() const;
      std::vector<TreeNode > getB_LeftFront() const;
      std::vector<TreeNode > getB_RightFront() const;
      std::vector<TreeNode > getB_TopFront() const;
      std::vector<TreeNode > getB_BottomFront() const;
      std::vector<TreeNode > getB_TopLeftFront() const;
      std::vector<TreeNode > getB_TopRightFront() const;
      std::vector<TreeNode > getB_BottomLeftFront() const;
      std::vector<TreeNode > getB_BottomRightFront() const;
      std::vector<TreeNode > getB_LeftBack() const;
      std::vector<TreeNode > getB_RightBack() const;
      std::vector<TreeNode > getB_TopBack() const;
      std::vector<TreeNode > getB_BottomBack() const;
      std::vector<TreeNode > getB_TopLeftBack() const;
      std::vector<TreeNode > getB_TopRightBack() const;
      std::vector<TreeNode > getB_BottomLeftBack() const;
      std::vector<TreeNode > getB_BottomRightBack() const;
      std::vector<std::vector<TreeNode> > getAllB_Neighbours() const;
      //@}

  };//end of the class definition.

}//end namespace

#include "TreeNode.txx"

namespace par {

  //Forward Declaration
  template <typename T>
    class Mpi_datatype;

/**
@author Rahul Sampath, rahul.sampath@gmail.com
@brief A template specialization of the abstract class "Mpi_datatype" for communicating messages of type "ot::TreeNode".
*/
  template <>
    class Mpi_datatype< ot::TreeNode > {

      static void Node_MAX_LEVEL(void *in, void *inout, int* len, MPI_Datatype * dptr) {
        for(int i = 0; i < (*len); i++) {
          ot::TreeNode first = (static_cast<ot::TreeNode*>(in))[i];
          ot::TreeNode second = (static_cast<ot::TreeNode*>(inout))[i];
          (static_cast<ot::TreeNode*>(inout))[i] =
            ( ( (first.getLevel()) > (second.getLevel()) )? first : second );
        }//end for	
      }//end function

      static void Node_MAX(void *in, void *inout, int* len, MPI_Datatype * dptr) {
        for(int i = 0; i < (*len); i++) {
          ot::TreeNode first = (static_cast<ot::TreeNode*>(in))[i];
          ot::TreeNode second = (static_cast<ot::TreeNode*>(inout))[i];
          (static_cast<ot::TreeNode*>(inout))[i] = ( ( first > second )? first : second );
        }//end for	
      }//end function

      static void Node_MIN(void *in, void *inout, int* len, MPI_Datatype * dptr) {
        for(int i = 0; i < (*len); i++) {
          ot::TreeNode first = (static_cast<ot::TreeNode*>(in))[i];
          ot::TreeNode second = (static_cast<ot::TreeNode*>(inout))[i];
          (static_cast<ot::TreeNode*>(inout))[i] = ( ( first < second )? first : second );
        }//end for	
      }//end function

      static void Node_NCA(void *in, void *inout, int* len, MPI_Datatype * dptr) {
        for(int i = 0; i < (*len); i++) {
          ot::TreeNode first = (static_cast<ot::TreeNode*>(in))[i];
          ot::TreeNode second = (static_cast<ot::TreeNode*>(inout))[i];
          if( first != second ) {
            (static_cast<ot::TreeNode*>(inout))[i] = ot::getNCA(first, second);
          }//end if
        }//end for	
      }//end function

      public:
      /** 
        @brief User defined MPI_Operation that sets second[i] to first[i] if first[i] is at a greater level than second[i]. 
        @remark first and second are 2 arrays of type TreeNode. 
      **/
      static MPI_Op MAX_LEVEL() {
        static bool         first = true;
        static MPI_Op maxLev;
        if (first) {
          first = false;
          MPI_Op_create(Mpi_datatype<ot::TreeNode>::Node_MAX_LEVEL ,true ,&maxLev);
        }
        return maxLev;
      }

      /** 
        @brief User defined MPI_Operation that computes: second[i] = Max(first[i], second[i]), 
        @remark first and second are 2 arrays of type TreeNode. 
        "MAX" is a macro
      **/
      static MPI_Op _MAX() {
        static bool         first = true;
        static MPI_Op max;
        if (first) {
          first = false;
          MPI_Op_create(Mpi_datatype<ot::TreeNode>::Node_MAX ,true ,&max);
        }
        return max;
      }

      /** 
        @brief User defined MPI_Operation that computes: second[i] = Min(first[i], second[i]), 
        @remark first and second are 2 arrays of type TreeNode. 
        "MIN" is a macro
      **/
      static MPI_Op _MIN() {
        static bool         first = true;
        static MPI_Op min;
        if (first) {
          first = false;
          MPI_Op_create(Mpi_datatype<ot::TreeNode>::Node_MIN ,true ,&min);
        }
        return min;
      }

      /** 
        @brief User defined MPI_Operation that computes: second[i] = NCA(first[i], second[i]), 
        @remark first and second are 2 arrays of type TreeNode and 
        NCA returns the nearest common ancestor of its 2 arguments. 
      **/
      static MPI_Op NCA() {
        static bool         first = true;
        static MPI_Op nca;
        if (first) {
          first = false;
          MPI_Op_create(Mpi_datatype<ot::TreeNode>::Node_NCA ,true ,&nca);
        }
        return nca;
      }

      /**
       @return The MPI_Datatype corresponding to the datatype "ot::TreeNode".
     */
      static MPI_Datatype value()
      {
        static bool         first = true;
        static MPI_Datatype datatype;

        if (first)
        {
          first = false;
          MPI_Type_contiguous(sizeof(ot::TreeNode), MPI_BYTE, &datatype);
          MPI_Type_commit(&datatype);
        }

        return datatype;
      }

    };

}//end namespace par

#endif /*NODE_H_*/

