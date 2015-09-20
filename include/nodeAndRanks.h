
/**
  @file nodeAndRanks.h
  @brief A small helper class that pairs TreeNode with a list of processors.
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  */

#ifndef _NODE_AND_RANKS_H_
#define _NODE_AND_RANKS_H_

#include <vector>
#include "TreeNode.h"

namespace ot {

  /**
    @brief A small helper class that pairs an octant with a list of processors
    @author Rahul Sampath
    */
  class NodeAndRanks {
    public:

      ot::TreeNode node; /**< The octant */
      std::vector<int> ranks; /**< The ranks of the processors */

      /** @name Overload operators */
      //@{
      bool operator == ( NodeAndRanks  const & other)  const{
        return ((this->node) == other.node);
      }//end fn.

      bool operator != ( NodeAndRanks const  & other)  const{
        return ((this->node) != other.node);
      }//end fn.

      bool operator  < ( NodeAndRanks  const & other)  const{
        return ((this->node) < other.node);
      }//end function

      bool operator  <= ( NodeAndRanks  const  & other)  const{
        return ((this->node) <= other.node);
      }//end fn.

      bool operator  > ( NodeAndRanks const  & other)  const{
        return ((this->node) > other.node);
      }//end fn.

      bool operator  >= ( NodeAndRanks  const & other)  const{    
        return ((this->node) >= other.node);
      }//end fn.

      //Asignment Operator
      NodeAndRanks & operator = ( NodeAndRanks const & other) {
        if(this == (&other)) {return *this;}	
        this->node = other.node;
        this->ranks = other.ranks;    
        return *this;
      }//end fn.

      //@}

      /** @name Constructors */
      //@{
      NodeAndRanks () {      }

      //copy constructor	
      NodeAndRanks  (const NodeAndRanks  &other) {
        node = other.node;
        ranks = other.ranks;    
      }
      //@}

  };

}//end namespace


#endif


