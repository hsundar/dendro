
/**
  @file testUtils.h
  @brief  A Set of utilities to test octrees.	
  @author	Rahul S. Sampath, rahul.sampath@gmail.com  
  */ 

#ifndef _TESTUTILS_H_
#define _TESTUTILS_H_

#include <mpi.h>
#include <vector>

namespace seq {

  /**
    @namespace test
    @author Rahul Sampath
    @brief A collection of functions for debugging
    */
  namespace test {

    /**
      @fn
      @param nodes[in] The vector of nodes that have to be tested.
      @return true if it is sorted and false otherwise
    **/ 
    template<typename T>
      bool isSorted(const std::vector<T >& nodes);

    template<typename T>
      bool isSorted(T * nodes, unsigned int sz);

    /**
      @fn
      @param nodes[in] The vector of nodes that have to be tested.
      @return true if it is sorted and unique and false otherwise
    **/ 
    template<typename T>
      bool isUniqueAndSorted(const std::vector<T >& nodes);

    template<typename T>
      bool isUniqueAndSorted(T * nodes, unsigned int sz);

  }//end namespace
}//end namespace

namespace par {

  /**
    @namespace test
    @author Rahul Sampath
    @brief A collection of functions for debugging
    */
  namespace test {

    template<typename T>
      bool isSorted(const std::vector<T>& nodes, MPI_Comm comm);

    template<typename T>
      bool isUniqueAndSorted(const std::vector<T >& nodes,MPI_Comm comm) ;

  }//end namespace
}//end namespace

namespace ot { 

  class TreeNode;

  namespace test {
    /**
      @fn
      @param nodes[in] The vector of nodes that have to be tested.
      @return false if any ancestor of any element in the vector is also present in the vector and true otherwise.
      @remark The function uses Binary Search and hence expects the vector of nodes to be sorted.
    **/ 
    bool isLinear(const std::vector<TreeNode >& nodes) ;

    /**
      @fn
      @param nodes[in] The vector of nodes that have to be tested.
      @return true if the sum of the volumes of the nodes in the vector equals the volume of the domain and false otherwise.
      @remark The vector of nodes is expected to be unique and linear.
    **/ 
    bool isComplete(const std::vector<TreeNode >& nodes) ;

    /**
      @fn
      @param dim[in] The dimension of the tree.
      @param maxDepth[in] The maximum depth of the tree.
      @param nodes[in] The vector of nodes that have to be tested.
      @return true if for every element in the vector, all the neighbours that balanced it are also present in the vector and false otherwise.
      @remark The function uses Binary Search and hence expects the vector of nodes to be sorted.
    **/ 
    bool isBalanced(unsigned int dim, unsigned int maxDepth, char* failFileName,
        const std::vector<TreeNode >& nodes, bool incCorn, unsigned int maxLevDiff = 1) ;

    bool isBalancedInternal(unsigned int dim, unsigned int maxDepth,
        char*failFileName,	const std::vector<TreeNode > & nodes,
        TreeNode holder, bool incCorn, unsigned int maxLevDiff = 1) ;

  }//end namespace
}//end namespace

#include "testUtils.tcc"

#endif

