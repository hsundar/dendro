
#ifndef __MAT_RECORD_H__
#define __MAT_RECORD_H__

#include "petsc.h"

namespace ot {

/**
    @brief A class to be used while setting values into a Matrix.
    @author Rahul Sampath , rahul.sampath@gmail.com
    */
  class MatRecord {
    public:

      unsigned int rowIdx; /**< The local index of the node corresponding to the row
                             of the matrix. This can be got using the function ot::DA::getNodeIndices() */

      unsigned int colIdx; /**< The local index of the node corresponding to the column 
                             of the matrix. This can be got using the function ot::DA::getNodeIndices(). */

      unsigned int rowDim; /**< The degree of freedom (0-based) of the node
                             corresponding to  the row of the matrix. */

      unsigned int colDim; /**< The degree of freedom (0-based) of the node 
                             corresponding to  the column of the matrix. */

      PetscScalar val; /**< The matrix entry. */

      /** @name constructors */
      //@{

      /**
        The default constructor
        */
      MatRecord() {
        rowIdx = colIdx = rowDim = colDim = static_cast<unsigned int>(-1);
        val = 0;
      }    

      /**
        The copy constructor          
        */
      MatRecord(const MatRecord & other) {
        rowIdx = other.rowIdx;
        colIdx = other.colIdx;
        rowDim = other.rowDim;
        colDim = other.colDim;
        val = other.val;
      }
      //@}

      /**
        The assignment operator
        */
      MatRecord & operator = (MatRecord const  & other) {
        if(this == (&other)) {return *this;}	
        rowIdx = other.rowIdx;
        colIdx = other.colIdx;
        rowDim = other.rowDim;
        colDim = other.colDim;
        val = other.val;
        return *this;
      }//end fn.

      /** @name Overloaded Operators */
      //@{
      bool  operator == ( MatRecord const  &other) const {
        return ( (rowIdx == other.rowIdx) && (colIdx == other.colIdx)
            && (val == other.val) && (rowDim == other.rowDim) && (colDim == other.colDim) );
      }

      bool  operator != (MatRecord const  &other) const {
        return (!((*this) == other));
      }

      bool  operator < (MatRecord const  &other) const {
        if(rowIdx < other.rowIdx) {
          return true;
        }else if(rowIdx == other.rowIdx) {
          if(rowDim < other.rowDim) {
            return true;
          }else if(rowDim == other.rowDim) {
            if(colIdx < other.colIdx) {
              return true;
            }else if (colIdx == other.colIdx) {
              if(colDim < other.colDim) {
                return true;
              }else if(colDim == other.colDim) {
                return (val < other.val);
              }else {
                return false;
              }
            }else {
              return false;
            }
          }else {
            return false;
          }
        }else {
          return false;
        }
      }

      bool  operator > (MatRecord const  &other) const {
        return ( (!((*this) < other)) && ((*this) != other) );
      }

      bool  operator <= (MatRecord const  &other) const {
        return ( ((*this) < other) || ((*this) == other) );
      }

      bool  operator >= (MatRecord const  &other) const {
        return (!((*this) < other)) ;
      }
      //@}
  };

} //end namespace

#endif

