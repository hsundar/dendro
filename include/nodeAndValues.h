
/**
  @file nodeAndValues.h
  @brief A small helper class that pairs TreeNode with an array of values.
  The class is templated on the length and type of the array. This class has an 
  MPI_Datatype associated with it and it can be communicated between processors. 
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  */

#ifndef _NODE_AND_VALUES_H_
#define _NODE_AND_VALUES_H_

#include "TreeNode.h"

namespace ot {

  /**
    @brief A small helper class that pairs an octant with some values. 
    The class is templated on the length and type of the array. This class has an 
    MPI_Datatype associated with it and it can be communicated between processors. 
    @author Rahul Sampath
    */
  template<typename T, unsigned int ARR_LEN>
    class NodeAndValues {

    public:
        ot::TreeNode node; /**< The octant */
        T values[ARR_LEN]; /**< The values */

        /** @name Overload operators */
        //@{
        bool Equals(NodeAndValues<T, ARR_LEN> const & other) const {
          bool result = ((this->node) == other.node);

          for(unsigned int i = 0; (i < ARR_LEN) && result; i++) {
            result = ( (this->values[i]) == (other.values[i]) );
          }

          return result;
        }

        bool operator == ( NodeAndValues<T, ARR_LEN>  const & other)  const{
          return ((this->node) == other.node);
        }//end fn.

        bool operator != ( NodeAndValues<T, ARR_LEN> const  & other)  const{
          return ((this->node) != other.node);
        }//end fn.

        bool operator  < ( NodeAndValues<T, ARR_LEN>  const & other)  const{
          return ((this->node) < other.node);
        }//end function

        bool operator  <= ( NodeAndValues<T, ARR_LEN>  const  & other)  const{
          return ((this->node) <= other.node);
        }//end fn.

        bool operator  > ( NodeAndValues<T, ARR_LEN> const  & other)  const{
          return ((this->node) > other.node);
        }//end fn.

        bool operator  >= ( NodeAndValues<T, ARR_LEN>  const & other)  const{    
          return ((this->node) >= other.node);
        }//end fn.

        //Asignment Operator
       NodeAndValues<T, ARR_LEN>& operator = ( NodeAndValues<T, ARR_LEN> const & other) {
          if(this == (&other)) {return *this;}	
          this->node = other.node;
          for(unsigned int i = 0; i < ARR_LEN; i++) {
            this->values[i] = other.values[i];
          }
          return *this;
        }//end fn.



        friend std::ostream& operator<<(std::ostream& os, NodeAndValues<T, ARR_LEN> const& other) {
	  return (os << other.node.getX() << " " << other.node.getY() << " " << other.node.getZ() << " " << other.node.getLevel());
	} //end fn.
        
        
        //@}

        /** @name Constructors */
        //@{
        NodeAndValues () {  }

        //copy constructor	
        NodeAndValues  (const NodeAndValues<T, ARR_LEN> & other) {
          this->node = other.node;
          for(unsigned int i = 0; i < ARR_LEN; i++) {
            this->values[i] = other.values[i];
          }
        }
        //@}

    };//end class definition

}//end namespace

namespace par {

  //Forward Declaration
  template <typename T>
    class Mpi_datatype;

  /**
    @author Rahul Sampath, rahul.sampath@gmail.com
    @brief A template specialization of the abstract class "Mpi_datatype" for
    communicating messages of type "ot::NodeAndValues<T, unsigned int>".
    */
  template <typename T, unsigned int ARR_LEN>
    class Mpi_datatype< ot::NodeAndValues<T, ARR_LEN> > {

      public:

        /**
          @return The MPI_Datatype corresponding to the datatype "ot::NodeAndValues".
          */
        static MPI_Datatype value()
        {
          static bool         first = true;
          static MPI_Datatype datatype;

          if (first)
          {
            first = false;
            MPI_Type_contiguous(sizeof(ot::NodeAndValues<T, ARR_LEN>), MPI_BYTE, &datatype);
            MPI_Type_commit(&datatype);
          }

          return datatype;
        }

    };

}//end namespace par


#endif


