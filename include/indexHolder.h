
/**
  @file indexHolder.h
  @brief A small helper class used while sorting a pair of value and its index in an array/vector. 
  @author Rahul S. Sampath, rahul.sampath@gmail.com

  This is useful when you want to map the sorted indices with the unsorted indices.
  */

#ifndef __INDEX_HOLDER_H_
#define __INDEX_HOLDER_H_

namespace seq {

  /**
    @brief A small helper class used while sorting a pair of value and its index in an array/vector. 
    @author Rahul Sampath
    */
  template <typename T> 
    class IndexHolder {
      public:
        T			*value;
        unsigned int		index;

        /** @name Constructors and Destructor */
        //@{
        IndexHolder():value(NULL), index(0) {};
        IndexHolder(T *v, unsigned int i) { value = v; index = i;};
        ~IndexHolder() {};
        //@}

        /**
          @return a < b
          */
        static bool lessThan (const IndexHolder<T>& a, const IndexHolder<T>& b) {
          return a < b;
        }

        /** @name Overload Operators */
        //@{
        bool  operator < ( IndexHolder<T> const  &other) const {
          return (*value < *(other.value) );
        }

        bool  operator > ( IndexHolder<T> const  &other) const {
          return (*value > *(other.value) );
        }

        bool  operator <= ( IndexHolder<T> const  &other) const {
          return (*value <= *(other.value) );
        }

        bool  operator >= ( IndexHolder<T> const  &other) const {
          return (*value >= *(other.value) );
        }

        bool  operator == ( IndexHolder<T> const  &other) const {
          return (*value == *(other.value) );
        }

        bool  operator != ( IndexHolder<T> const  &other) const {
          return (*value != *(other.value) );
        }
        //@}

    };

}//end namespace

#endif

