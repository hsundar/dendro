#ifndef __FE_MAT_H_
#define __FE_MAT_H_

#include "petscdmda.h"
#include "oda.h"

namespace ot {
  namespace fem {

    class feMat {
      public:
        /// typedefs and enums
        enum daType {
          PETSC, OCT
        }; 

        /// Contructors 
        feMat() { };
        feMat(daType da) {
#ifdef __DEBUG__
          assert ( ( da == PETSC ) || ( da == OCT ) );
#endif
          m_daType = da;


        }
        ~feMat() {

        }

        void setDA (DM da) { m_DA = da; }
        void setDA (ot::DA* da) { m_octDA = da; }

        DM  getDA() { return m_DA; }
        ot::DA* getOctDA() { return m_octDA; }
        /**
         * 	@brief		The matrix-vector multiplication routine that is used by
         * 				matrix-free methods. 
         * 	@param		_in	PETSc Vec which is the input vector with whom the 
         * 				product is to be calculated.
         * 	@param		_out PETSc Vec, the output of M*_in
         * 	@return		bool true if successful, false otherwise.
         * 
         *  The matrix-vector multiplication routine that is used by matrix-free 
         * 	methods. The product is directly calculated from the elemental matrices,
         *  which are computed by the ElementalMatVec() function. Use the Assemble()
         *  function for matrix based methods.
        **/ 
        virtual bool MatVec(Vec _in, Vec _out, double scale=1.0) = 0;
        virtual bool MatGetDiagonal(Vec _diag, double scale=1.0) = 0;


        void setProblemDimensions(double x, double y, double z) {
          m_dLx = x;
          m_dLy = y;
          m_dLz = z;
        }
      protected:

        daType          m_daType;

        DM              m_DA;
        ot::DA*         m_octDA;
        /// The dimensions of the problem.
        double m_dLx;
        double m_dLy;
        double m_dLz;
    };

  } // end namespace fem
} // end namespace ot

#endif

