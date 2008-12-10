#ifndef __FE_MATRIX_H_
#define __FE_MATRIX_H_

#include <string>
#include "feMat.h"
#include "timeInfo.h"

namespace ot {
  namespace fem {

    template <typename T>
      class feMatrix : public feMat {
        public:  

          enum stdElemType {
            ST_0,ST_1,ST_2,ST_3,ST_4,ST_5,ST_6,ST_7
          };

          enum exhaustiveElemType {
            //Order: 654321
            //YZ ZX Z XY Y X
            ET_N = 0,
            ET_Y = 2,
            ET_X = 1,
            ET_XY = 3,
            ET_Z = 8,
            ET_ZY = 10,
            ET_ZX = 9,
            ET_ZXY = 11,
            ET_XY_XY = 7,
            ET_XY_ZXY = 15,
            ET_YZ_ZY = 42,
            ET_YZ_ZXY = 43,
            ET_YZ_XY_ZXY = 47,
            ET_ZX_ZX = 25,
            ET_ZX_ZXY = 27,
            ET_ZX_XY_ZXY = 31,
            ET_ZX_YZ_ZXY = 59,
            ET_ZX_YZ_XY_ZXY = 63
          };

          feMatrix();  
          feMatrix(daType da);
          ~feMatrix();

          // operations.
          void setStencil(void* stencil);

          void setName(std::string name);

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
          virtual bool MatVec(Vec _in, Vec _out, double scale=1.0);

          virtual bool MatGetDiagonal(Vec _diag, double scale=1.0);

          /**
           * 	@brief		The elemental matrix-vector multiplication routine that is used
           *				by matrix-free methods. 
           * 	@param		_in	PETSc Vec which is the input vector with whom the 
           * 				product is to be calculated.
           * 	@param		_out PETSc Vec, the output of M*_in
           * 	@return		bool true if successful, false otherwise.
           *  @todo		Might have to change _in and _out to std. C arrays for speed.
           *
           *  The implementation for this function shall be in derived classes, based on
           * 	the problem formulation. Look at MassMatrix and StiffnessMatrix for standard
           * 	implementations. 
          **/ 
          inline bool ElementalMatVec(int i, int j, int k, PetscScalar ***in, PetscScalar ***out, double scale) {
            return asLeaf().ElementalMatVec(i,j,k,in,out,scale);  
          }

          inline bool ElementalMatGetDiagonal(int i, int j, int k, PetscScalar ***diag, double scale) {
            return asLeaf().ElementalMatGetDiagonal(i,j,k,diag,scale);  
          }

          /**
           * 	@brief		The elemental matrix-vector multiplication routine that is used
           *				by matrix-free methods. 
           * 	@param		_in	PETSc Vec which is the input vector with whom the 
           * 				product is to be calculated.
           * 	@param		_out PETSc Vec, the output of M*_in
           * 	@return		bool true if successful, false otherwise.
           *  @todo		Might have to change _in and _out to std. C arrays for speed.
           *
           *  The implementation for this function shall be in derived classes, based on
           * 	the problem formulation. Look at MassMatrix and StiffnessMatrix for standard
           * 	implementations. 
          **/ 
          inline bool ElementalMatVec(unsigned int idx, PetscScalar *in, PetscScalar *out, double scale) {
            return asLeaf().ElementalMatVec(idx, in, out, scale);  
          }

          inline bool ElementalMatGetDiagonal(unsigned int idx, PetscScalar *diag, double scale) {
            return asLeaf().ElementalMatGetDiagonal(idx, diag, scale);  
          }

          // PetscErrorCode matVec(Vec in, Vec out, timeInfo info);

          /**
           * @brief		Allows static polymorphism access to the derived class. Using the Barton Nackman trick.
           * 
          **/

          T& asLeaf() { return static_cast<T&>(*this);}  

          bool initStencils() {
            return asLeaf().initStencils();
          }

          bool preMatVec() {
            return asLeaf().preMatVec();
          }

          bool postMatVec() {
            return asLeaf().postMatVec();
          }

          void setDof(unsigned int dof) { m_uiDof = dof; }
          unsigned int getDof() { return m_uiDof; }

          void setTimeInfo(timeInfo *t) { m_time =t; }
          timeInfo* getTimeInfo() { return m_time; }

          void initOctLut();


          inline PetscErrorCode alignElementAndVertices(ot::DA * da,
              stdElemType & sType, unsigned int* indices);
          inline PetscErrorCode mapVtxAndFlagsToOrientation(int childNum, 
              unsigned int* indices, unsigned char & mask);
          inline PetscErrorCode reOrderIndices(unsigned char eType,
              unsigned int* indices);

        protected:
          void *          	m_stencil;

          std::string     	m_strMatrixType;

          timeInfo		*m_time;

          unsigned int		m_uiDof;

          // Octree specific stuff ...
          unsigned char **	m_ucpLut;
      };

#include "feMatrix.txx"

  } // end namespace fem
} // end namespace ot
#endif
