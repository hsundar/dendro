#ifndef __FE_VECTOR_H_
#define __FE_VECTOR_H_

#include <string>
#include "feVec.h"
#include "timeInfo.h"

namespace ot {
  namespace fem {

    template <typename T>
      class feVector : public feVec {
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

          feVector();  
          feVector(daType da);
          ~feVector();

          // operations.
          void setStencil(void* stencil);

          void setName(std::string name);

          virtual bool addVec(Vec _in, double scale=1.0, int indx = -1);

          virtual bool computeVec(Vec _in, Vec _out, double scale = 1.0);

          inline bool ElementalAddVec(int i, int j, int k, PetscScalar ***in, double scale) {
            return asLeaf().ElementalAddVec(i,j,k,in,scale);  
          }

          inline bool ElementalAddVec(unsigned int idx, PetscScalar *in, double scale) {
            return asLeaf().ElementalAddVec(idx, in, scale);  
          }

          inline bool ComputeNodalFunction(int i, int j, int k, PetscScalar ***in, PetscScalar ***out,double scale) {
            return asLeaf().ComputeNodalFunction(i,j,k,in,out,scale);  
          }

          inline bool ComputeNodalFunction(PetscScalar *in, PetscScalar *out,double scale) {
            return asLeaf().ComputeNodalFunction(in, out,scale);  
          }

          /**
           * @brief		Allows static polymorphism access to the derived class. Using the Barton Nackman trick.
           * 
          **/

          T& asLeaf() { return static_cast<T&>(*this);}  

          bool initStencils() {
            return asLeaf().initStencils();
          }

          bool preAddVec() {
            return asLeaf().preAddVec();
          }

          bool postAddVec() {
            return asLeaf().postAddVec();
          }

          bool preComputeVec() {
            return asLeaf().preComputeVec();
          }

          bool postComputeVec() {
            return asLeaf().postComputeVec();
          }

          void setDof(unsigned int dof) { m_uiDof = dof; }
          unsigned int getDof() { return m_uiDof; }

          void setTimeInfo(timeInfo *t) { m_time =t; }
          timeInfo* getTimeInfo() { return m_time; }

          void initOctLut();


          inline PetscErrorCode alignElementAndVertices(ot::DA * da,
              stdElemType & sType, ot::index* indices);
          inline PetscErrorCode mapVtxAndFlagsToOrientation(int childNum, 
              ot::index* indices, unsigned char & mask);
          inline PetscErrorCode reOrderIndices(unsigned char eType, 
              ot::index* indices);

        protected:
          void *          	m_stencil;

          std::string     	m_strVectorType;

          timeInfo		*m_time;

          unsigned int		m_uiDof;

          // Octree specific stuff ...
          unsigned char **	m_ucpLut;
      };

#include "feVector.txx"

  } // end namespace fem
} // end namespace ot

#endif
