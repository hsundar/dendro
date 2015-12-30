
/**
  @file oda.C
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  @author Hari Sundar, hsundar@gmail.com
  */

#include "oda.h"
#include "parUtils.h"
#include "colors.h"
#include "testUtils.h"
#include "dendro.h"

#ifdef __DEBUG__
#ifndef __DEBUG_DA__
#define __DEBUG_DA__
#endif
#endif

#ifdef __DEBUG_DA__
#ifndef __DEBUG_DA_PUBLIC__
#define __DEBUG_DA_PUBLIC__
#endif
#endif

#ifdef __DEBUG_DA_PUBLIC__
#ifndef __MEASURE_DA__
#define __MEASURE_DA__
#endif
#endif

namespace ot {

  int DA::computeLocalToGlobalElemMappings() {
    DendroIntL localElemSize = getElementSize();
    DendroIntL off1, globalOffset;
    MPI_Request sendRequest;
    MPI_Status status;
    if(m_bIamActive) {
      par::Mpi_Scan<DendroIntL>(&localElemSize, &off1, 1, MPI_SUM, m_mpiCommActive); 
      if(m_iRankActive < (m_iNpesActive - 1)) {
        par::Mpi_Issend<DendroIntL>(&off1, 1, m_iRankActive+1, 0, m_mpiCommActive, &sendRequest);
      }

      if(m_iRankActive) {
        par::Mpi_Recv<DendroIntL>(&globalOffset, 1, m_iRankActive-1, 0, m_mpiCommActive, &status );
      }else {
        globalOffset = 0;
      }
    }

    //Equivalent to createVector: elemental, non-ghosted, 1 dof
    std::vector<DendroIntL> gNumNonGhostElems(localElemSize); 

    for(DendroIntL i = 0; i < localElemSize; i++) {
      gNumNonGhostElems[i] = (i+globalOffset);   
    }

    vecGetBuffer<DendroIntL>(gNumNonGhostElems,
        m_dilpLocalToGlobalElems, true, false, true, 1);

    if( m_bIamActive && (m_iRankActive < (m_iNpesActive-1)) ) {
      MPI_Status statusWait;
      MPI_Wait(&sendRequest, &statusWait);
    }

    ReadFromGhostElemsBegin<DendroIntL>(m_dilpLocalToGlobalElems,1);
    ReadFromGhostElemsEnd<DendroIntL>(m_dilpLocalToGlobalElems);

    gNumNonGhostElems.clear();
    m_bComputedLocalToGlobalElems = true;

    return 0;
  }//end function

  int DA::computeLocalToGlobalMappings() {
    DendroIntL localNodeSize = getNodeSize();
    DendroIntL off1, globalOffset;
    MPI_Request sendRequest;
    MPI_Status status;
    if(m_bIamActive) {
      par::Mpi_Scan<DendroIntL>(&localNodeSize, &off1, 1, MPI_SUM, m_mpiCommActive); 
      if(m_iRankActive < (m_iNpesActive-1)) {
        par::Mpi_Issend<DendroIntL>(&off1, 1, m_iRankActive+1, 0, m_mpiCommActive, &sendRequest);
      }

      if(m_iRankActive) {
        par::Mpi_Recv<DendroIntL>(&globalOffset, 1, m_iRankActive-1, 0, m_mpiCommActive, &status);
      }else {
        globalOffset = 0;
      }
    }

    std::vector<DendroIntL> gNumNonGhostNodes(localNodeSize); 
    for(DendroIntL i = 0; i < localNodeSize; i++) {
      gNumNonGhostNodes[i] = (i+globalOffset);   
    }

    vecGetBuffer<DendroIntL>(gNumNonGhostNodes, m_dilpLocalToGlobal,
        false, false, true, 1);

    if(m_bIamActive && (m_iRankActive < (m_iNpesActive-1))) {
      MPI_Status statusWait;
      MPI_Wait(&sendRequest, &statusWait);
    }

    ReadFromGhostsBegin<DendroIntL>(m_dilpLocalToGlobal,1);
    ReadFromGhostsEnd<DendroIntL>(m_dilpLocalToGlobal);

    gNumNonGhostNodes.clear();
    m_bComputedLocalToGlobal = true;

    return 0;
  }//end function

  int DA::setValuesInMatrix(Mat mat, std::vector<ot::MatRecord>& records, unsigned int dof, InsertMode mode) {

    PROF_SET_MAT_VALUES_BEGIN 

      assert(m_bComputedLocalToGlobal);
    std::vector<PetscScalar> values;
    std::vector<PetscInt> colIndices;

    //Can make it more efficient later.
    if(!records.empty()) {
      //Sort Order: row first, col next, val last
      std::sort(records.begin(), records.end());

      unsigned int currRecord = 0;
      while(currRecord < (records.size()-1)) {
        values.push_back(records[currRecord].val);
        colIndices.push_back( static_cast<PetscInt>(
              (dof*m_dilpLocalToGlobal[records[currRecord].colIdx]) +
              records[currRecord].colDim) );
        if( (records[currRecord].rowIdx != records[currRecord+1].rowIdx) ||
            (records[currRecord].rowDim != records[currRecord+1].rowDim) ) {
          PetscInt rowId = static_cast<PetscInt>(
              (dof*m_dilpLocalToGlobal[records[currRecord].rowIdx]) + 
              records[currRecord].rowDim);
          MatSetValues(mat,1,&rowId,colIndices.size(),(&(*colIndices.begin())),
              (&(*values.begin())),mode);
          colIndices.clear();
          values.clear();
        }
        currRecord++;
      }//end while
      PetscInt rowId = static_cast<PetscInt>(
          (dof*m_dilpLocalToGlobal[records[currRecord].rowIdx]) +
          records[currRecord].rowDim);
      if(values.empty()) {
        //Last row is different from the previous row
        PetscInt colId = static_cast<PetscInt>(
            (dof*m_dilpLocalToGlobal[records[currRecord].colIdx]) + 
            records[currRecord].colDim);
        PetscScalar value = records[currRecord].val;
        MatSetValues(mat,1,&rowId,1,&colId,&value,mode);
      }else {
        //Last row is same as the previous row
        values.push_back(records[currRecord].val);
        colIndices.push_back( static_cast<PetscInt>(
              (dof*m_dilpLocalToGlobal[records[currRecord].colIdx]) + 
              records[currRecord].colDim) );
        MatSetValues(mat,1,&rowId,colIndices.size(),(&(*colIndices.begin())),
            (&(*values.begin())),mode);
        colIndices.clear();
        values.clear();
      }
      records.clear();
    }

    PROF_SET_MAT_VALUES_END
  }//end function

  //***************Constructor*****************//
  DA::DA(std::vector<ot::TreeNode> &in, MPI_Comm comm, MPI_Comm activeInputComm,
      bool compressLut, const std::vector<ot::TreeNode>* blocksPtr, bool* iAmActive ) {

#ifdef __PROF_WITH_BARRIER__
    MPI_Barrier(comm);
#endif
    PROF_BUILD_DA_BEGIN

      //@milinda
      int rank;
      MPI_Comm_rank(comm,&rank);

      DA_FactoryPart0(in, comm, activeInputComm, compressLut, iAmActive);

      if(!rank)
          std::cout<<"ODA Part 0 completed"<<std::endl;

    if(m_bIamActive) {
      DA_FactoryPart1(in);
      if(!rank)
            std::cout<<"ODA Part 1 completed"<<std::endl;

       DA_FactoryPart2(in);
       if(!rank)
           std::cout<<"ODA Part 2 completed"<<std::endl;

        std::vector<ot::TreeNode> tmpIn;
        par::sampleSort(in,tmpIn,activeInputComm);
        in=tmpIn;
        tmpIn.clear();


      DA_FactoryPart3(in, comm, compressLut, blocksPtr, iAmActive);

       if(!rank)
           std::cout<<"ODA Part 3 completed"<<std::endl;

    }

    PROF_BUILD_DA_END
  }//end constructor

  DA::DA(unsigned int dummy, std::vector<ot::TreeNode> &in, MPI_Comm comm,
      MPI_Comm activeInputComm, bool compressLut, 
      const std::vector<ot::TreeNode> * blocksPtr, bool* iAmActive ) {

#ifdef __PROF_WITH_BARRIER__
    MPI_Barrier(comm);
#endif
    PROF_BUILD_DA_BEGIN 

      DA_FactoryPart0(in, comm, activeInputComm, compressLut, iAmActive);

    if(m_bIamActive) {
      DA_FactoryPart3(in, comm, compressLut, blocksPtr, iAmActive);
    }

    PROF_BUILD_DA_END
  }//end constructor

  DA::~DA() {
    if (m_ucpOctLevels != NULL) {
      delete [] m_ucpOctLevels;
      m_ucpOctLevels = NULL;
    }

    if(m_dilpLocalToGlobal != NULL) {
      delete [] m_dilpLocalToGlobal;
      m_dilpLocalToGlobal = NULL;
    }

    if(m_dilpLocalToGlobalElems != NULL) {
      delete [] m_dilpLocalToGlobalElems;
      m_dilpLocalToGlobalElems = NULL;
    }

    m_ucpLutRemainders.clear();
    m_uspLutQuotients.clear();
    m_ucpLutMasks.clear();
    m_ucpSortOrders.clear();
    m_uiNlist.clear();
  }

  /************** Domain Access ****************/
  Point DA::getNextOffset(Point p, unsigned char d) {
#ifdef __DEBUG_DA_PUBLIC__
    assert(m_bIamActive);
#endif

    unsigned int len = (unsigned int)(1u<<( m_uiMaxDepth - d ) );
    unsigned int len_par = (unsigned int)(1u<<( m_uiMaxDepth - d +1 ) );

    unsigned int i,j,k;

    i = p.xint(); i %= len_par;
    j = p.yint(); j %= len_par;
    k = p.zint(); k %= len_par;
    i /= len;
    j /= len;
    k /= len;

    unsigned int childNum = 4*k + 2*j + i;

    Point p2;
    switch (childNum) {
      case 7:
        p2.x() = p.x() -len; p2.y() = p.y() - len; p2.z() = p.z() -len;
        return getNextOffset(p2, d-1);
      case 0:
        p2.x() = p.x() +len; p2.y() = p.y(); p2.z() = p.z();
        break;
      case 1:
        p2.x() = p.x() -len; p2.y() = p.y() +len; p2.z() = p.z();
        break;
      case 2:
        p2.x() = p.x() +len; p2.y() = p.y(); p2.z() = p.z();
        break;
      case 3:
        p2.x() = p.x() -len; p2.y() = p.y() - len; p2.z() = p.z() +len;
        break;
      case 4:
        p2.x() = p.x() +len; p2.y() = p.y(); p2.z() = p.z();
        break;
      case 5:
        p2.x() = p.x() -len; p2.y() = p.y()+len; p2.z() = p.z();
        break;
      case 6:
        p2.x() = p.x() +len; p2.y() = p.y(); p2.z() = p.z();
        break;
      default:
        std::cerr << "Wrong child number in " << __func__ << std::endl;
        assert(false);
        break;
    } // switch (childNum)

    return p2;
  }

  void DA::incrementCurrentOffset() {
#ifdef __DEBUG_DA_PUBLIC__
    assert(m_bIamActive);
#endif

    // if it is the first element, simply return the stored offset ...
    if ( m_uiCurrent == (m_uiElementBegin-1)) {
      m_ptCurrentOffset = m_ptOffset;
      return; 
    }

#ifdef __DEBUG_DA_PUBLIC__
    if ( m_ucpOctLevels[m_uiCurrent] & ot::TreeNode::BOUNDARY ) {
      std::cerr << RED "ERROR, Boundary eleme in incre Curr offset" NRM << std::endl;
      assert(false);
    }
#endif

    unsigned char d = (m_ucpOctLevels[m_uiCurrent] & ot::TreeNode::MAX_LEVEL );
    unsigned int len = (unsigned int)(1u<<( m_uiMaxDepth - d ) );
    unsigned int len_par = (unsigned int)(1u<<( m_uiMaxDepth - d +1 ) );

    unsigned int i,j,k;

    i = m_ptCurrentOffset.xint(); 
    j = m_ptCurrentOffset.yint(); 
    k = m_ptCurrentOffset.zint(); 

    i %= len_par;
    j %= len_par;
    k %= len_par;

    i /= len;
    j /= len;
    k /= len;

    unsigned int childNum = 4*k + 2*j + i;

    Point p = m_ptCurrentOffset;
    Point p2;
    switch (childNum) {
      case 7:
        p2.x() = p.x() -len; p2.y() = p.y() - len; p2.z() = p.z() -len;
        p2 = getNextOffset(p2, d-1);
        break;
      case 0:
        p2.x() = p.x() +len; p2.y() = p.y(); p2.z() = p.z();
        break;
      case 1:
        p2.x() = p.x() -len; p2.y() = p.y() +len; p2.z() = p.z();
        break;
      case 2:
        p2.x() = p.x() +len; p2.y() = p.y(); p2.z() = p.z();
        break;
      case 3:
        p2.x() = p.x() -len; p2.y() = p.y() - len; p2.z() = p.z() +len;
        break;
      case 4:
        p2.x() = p.x() +len; p2.y() = p.y(); p2.z() = p.z();
        break;
      case 5:
        p2.x() = p.x() -len; p2.y() = p.y()+len; p2.z() = p.z();
        break;
      case 6:
        p2.x() = p.x() +len; p2.y() = p.y(); p2.z() = p.z();
        break;
      default:
        std::cerr << "Wrong child number in " << __func__ << std::endl;
        assert(false);
        break;
    } // switch (childNum)

    m_ptCurrentOffset = p2;
  }

  //This is for real octants only, pseudo-boundary octants can not be tested
  //using this. 
  bool DA::isBoundaryOctant(unsigned char *flags) {
#ifdef __DEBUG_DA_PUBLIC__
    assert(m_bIamActive);
#endif

    unsigned char _flags = 0;
    Point pt = getCurrentOffset();
    unsigned int x = pt.xint();
    unsigned int y = pt.yint();
    unsigned int z = pt.zint();
    unsigned int d = getLevel(curr())-1;
    unsigned int maxD = getMaxDepth()-1;
    unsigned int len  = (unsigned int)(1u<<(maxD - d) );
    unsigned int blen = (unsigned int)(1u << maxD);

    if (!x) _flags |= ot::TreeNode::X_NEG_BDY;  
    if (!y) _flags |=  ot::TreeNode::Y_NEG_BDY; 
    if (!z) _flags |=   ot::TreeNode::Z_NEG_BDY;

    if ( (x+len) == blen )  _flags |= ot::TreeNode::X_POS_BDY;
    if ( (y+len) == blen )  _flags |= ot::TreeNode::Y_POS_BDY;
    if ( (z+len) == blen )  _flags |= ot::TreeNode::Z_POS_BDY;

    if(flags) {
      *flags = _flags;
    }
    return _flags;
  }//end function

  /***************** Array Access ********************/
  int DA::createMatrix(Mat &M, MatType mtype, unsigned int dof) {
    // first determine the size ...
    unsigned int sz = 0;
    if(m_bIamActive) {
      sz = dof*(m_uiNodeSize + m_uiBoundaryNodeSize);
    }//end if active

    // now create the PETSc Mat
    // The "parallel direct solver" matrix types like MATAIJSPOOLES are ALL gone in petsc-3.0.0
    // Thus, I (Ilya Lashuk) "delete" all such checks for matrix type.  Hope it is reasonable thing to do.
    PetscBool isAij, isAijSeq, isAijPrl, isSuperLU, isSuperLU_Dist;
    PetscStrcmp(mtype,MATAIJ,&isAij);
    PetscStrcmp(mtype,MATSEQAIJ,&isAijSeq);
    PetscStrcmp(mtype,MATMPIAIJ,&isAijPrl);
    isSuperLU = PETSC_FALSE; // PetscStrcmp(mtype,MATSUPERLU,&isSuperLU);
    isSuperLU_Dist = PETSC_FALSE; // PetscStrcmp(mtype,MATSUPERLU_DIST,&isSuperLU_Dist);

    MatCreate(m_mpiCommAll, &M);
    MatSetSizes(M, sz,sz, PETSC_DECIDE, PETSC_DECIDE);
    MatSetType(M,mtype);

    if(isAij || isAijSeq || isAijPrl || isSuperLU || isSuperLU_Dist) {
      if(m_iNpesAll > 1) {
        MatMPIAIJSetPreallocation(M, 53*dof , PETSC_NULL, 53*dof , PETSC_NULL);
      }else {
        MatSeqAIJSetPreallocation(M, 53*dof , PETSC_NULL);
      }
    }

    return 0;
  }//end function

  int DA::createActiveMatrix(Mat &M, MatType mtype, unsigned int dof) {
    // first determine the size ...
    unsigned int sz = 0;
    if(m_bIamActive) {
      sz = dof*(m_uiNodeSize + m_uiBoundaryNodeSize);

      // now create the PETSc Mat
      PetscBool isAij, isAijSeq, isAijPrl, isSuperLU, isSuperLU_Dist;
      PetscStrcmp(mtype,MATAIJ,&isAij);
      PetscStrcmp(mtype,MATSEQAIJ,&isAijSeq);
      PetscStrcmp(mtype,MATMPIAIJ,&isAijPrl);
      isSuperLU = PETSC_FALSE; //PetscStrcmp(mtype,MATSUPERLU,&isSuperLU);
      isSuperLU_Dist = PETSC_FALSE; //PetscStrcmp(mtype,MATSUPERLU_DIST,&isSuperLU_Dist);

      MatCreate(m_mpiCommActive, &M);
      MatSetSizes(M, sz,sz, PETSC_DECIDE, PETSC_DECIDE);
      MatSetType(M,mtype);

      if(isAij || isAijSeq || isAijPrl || isSuperLU || isSuperLU_Dist) {
        if(m_iNpesActive > 1) {
          MatMPIAIJSetPreallocation(M, 53*dof , PETSC_NULL, 53*dof , PETSC_NULL);
        }else {
          MatSeqAIJSetPreallocation(M, 53*dof , PETSC_NULL);
        }
      }
    }//end if active

    return 0;
  }//end function

  int DA::createVector(Vec &arr, bool isElemental, bool isGhosted, unsigned int dof) {
    // first determine the length of the vector ...
    unsigned int sz = 0;
    if(m_bIamActive) {
      if (isElemental) {
        sz = m_uiElementSize;
        if (isGhosted) {
          sz += (m_uiPreGhostElementSize);
        }
      } else {
        sz = m_uiNodeSize + m_uiBoundaryNodeSize;
        if (isGhosted) {
          sz += (m_uiPreGhostNodeSize + m_uiPreGhostBoundaryNodeSize + m_uiPostGhostNodeSize);
        }
      }
      // now for dof ...
      sz *= dof;
    }//end if active

    // now create the PETSc Vector
    VecCreate(m_mpiCommAll, &arr);
    VecSetSizes(arr, sz, PETSC_DECIDE);
    if (m_iNpesAll > 1) {
      VecSetType(arr,VECMPI);
    } else {
      VecSetType(arr,VECSEQ);
    }    
    return 0;
  }

  int DA::createActiveVector(Vec &arr, bool isElemental, bool isGhosted, unsigned int dof) {
    // first determine the length of the vector ...
    unsigned int sz = 0;
    if(m_bIamActive) {
      if (isElemental) {
        sz = m_uiElementSize;
        if (isGhosted) {
          sz += (m_uiPreGhostElementSize);
        }
      } else {
        sz = m_uiNodeSize + m_uiBoundaryNodeSize;
        if (isGhosted) {
          sz += (m_uiPreGhostNodeSize + m_uiPreGhostBoundaryNodeSize + m_uiPostGhostNodeSize);
        }
      }
      // now for dof ...
      sz *= dof;

      // now create the PETSc Vector
      VecCreate(m_mpiCommActive, &arr);
      VecSetSizes(arr, sz, PETSC_DECIDE);
      if (m_iNpesActive > 1) {
        VecSetType(arr, VECMPI);
      } else {
        VecSetType(arr, VECSEQ);
      }    
    }//end if active

    return 0;
  }

  // Obtains a ot::index aligned buffer of the Vector
  int DA::vecGetBuffer(Vec in, PetscScalar* &out, bool isElemental, bool isGhosted,
      bool isReadOnly, unsigned int dof) {
    // Some error checks ... make sure the size of Vec in matches those implied
    // by the other params ...
    unsigned int sz = 0;
    if(m_bIamActive) {
      if (isElemental) {
        sz = m_uiElementSize;
        if (isGhosted) {
          sz += m_uiPreGhostElementSize;
        }
      } else {
        sz = m_uiNodeSize + m_uiBoundaryNodeSize;
        if (isGhosted) {
          sz += (m_uiPreGhostNodeSize + m_uiPreGhostBoundaryNodeSize + m_uiPostGhostNodeSize);
        }
      }
      // now for dof ...
      sz *= dof;
    }//end if active

    PetscInt vecSz=0;
    VecGetLocalSize(in, &vecSz);

    if ( sz != vecSz) {
      std::cerr  << m_iRankAll << ": In function " << __func__ << " sizes are unequal, sz is  " 
        << sz << " and vecSz is " << vecSz << std::endl;
      std::cerr << "Params are: isElem " << isElemental << " isGhosted " << isGhosted << std::endl;
      assert(false);
      return -1;; 
    };

    if(!m_bIamActive) {
      assert(m_uiLocalBufferSize == 0);
      assert(m_uiElementBegin == 0);
      assert(m_uiElementEnd == 0);
      assert(m_uiPostGhostBegin == 0);
    }

    // get the local Petsc Arrray,
    PetscScalar *array = NULL;
    VecGetArray(in, &array);

    // allocate except for the case of ghosted-elemental vectors...
    if(isGhosted && isElemental) {
      //simply copy the pointer
      //This is the only case where the buffer will not be the size of the
      //fullLocalBufferSize. 
      out = array;
    }else {
      // First let us allocate for the buffer ... the local buffer will be of full
      // length.
      sz = dof*m_uiLocalBufferSize;

      if(sz) {
        out = new PetscScalar[sz];
        assert(out);
      }else {
        out = NULL;
      }

      //Zero Entries first if you plan to modify the buffer 
      if(!isReadOnly) {
        for(unsigned int i = 0; i < sz; i++) {
          out[i] = 0.0;
        }
      }
    }

    unsigned int vecCnt=0;
    // Now we can populate the out buffer ... and that needs a loop through the
    // elements ...
    if (isGhosted) {
      if (isElemental) {
        //Nothing to be done here.
      } else {
        // now copy ...
        for (unsigned int i=0; i<m_uiLocalBufferSize; i++) {
          // skip the ones that are not nodes ...
          if ( ! (m_ucpOctLevels[i] & ot::TreeNode::NODE ) ) {
            continue;
          }
          for (unsigned int j=0; j<dof; j++) {
            out[dof*i+j] = array[dof*vecCnt + j];
          }
          vecCnt++;
        }//end for i
      }//end if elemental
    } else {
      if (isElemental) {
        // is a simple copy ...
        for (unsigned int i = m_uiElementBegin; i < m_uiElementEnd; i++) {
          for (unsigned int j = 0; j < dof; j++) {
            out[dof*i+j] = array[dof*vecCnt + j];
          }
          vecCnt++;
        }//end for i
      } else {
        for (unsigned int i = m_uiElementBegin; i < m_uiElementEnd; i++) {
          if ( ! (m_ucpOctLevels[i] & ot::TreeNode::NODE ) ) {
            continue;
          }
          for (unsigned int j=0; j<dof; j++) {
            out[dof*i+j] = array[dof*vecCnt + j];
          }
          vecCnt++;
        }//end for i
        for (unsigned int i = m_uiElementEnd; i < m_uiPostGhostBegin; i++) {
          // add the remaining boundary nodes ...
          if ( ! ( (m_ucpOctLevels[i] & ot::TreeNode::NODE ) &&
                (m_ucpOctLevels[i] & ot::TreeNode::BOUNDARY ) ) ) {
            continue;
          }
          for (unsigned int j=0; j<dof; j++) {
            out[dof*i+j] = array[dof*vecCnt + j];
          }
          vecCnt++;
        }//end for i
      }
    }

    if(!(isGhosted && isElemental)) {
      VecRestoreArray(in, &array);
    }

    return 0;
  }

  int DA::vecRestoreBuffer(Vec in, PetscScalar* out, bool isElemental, bool isGhosted, bool isReadOnly, unsigned int dof) {
    // Some error checks ... make sure the size of Vec in matches those implied
    // by the other params ...
    unsigned int sz = 0;
    if(m_bIamActive) {
      if (isElemental) {
        sz = m_uiElementSize;
        if (isGhosted) {
          sz += m_uiPreGhostElementSize;
        }
      } else {
        sz = m_uiNodeSize + m_uiBoundaryNodeSize;
        if (isGhosted) {
          sz += (m_uiPreGhostNodeSize + m_uiPreGhostBoundaryNodeSize + m_uiPostGhostNodeSize);
        }
      }
      // now for dof ...
      sz *= dof;
    }//end if active

    PetscInt vecSz=0;
    VecGetLocalSize(in, &vecSz);

    if ( sz != vecSz) {
      std::cerr  << RED<<"In function PETSc::" << __func__ <<
        NRM<<" sizes are unequal, sz is  " << sz << " and vecSz is " << vecSz << std::endl;
      std::cerr << "Params are: isElem " << isElemental << " isGhosted " << isGhosted << std::endl;
      assert(false);
      return -1;;
    }

    if(!m_bIamActive) {
      assert(m_uiLocalBufferSize == 0);
      assert(m_uiElementBegin == 0);
      assert(m_uiElementEnd == 0);
      assert(m_uiPostGhostBegin == 0);
    }

    unsigned int vecCnt=0;

    if(isGhosted && isElemental) {
      //If it is ghosted and elemental, simply restore the array.
      //out was not allocated expicitly in this case. It was just a copy of the
      //array's pointer. The readOnly flag is immaterial for this case.
      VecRestoreArray(in, &out);
      out = NULL;
    }  else if ( isReadOnly ) {
      // no need to write back ... simply clean up and return
      //Since this is not an elemental and ghosted vector, out was allocated
      //explicitly 
      if(out) {
        delete [] out;
        out = NULL;
      }
    } else {
      //ghosted and elemental is already taken care of. So only need to tackle
      //the other 3 cases.
      // need to write back ...
      // get the local Petsc Arrray,
      PetscScalar *array;
      VecGetArray(in, &array);

      if ( isElemental ) {
        //non-ghosted, elemental
        for (unsigned int i = m_uiElementBegin; i < m_uiElementEnd; i++) {
          for (unsigned int j=0; j<dof; j++) {
            array[dof*vecCnt + j] = out[dof*i+j];
          }
          vecCnt++;
        }
      } else if ( isGhosted ) {
        // nodal and ghosted ...
        for (unsigned int i=0; i<sz; i++) {
          // skip the ones that are not nodes ...
          if ( ! (m_ucpOctLevels[i] & ot::TreeNode::NODE ) ) {
            continue;
          }
          for (unsigned int j=0; j<dof; j++) {
            array[dof*vecCnt + j] = out[dof*i+j];
          }
          vecCnt++;
        }
      } else {
        // nodal non ghosted ...
        for (unsigned int i = m_uiElementBegin; i < m_uiElementEnd; i++) {
          if ( ! (m_ucpOctLevels[i] & ot::TreeNode::NODE ) ) {
            continue;
          }
          for (unsigned int j=0; j<dof; j++) {
            array[dof*vecCnt + j] = out[dof*i+j];
          }
          vecCnt++;
        }
        for (unsigned int i = m_uiElementEnd; i < m_uiPostGhostBegin; i++) {
          // add the remaining boundary nodes ...
          if ( ! ( (m_ucpOctLevels[i] & ot::TreeNode::NODE ) &&
                (m_ucpOctLevels[i] & ot::TreeNode::BOUNDARY ) ) ) {
            continue;
          }
          for (unsigned int j=0; j<dof; j++) {
            array[dof*vecCnt + j] = out[dof*i+j];
          }
          vecCnt++;
        }
      }

      VecRestoreArray(in, &array);
      //Since this is not an elemental and ghosted vector, out was allocated
      //explicitly 
      if(out) {
        delete [] out;
        out = NULL;
      }
    }
    return 0;
  }

  void DA::updateQuotientCounter() {
#ifdef __DEBUG_DA_PUBLIC__
    assert(m_bIamActive);
#endif

    // m_ucpLutRemainders, m_uspLutQuotients.
    unsigned char _mask = m_ucpLutMasksPtr[2*m_uiCurrent];

    // first let us get the offsets ..
    for (int j=0; j < 8; j++) {
      if ( _mask & (1 << j ) ) {
        m_uiQuotientCounter++;
      }
    }
  }

  unsigned char DA::getHangingNodeIndex(unsigned int i) {
#ifdef __DEBUG_DA_PUBLIC__
    assert(m_bIamActive);
#endif
    return m_ucpLutMasks[2*i + 1];
  }

  void DA::incrementPreGhostOffset() {
#ifdef __DEBUG_DA_PUBLIC__
    assert(m_bIamActive);
#endif

    unsigned char c = m_ucpPreGhostConnectivity[m_uiCurrent];
    if ( c ) {
      unsigned char nd = m_ucpOctLevels[m_uiCurrent+1];
      unsigned char cd = m_ucpOctLevels[m_uiCurrent];
      unsigned int ns = (unsigned int)(1u << ( m_uiMaxDepth - nd ) );
      unsigned int cs = (unsigned int)(1u << ( m_uiMaxDepth - cd ) );

      Point curr = m_ptCurrentOffset;

      unsigned int cx = curr.xint();
      unsigned int cy = curr.yint();
      unsigned int cz = curr.zint();
      unsigned int nx = cx;
      unsigned int ny = cy;
      unsigned int nz = cz;

      //_zzyyxxT
      unsigned char xFlag = ((c & (3<<1) ) >> 1);
      unsigned char yFlag = ((c & (3<<3) ) >> 3);
      unsigned char zFlag = ((c & (3<<5) ) >> 5);

      switch (xFlag) {
        case 0: nx = cx;
                break;
        case 1: nx = (cx - ns);
                break;
        case 2: nx = (cx + cs); 
                break;
        case 3: nx = (cx + cs - ns);
                break;
        default: assert(false);
      }

      switch (yFlag) {
        case 0: ny = cy;
                break;
        case 1: ny = (cy - ns);
                break;
        case 2: ny = (cy + cs); 
                break;
        case 3: ny = (cy + cs - ns); 
                break;
        default: assert(false);
      }

      switch (zFlag) {
        case 0: nz = cz;
                break;
        case 1: nz = (cz - ns); 
                break;
        case 2: nz = (cz + cs); 
                break;
        case 3: nz = (cz + cs - ns); 
                break;
        default: assert(false);
      }

      m_ptCurrentOffset = Point(nx,ny,nz);    

    } else {
      m_ptCurrentOffset = m_ptsPreGhostOffsets[m_uiPreGhostQuotientCnt++];     
    }

  }//end function

} // end namespace ot


