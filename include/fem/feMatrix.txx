
template <typename T>
feMatrix<T>::feMatrix() {
  m_daType = PETSC;
  m_DA 		= NULL;
  m_octDA 	= NULL;
  m_stencil	= NULL;
  m_uiDof	= 1;
  m_ucpLut	= NULL;

  // initialize the stencils ...
  initStencils();
}

template <typename T>
feMatrix<T>::feMatrix(daType da) {
#ifdef __DEBUG__
  assert ( ( da == PETSC ) || ( da == OCT ) );
#endif
  m_daType = da;
  m_DA 		= NULL;
  m_octDA 	= NULL;
  m_stencil	= NULL;
  m_ucpLut	= NULL;

  // initialize the stencils ...
  initStencils();
  if (da == OCT)
    initOctLut();
}

template <typename T>
void feMatrix<T>::initOctLut() {
  //Note: It is not symmetric.
  unsigned char tmp[8][8]={
    {0,1,2,3,4,5,6,7},
    {1,3,0,2,5,7,4,6},
    {2,0,3,1,6,4,7,5},
    {3,2,1,0,7,6,5,4},
    {4,5,0,1,6,7,2,3},
    {5,7,1,3,4,6,0,2},
    {6,4,2,0,7,5,3,1},
    {7,6,3,2,5,4,1,0}
  };

  //Is Stored in  ROW_MAJOR Format.  
  typedef unsigned char* charPtr;
  m_ucpLut = new charPtr[8];
  for (int i=0;i<8;i++) {
    m_ucpLut[i] = new unsigned char[8]; 
    for (int j=0;j<8;j++) {
      m_ucpLut[i][j] = tmp[i][j];
    }
  }
}

template <typename T>
feMatrix<T>::~feMatrix() {
}


#undef __FUNCT__
#define __FUNCT__ "feMatrix_MatGetDiagonal"
template <typename T>
bool feMatrix<T>::MatGetDiagonal(Vec _diag, double scale){
  PetscFunctionBegin;
#ifdef __DEBUG__
  assert ( ( m_daType == PETSC ) || ( m_daType == OCT ) );
#endif

  int ierr;

  // PetscScalar zero=0.0;

  if (m_daType == PETSC) {

    PetscInt x,y,z,m,n,p;
    PetscInt mx,my,mz;
    int xne,yne,zne;

    PetscScalar ***diag;
    Vec diagLocal;

    /* Get all corners*/
    if (m_DA == NULL)
      std::cerr << "Da is null" << std::endl;
    ierr = DAGetCorners(m_DA, &x, &y, &z, &m, &n, &p); CHKERRQ(ierr); 
    /* Get Info*/
    ierr = DAGetInfo(m_DA,0, &mx, &my, &mz, 0,0,0,0,0,0,0); CHKERRQ(ierr); 

    if (x+m == mx) {
      xne=m-1;
    } else {
      xne=m;
    }
    if (y+n == my) {
      yne=n-1;
    } else {
      yne=n;
    }
    if (z+p == mz) {
      zne=p-1;
    } else {
      zne=p;
    }

    ierr = DAGetLocalVector(m_DA, &diagLocal); CHKERRQ(ierr);
    ierr = VecZeroEntries(diagLocal);

    // ierr = DAGlobalToLocalBegin(m_DA, _diag, INSERT_VALUES, diagLocal); CHKERRQ(ierr);
    // ierr = DAGlobalToLocalEnd(m_DA, _diag, INSERT_VALUES, diagLocal); CHKERRQ(ierr);

   
    ierr = DAVecGetArray(m_DA, diagLocal, &diag);

    // Any derived class initializations ...
    preMatVec();

    // loop through all elements ...
    for (int k=z; k<z+zne; k++){
      for (int j=y; j<y+yne; j++){
        for (int i=x; i<x+xne; i++){
         ElementalMatGetDiagonal(i, j, k, diag, scale);
        } // end i
      } // end j
    } // end k

    postMatVec();

    ierr = DAVecRestoreArray(m_DA, diagLocal, &diag);	CHKERRQ(ierr);  

   
    ierr = DALocalToGlobalBegin(m_DA, diagLocal, _diag); CHKERRQ(ierr);  
    ierr = DALocalToGlobalEnd(m_DA, diagLocal, _diag); CHKERRQ(ierr);  
    
    ierr = DARestoreLocalVector(m_DA, &diagLocal); CHKERRQ(ierr);  

    
  } else {
    // loop for octree DA.
    PetscScalar *diag=NULL;

    // get Buffers ...
    //Nodal,Non-Ghosted,Read,1 dof, Get in array and get ghosts during computation
    m_octDA->vecGetBuffer(_diag, diag, false, false, false, m_uiDof);
    
    preMatVec();

    // loop through all elements ...
    for ( m_octDA->init<ot::DA_FLAGS::ALL>(); m_octDA->curr() < m_octDA->end<ot::DA_FLAGS::ALL>(); m_octDA->next<ot::DA_FLAGS::ALL>() ) {
      ElementalMatGetDiagonal( m_octDA->curr(), diag, scale); 
    }//end 

    postMatVec();

    // Restore Vectors ..
    m_octDA->vecRestoreBuffer(_diag, diag, false, false, false, m_uiDof);
  }

  PetscFunctionReturn(0);
}



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
 *  which are computed by the ElementalMatrix() function. Use the Assemble()
 *  function for matrix based methods.
 **/
#undef __FUNCT__
#define __FUNCT__ "feMatrix_MatVec"
template <typename T>
bool feMatrix<T>::MatVec(Vec _in, Vec _out, double scale){
  PetscFunctionBegin;

#ifdef __DEBUG__
  assert ( ( m_daType == PETSC ) || ( m_daType == OCT ) );
#endif

  int ierr;
  // PetscScalar zero=0.0;

  if (m_daType == PETSC) {

    PetscInt x,y,z,m,n,p;
    PetscInt mx,my,mz;
    int xne,yne,zne;

    PetscScalar ***in, ***out;
    Vec inlocal, outlocal;

    /* Get all corners*/
    if (m_DA == NULL)
      std::cerr << "Da is null" << std::endl;
    ierr = DAGetCorners(m_DA, &x, &y, &z, &m, &n, &p); CHKERRQ(ierr); 
    /* Get Info*/
    ierr = DAGetInfo(m_DA,0, &mx, &my, &mz, 0,0,0,0,0,0,0); CHKERRQ(ierr); 

    if (x+m == mx) {
      xne=m-1;
    } else {
      xne=m;
    }
    if (y+n == my) {
      yne=n-1;
    } else {
      yne=n;
    }
    if (z+p == mz) {
      zne=p-1;
    } else {
      zne=p;
    }

    // std::cout << x << "," << y << "," << z << " + " << xne <<","<<yne<<","<<zne<<std::endl;

    // Get the local vector so that the ghost nodes can be accessed
    ierr = DAGetLocalVector(m_DA, &inlocal); CHKERRQ(ierr);
    ierr = DAGetLocalVector(m_DA, &outlocal); CHKERRQ(ierr);
    // ierr = VecDuplicate(inlocal, &outlocal); CHKERRQ(ierr);

    ierr = DAGlobalToLocalBegin(m_DA, _in, INSERT_VALUES, inlocal); CHKERRQ(ierr);
    ierr = DAGlobalToLocalEnd(m_DA, _in, INSERT_VALUES, inlocal); CHKERRQ(ierr);
    // ierr = DAGlobalToLocalBegin(m_DA, _out, INSERT_VALUES, outlocal); CHKERRQ(ierr);
    // ierr = DAGlobalToLocalEnd(m_DA, _out, INSERT_VALUES, outlocal); CHKERRQ(ierr);

    ierr = VecZeroEntries(outlocal);

    ierr = DAVecGetArray(m_DA, inlocal, &in);
    ierr = DAVecGetArray(m_DA, outlocal, &out);

    // Any derived class initializations ...
    preMatVec();

    // loop through all elements ...
    for (int k=z; k<z+zne; k++){
      for (int j=y; j<y+yne; j++){
        for (int i=x; i<x+xne; i++){
          // std::cout << i <<"," << j << "," << k << std::endl;
          ElementalMatVec(i, j, k, in, out, scale);
        } // end i
      } // end j
    } // end k

    postMatVec();

    ierr = DAVecRestoreArray(m_DA, inlocal, &in); CHKERRQ(ierr);  
    ierr = DAVecRestoreArray(m_DA, outlocal, &out);	CHKERRQ(ierr);  

    ierr = DALocalToGlobalBegin(m_DA, outlocal, _out); CHKERRQ(ierr);  
    ierr = DALocalToGlobalEnd(m_DA, outlocal, _out); CHKERRQ(ierr);  
    
    ierr = DARestoreLocalVector(m_DA, &inlocal); CHKERRQ(ierr);  
    ierr = DARestoreLocalVector(m_DA, &outlocal); CHKERRQ(ierr);  
    // ierr = VecDestroy(outlocal); CHKERRQ(ierr);  

  } else {
    // loop for octree DA.
    

    PetscScalar *out=NULL;
    PetscScalar *in=NULL; 

    // get Buffers ...
    //Nodal,Non-Ghosted,Read,1 dof, Get in array and get ghosts during computation
    m_octDA->vecGetBuffer(_in,   in, false, false, true,  m_uiDof);
    m_octDA->vecGetBuffer(_out, out, false, false, false, m_uiDof);
    
    // start comm for in ...
    //m_octDA->updateGhostsBegin<PetscScalar>(in, false, m_uiDof);
    // m_octDA->ReadFromGhostsBegin<PetscScalar>(in, false, m_uiDof);
    m_octDA->ReadFromGhostsBegin<PetscScalar>(in, m_uiDof);
    preMatVec();

    // Independent loop, loop through the nodes this processor owns..
    for ( m_octDA->init<ot::DA_FLAGS::INDEPENDENT>(), m_octDA->init<ot::DA_FLAGS::WRITABLE>(); m_octDA->curr() < m_octDA->end<ot::DA_FLAGS::INDEPENDENT>(); m_octDA->next<ot::DA_FLAGS::INDEPENDENT>() ) {
      ElementalMatVec( m_octDA->curr(), in, out, scale); 
    }//end INDEPENDENT

    // Wait for communication to end.
    //m_octDA->updateGhostsEnd<PetscScalar>(in);
	 m_octDA->ReadFromGhostsEnd<PetscScalar>(in);

    // Dependent loop ...
    for ( m_octDA->init<ot::DA_FLAGS::DEPENDENT>(), m_octDA->init<ot::DA_FLAGS::WRITABLE>(); m_octDA->curr() < m_octDA->end<ot::DA_FLAGS::DEPENDENT>(); m_octDA->next<ot::DA_FLAGS::DEPENDENT>() ) {
      ElementalMatVec( m_octDA->curr(), in, out, scale); 
    }//end DEPENDENT

    postMatVec();

    // Restore Vectors ...
    m_octDA->vecRestoreBuffer(_in,   in, false, false, true,  m_uiDof);
    m_octDA->vecRestoreBuffer(_out, out, false, false, false, m_uiDof);

  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "alignElementAndVertices"
template <typename T>
PetscErrorCode feMatrix<T>::alignElementAndVertices(ot::DA * da, stdElemType & sType, unsigned int* indices) {
  PetscFunctionBegin;
  
  sType = ST_0;
  da->getNodeIndices(indices); 

  // not required ....
  // int rank;
  // MPI_Comm_rank(da->getComm(), &rank);

  if (da->isHanging(da->curr())) {

    int childNum = da->getChildNumber();
    Point pt = da->getCurrentOffset();   

    unsigned char hangingMask = da->getHangingNodeIndex(da->curr());    

    //Change HangingMask and indices based on childNum
    mapVtxAndFlagsToOrientation(childNum, indices, hangingMask);    

    unsigned char eType = ((126 & hangingMask)>>1);

    reOrderIndices(eType, indices);
  }//end if hangingElem.
  PetscFunctionReturn(0);
}//end function.

#undef __FUNCT__
#define __FUNCT__ "mapVtxAndFlagsToOrientation"
template <typename T>
PetscErrorCode feMatrix<T>::mapVtxAndFlagsToOrientation(int childNum, 
 unsigned int* indices, unsigned char & mask) {
  PetscFunctionBegin;
  unsigned int tmp[8];
  unsigned char tmpFlags = 0;
  for (int i = 0; i < 8; i++) {
    tmp[i] = indices[m_ucpLut[childNum][i]];
    tmpFlags = ( tmpFlags | ( ( (1 << (m_ucpLut[childNum][i])) & mask ) ? (1 << i) : 0 ) );
  }
  for (int i=0;i<8;i++) {
    indices[i] = tmp[i];
  }
  mask = tmpFlags;
  PetscFunctionReturn(0);
}//end function

#undef __FUNCT__
#define __FUNCT__ "reOrderIndices"
template <typename T>
PetscErrorCode feMatrix<T>::reOrderIndices(unsigned char eType, unsigned int* indices) {
#ifdef __DEBUG_1
  std::cout << "Entering " << __func__ << std::endl;
#endif
  PetscFunctionBegin;
  unsigned int tmp;
  switch (eType) {
  case  ET_N: 
    break;
  case  ET_Y:
    break;
  case  ET_X:
    //Swap 1 & 2, Swap 5 & 6
    tmp = indices[1];
    indices[1] = indices[2];
    indices[2] = tmp;
    tmp = indices[5];
    indices[5] = indices[6];
    indices[6] = tmp;
    break;
  case  ET_XY:
    break;
  case  ET_Z:
    //Swap 2 & 4, Swap 3 & 5
    tmp = indices[2];
    indices[2] = indices[4];
    indices[4] = tmp;
    tmp = indices[3];
    indices[3] = indices[5];
    indices[5] = tmp;
    break;
  case  ET_ZY:
    //Swap 1 & 4, Swap 3 & 6
    tmp = indices[1];
    indices[1] = indices[4];
    indices[4] = tmp;
    tmp = indices[3];
    indices[3] = indices[6];
    indices[6] = tmp;
    break;
  case  ET_ZX:
    //Swap 2 & 4, Swap 3 & 5
    tmp = indices[2];
    indices[2] = indices[4];
    indices[4] = tmp;
    tmp = indices[3];
    indices[3] = indices[5];
    indices[5] = tmp;
    break;
  case  ET_ZXY:
    break;
  case  ET_XY_XY:
    break;
  case  ET_XY_ZXY:
    break;
  case  ET_YZ_ZY:
    //Swap 1 & 4, Swap 3 & 6
    tmp = indices[1];
    indices[1] = indices[4];
    indices[4] = tmp;
    tmp = indices[3];
    indices[3] = indices[6];
    indices[6] = tmp;
    break;
  case  ET_YZ_ZXY:
    //Swap 1 & 4, Swap 3 & 6
    tmp = indices[1];
    indices[1] = indices[4];
    indices[4] = tmp;
    tmp = indices[3];
    indices[3] = indices[6];
    indices[6] = tmp;
    break;
  case  ET_YZ_XY_ZXY:
    break;
  case  ET_ZX_ZX:
    //Swap 2 & 4, Swap 3 & 5
    tmp = indices[2];
    indices[2] = indices[4];
    indices[4] = tmp;
    tmp = indices[3];
    indices[3] = indices[5];
    indices[5] = tmp;
    break;
  case  ET_ZX_ZXY:
    //Swap 2 & 4, Swap 3 & 5
    tmp = indices[2];
    indices[2] = indices[4];
    indices[4] = tmp;
    tmp = indices[3];
    indices[3] = indices[5];
    indices[5] = tmp;
    break;
  case  ET_ZX_XY_ZXY:
    //Swap 1 & 2, Swap 5 & 6
    tmp = indices[1];
    indices[1] = indices[2];
    indices[2] = tmp;
    tmp = indices[5];
    indices[5] = indices[6];
    indices[6] = tmp;
    break;
  case  ET_ZX_YZ_ZXY:
    //Swap 2 & 4, Swap 3 & 5
    tmp = indices[2];
    indices[2] = indices[4];
    indices[4] = tmp;
    tmp = indices[3];
    indices[3] = indices[5];
    indices[5] = tmp;
    break;
  case  ET_ZX_YZ_XY_ZXY:
    break;
  default:
    std::cout<<"in reOrder Etype: "<< (int) eType << std::endl;
    assert(false);
  }
#ifdef __DEBUG_1
  std::cout << "Leaving " << __func__ << std::endl;
#endif
  PetscFunctionReturn(0);
}

