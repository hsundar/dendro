/**
 *  @file	stiffnessMatrix.h
 *  @brief	Main class for finite element assembly of a generic
 *           stiffness matrix.
 *  @author	Hari Sundar
 *  @date	5/8/7
 *
 *  Main class for finite element assembly of a generic stiffness matrix. The Lame 
 *  parameter images, i.e., the material properties, \f$\mu\f$ and \f$\lambda\f$, need to
 *  be specified. MatVec functions are virtual and in most cases need not be 
 *  specified as they are obtained from the parent feMatrix class. In most cases 
 *  only the assembly of the element matrix needs to be  done within this or a 
 *  derived class.
 */
#ifndef _STIFFNESSMATRIX_H_
#define _STIFFNESSMATRIX_H_

#include "feMatrix.h"

/**
 *  @brief	Main class for finite element assembly of a generic stiffness matrix.
 *  @author	Hari Sundar
 *  @date	5/8/7
 *
 *  Main class for finite element assembly of a generic stiffness matrix. The Lame 
 *  parameter images, i.e., the material properties, \f$\mu\f$ and \f$\lambda\f$ need to
 *  be specified. MatVec functions are virtual and in most cases need not be 
 *  specified as they are obtained from the parent feMatrix class. In most cases 
 *  only the assembly of the element matrix needs to be  done within this or a 
 *  derived class.
 */

class stiffnessMatrix : public ot::fem::feMatrix<stiffnessMatrix>
{
  public:	
    stiffnessMatrix(daType da);
    /**
     * 	@brief		The elemental matrix-vector multiplication routine that is used
     *				by matrix-free methods. 
     * 	@param		_in	PETSc Vec which is the input vector with whom the 
     * 				product is to be calculated.
     * 	@param		_out PETSc Vec, the output of M*_in
     * 	@return		bool true if successful, false otherwise.
     *  @todo		Might have to change _in and _out to std. C arrays for speed.
     * 
     **/ 
    inline bool ElementalMatVec(int i, int j, int k, PetscScalar ***in, PetscScalar ***out, double scale);
    inline bool ElementalMatVec(unsigned int idx, PetscScalar *in, PetscScalar *out, double scale);

    inline bool ElementalMatGetDiagonal(int i, int j, int k, PetscScalar ***diag, double scale);
    inline bool ElementalMatGetDiagonal(unsigned int idx, PetscScalar *diag, double scale);
    
    inline bool initStencils();

    bool preMatVec();
    bool postMatVec();

    void setNuVec(Vec nv) {
      nuvec = nv;
    }

   private:
    void*      		m_nuarray; /* Diffusion coefficient array*/
    Vec                 nuvec;     

    double 		m_dHx;
	 double     m_nuval;
    double xFac, yFac, zFac;
    unsigned int maxD;

};

stiffnessMatrix::stiffnessMatrix(daType da) {
#ifdef __DEBUG__
  assert ( ( da == PETSC ) || ( da == OCT ) );
#endif
  m_daType = da;
  m_DA 		= NULL;
  m_octDA 	= NULL;
  m_stencil	= NULL;

  // initialize the stencils ...
  initStencils();
  if (da == OCT)
    initOctLut();
}

bool stiffnessMatrix::initStencils() {
  typedef int* int1Ptr;
  typedef int** int2Ptr;

  if ( m_daType == PETSC) {
    int Bjk[8][8] = {
      { 64,   0,   0, -16,   0, -16, -16, -16},
      {0,  64, -16,   0, -16,   0, -16, -16},
      {0, -16,  64,   0, -16, -16,   0, -16},
      { -16,   0,   0,  64, -16, -16, -16,   0},
      { 0, -16, -16, -16,  64,   0,   0, -16},
      { -16,   0, -16, -16,   0,  64, -16,   0},
      { -16, -16,   0, -16,   0, -16,  64,   0},
      { -16, -16, -16,   0, -16,   0,   0,  64}
    };
    
    int** Ajk = new int1Ptr[8];
    for (int j=0;j<8;j++) {
      Ajk[j] = new int[8];
      for (int k=0;k<8;k++) {
        Ajk[j][k] = Bjk[j][k];
      }//end k
    }//end j
    m_stencil = Ajk;

  } else {
    int Bijk[8][8][8] = {
      //Type-0: No Hanging
      {{ 64,   0,   0, -16,   0, -16, -16, -16},
        {0,  64, -16,   0, -16,   0, -16, -16},
        {0, -16,  64,   0, -16, -16,   0, -16},
        { -16,   0,   0,  64, -16, -16, -16,   0},
        { 0, -16, -16, -16,  64,   0,   0, -16},
        { -16,   0, -16, -16,   0,  64, -16,   0},
        { -16, -16,   0, -16,   0, -16,  64,   0},
        { -16, -16, -16,   0, -16,   0,   0,  64}},
      //Type-1: Y Hanging
      {{  80,  -8,  16, -16,  -8, -24, -16, -24},
        {  -8,  64,  -8,   0, -16,   0, -16, -16},
        {  16,  -8,  16,   0,  -8,  -8,   0,  -8},
        { -16,   0,   0,  64, -16, -16, -16,   0},
        {  -8, -16,  -8, -16,  64,   0,   0, -16},
        { -24,   0,  -8, -16,   0,  64, -16,   0},
        { -16, -16,   0, -16,   0, -16,  64,   0},
        { -24, -16,  -8,   0, -16,   0,   0,  64}},
      //Type-2: X and Y Hanging
      {{  88,  12,  12, -16, -16, -24, -24, -32},
        {  12,  16,  -4,   0,  -8,   0,  -8,  -8},
        {  12,  -4,  16,   0,  -8,  -8,   0,  -8},
        { -16,   0,   0,  64, -16, -16, -16,   0},
        { -16,  -8,  -8, -16,  64,   0,   0, -16},
        { -24,   0,  -8, -16,   0,  64, -16,   0},
        { -24,  -8,   0, -16,   0, -16,  64,   0},
        { -32,  -8,  -8,   0, -16,   0,   0,  64}},
      //Type-3: X, Y and Z Hanging
      {{  88,   8,   8, -24,   8, -24, -24, -40},
        {   8,  16,  -4,   0,  -4,   0,  -8,  -8},
        {   8,  -4,  16,   0,  -4,  -8,   0,  -8},
        { -24,   0,   0,  64,  -8, -16, -16,   0},
        {   8,  -4,  -4,  -8,  16,   0,   0,  -8},
        { -24,   0,  -8, -16,   0,  64, -16,   0},
        { -24,  -8,   0, -16,   0, -16,  64,   0},
        { -40,  -8,  -8,   0,  -8,   0,   0,  64}},
      //Type-4: XY and X and Y Hanging
      {{  84,  12,  12,   0, -20, -28, -28, -32},
        {  12,  20,   0,   4, -12,  -4, -12,  -8},
        {  12,   0,  20,   4, -12, -12,  -4,  -8},
        {   0,   4,   4,   4,  -4,  -4,  -4,   0},
        { -20, -12, -12,  -4,  64,   0,   0, -16},
        { -28,  -4, -12,  -4,   0,  64, -16,   0},
        { -28, -12,  -4,  -4,   0, -16,  64,   0},
        { -32,  -8,  -8,   0, -16,   0,   0,  64}},
      //Type-5: XY and X and Y and Z Hanging
      {{  80,   6,   6,  -2,   6, -28, -28, -40},
        {   6,  20,   0,   4,  -6,  -4, -12,  -8},
        {   6,   0,  20,   4,  -6, -12,  -4,  -8},
        {  -2,   4,   4,   4,  -2,  -4,  -4,   0},
        {   6,  -6,  -6,  -2,  16,   0,   0,  -8},
        { -28,  -4, -12,  -4,   0,  64, -16,   0},
        { -28, -12,  -4,  -4,   0, -16,  64,   0},
        { -40,  -8,  -8,   0,  -8,   0,   0,  64}},
      //Type-6:XY and YZ and X and Y and Z Hanging
      {{  70,   3,   2,  -3,   3, -32,  -3, -40},
        {   3,  20,  -3,   4,  -9,  -4,  -3,  -8},
        {   2,  -3,  22,   3,  -3, -16,   3,  -8},
        {  -3,   4,   3,   4,  -3,  -4,  -1,   0},
        {   3,  -9,  -3,  -3,  20,  -4,   4,  -8},
        { -32,  -4, -16,  -4,  -4,  64,  -4,   0},
        {  -3,  -3,   3,  -1,   4,  -4,   4,   0},
        { -40,  -8,  -8,   0,  -8,   0,   0,  64}},
      //Type-7:All 6 Hanging
      {{  58,  -2,  -2,  -4,  -2,  -4,  -4, -40},
        {  -2,  22,  -7,   3,  -7,   3,  -4,  -8},
        {  -2,  -7,  22,   3,  -7,  -4,   3,  -8},
        {  -4,   3,   3,   4,  -4,  -1,  -1,   0},
        {  -2,  -7,  -7,  -4,  22,   3,   3,  -8},
        {  -4,   3,  -4,  -1,   3,   4,  -1,   0},
        {  -4,  -4,   3,  -1,   3,  -1,   4,   0},
        { -40,  -8,  -8,   0,  -8,   0,   0,  64}}
    };
    int ***Aijk = new int2Ptr[8];
    for (int i=0;i<8;i++) {
      Aijk[i] = new int1Ptr[8];
      for (int j=0;j<8;j++) {
        Aijk[i][j] = new int[8];
        for (int k=0;k<8;k++) {
          Aijk[i][j][k] = Bijk[i][j][k];
        }//end k
      }//end j
    }//end i
    m_stencil = Aijk;
  }
  return true;
}

bool stiffnessMatrix::preMatVec() {
  // nuVec should be set directly into matrix outside the loop ...
  if (m_daType == PETSC) {
    PetscScalar ***nuarray; //
	 // int ierr;
	 // ierr = VecNorm(nuvec,NORM_INFINITY,&m_nuval); CHKERRQ(ierr);

    int ierr = DAVecGetArray(m_DA, nuvec, &nuarray);

    m_nuarray = nuarray;
    // compute Hx
    PetscInt mx,my,mz;
    ierr = DAGetInfo(m_DA,0, &mx, &my, &mz, 0,0,0,0,0,0,0); CHKERRQ(ierr); 
    CHKERRQ(ierr);
    
    m_dHx = m_dLx/(mx -1);

    m_dHx /= 192.0;
  } else {
    PetscScalar *nuarray; 
    // Get nuarray
    m_octDA->vecGetBuffer(nuvec, nuarray, false, false, true,m_uiDof);
    m_nuarray = nuarray;

    // compute Hx
    // For octree Hx values will change per element, so has to be 
    // computed inside the loop.
	 maxD = m_octDA->getMaxDepth();

    
    // Get the  x,y,z factors 
    xFac = 1.0/((double)(1u << (maxD-1)));
    if (m_octDA->getDimension() > 1) {
      yFac = 1.0/((double)(1u << (maxD-1)));
      if (m_octDA->getDimension() > 2) {
        zFac = 1.0/((double)(1u << (maxD-1)));
      }
    }
  }

  return true;
}

bool stiffnessMatrix::ElementalMatVec(unsigned int i, PetscScalar *in, PetscScalar *out, double scale) {
  unsigned int lev = m_octDA->getLevel(i);
  double hx = xFac*(1u << (maxD - lev));
  double hy = yFac*(1u << (maxD - lev));
  double hz = zFac*(1u << (maxD - lev));

  double fac11 = -hx*scale/192.0;
  
  stdElemType elemType;
  unsigned int idx[8];

  unsigned char hangingMask = m_octDA->getHangingNodeIndex(i);    //  alignElementAndVertices(m_octDA, elemType, idx);       
  int ***Aijk = (int ***)m_stencil;
  
  alignElementAndVertices(m_octDA, elemType, idx);       

  PetscScalar *nuarray = (PetscScalar *)m_nuarray;
  for (int k = 0;k < 8;k++) {
    for (int j=0;j<8;j++) {
		double fac1 = nuarray[idx[k]]*fac11;
      out[m_uiDof*idx[k]] +=  (fac1*(Aijk[elemType][k][j]))*in[m_uiDof*idx[j]];
		//if(hangingMask) std::cout << fac1*Aijk[elemType][k][j] << " " ;
    }//end for j
	 //	 if(hangingMask) std::cout << std::endl;
  }//end for k
  //  if(hangingMask) std::cout << " /* ********************************************************************** */" << std::endl;
  return true;
}

bool stiffnessMatrix::ElementalMatVec(int i, int j, int k, PetscScalar ***in, PetscScalar ***out, double scale){
  int dof= m_uiDof;
  int idx[8][3]={
    {k, j, dof*i},
    {k,j,dof*(i+1)},
    {k,j+1,dof*i},
    {k,j+1,dof*(i+1)},
    {k+1, j, dof*i},
    {k+1,j,dof*(i+1)},
    {k+1,j+1,dof*i},
    {k+1,j+1,dof*(i+1)}               
  };             

  double stencilScale;
  int **Ajk = (int **)m_stencil;

  PetscScalar ***nuarray = (PetscScalar ***)m_nuarray;
  
  for (int q = 0; q < 8; q++) {
    for (int r = 0; r < 8; r++) {
      stencilScale = -(nuarray[k][j][i]*m_dHx)*scale;
		// stencilScale = -(m_nuval*m_dHx)*scale;
		//		stencilScale = -(nuarray[k][j][dof*i]*m_dHx)*scale;
		//stencilScale = -(m_dHx)*scale;
      out[idx[q][0]][idx[q][1]][idx[q][2]] +=  stencilScale*Ajk[q][r]*in[idx[r][0]][idx[r][1]][idx[r][2]];
      if (dof > 1)	
        out[idx[q][0]][idx[q][1]][idx[q][2]+1] += 0.0;
    }
  }
#ifdef __DEBUG__
  //		PetscPrintf(0,"stencil scale in stiffness matrix = %f\n",stencilScale);
#endif

  return true;
}

bool stiffnessMatrix::postMatVec() {
  if ( m_daType == PETSC) {
	 PetscScalar ***nuarray = (PetscScalar ***)m_nuarray;
    int ierr = DAVecRestoreArray(m_DA, nuvec, &nuarray);
    CHKERRQ(ierr);
  } else {
    PetscScalar *nuarray = (PetscScalar *)m_nuarray;
    m_octDA->vecRestoreBuffer(nuvec, nuarray, false, false, true,m_uiDof);
  }

  return true;
}

bool stiffnessMatrix::ElementalMatGetDiagonal(int i, int j, int k, PetscScalar ***diag, double scale) {
  int dof= m_uiDof;
  int idx[8][3]={
    {k, j, dof*i},
    {k,j,dof*(i+1)},
    {k,j+1,dof*i},
    {k,j+1,dof*(i+1)},
    {k+1, j, dof*i},
    {k+1,j,dof*(i+1)},
    {k+1,j+1,dof*i},
    {k+1,j+1,dof*(i+1)}               
  };             

  double stencilScale;
  int **Ajk = (int **)m_stencil;

  PetscScalar ***nuarray = (PetscScalar ***)m_nuarray;
  
  for (int q = 0; q < 8; q++) {
	// stencilScale = -(m_nuval*m_dHx)*scale;
    stencilScale = -(nuarray[k][j][i]*m_dHx)*scale;
    diag [idx[q][0]][idx[q][1]][idx[q][2]] +=  stencilScale*Ajk[q][q] ;
	if (dof > 1)	
      diag[idx[q][0]][idx[q][1]][idx[q][2]+1] += 0.0;
  }

  // std::cout << "Nuval is " << nuarray[k][j][i] << std::endl;
  return true;
}

bool stiffnessMatrix::ElementalMatGetDiagonal(unsigned int i, PetscScalar *diag, double scale) {
  unsigned int lev = m_octDA->getLevel(i);
  double hx = xFac*(1u << (maxD - lev));
  double hy = yFac*(1u << (maxD - lev));
  double hz = zFac*(1u << (maxD - lev));

  double fac11 = -hx*scale/192.0;
  
  stdElemType elemType;
  unsigned int idx[8];

  unsigned char hangingMask = m_octDA->getHangingNodeIndex(i);    //  alignElementAndVertices(m_octDA, elemType, idx);       
  int ***Aijk = (int ***)m_stencil;
  
  alignElementAndVertices(m_octDA, elemType, idx);       

  PetscScalar *nuarray = (PetscScalar *)m_nuarray;
  for (int k = 0;k < 8;k++) {
    double fac1 = nuarray[idx[k]]*fac11;
    diag[m_uiDof*idx[k]] +=  (fac1*(Aijk[elemType][k][k]));
  }//end for k

  return true;
}


#endif /*_STIFFNESSMATRIX_H_*/

