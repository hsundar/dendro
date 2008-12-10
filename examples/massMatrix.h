
/**
*  @file	massMatrix.h
*  @brief	Main class for finite element assembly of a generic
*           mass matrix.
*  @author	Hari Sundar
*  @date	5/8/07
* 
*  Main class for finite element assembly of a generic mass matrix. The density 
*  image needs to be set and is used for the computation of the matrix. MatVec 
*  functions are virtual and in most cases need not be specified as they are 
*  obtained from the parent feMatrix class. In most cases only the assembly of 
*  the element matrix needs to be done within this or a derived class.
*/

#ifndef _MASSMATRIX_H_
#define _MASSMATRIX_H_

/**
*  @brief	Main class for finite element assembly of a generic mass matrix.
*  @author	Hari Sundar
*  @date	5/8/07
*
*  Main class for finite element assembly of a generic mass matrix. The density 
*  image needs to be set and is used for the computation of the matrix. MatVec 
*  functions are virtual and in most cases need not be specified as they are 
*  obtained from the parent feMatrix class. In most cases only the assembly of 
*  the element matrix needs to be done within this or a derived class.
*/

#include "feMatrix.h"

class massMatrix : public ot::fem::feMatrix<massMatrix>
{
  public:
  massMatrix(daType da);
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
    PetscScalar***      nuarray; /* Diffusion coefficient array*/
    Vec                 nuvec;     

    double 		m_dHx;
    
    double xFac, yFac, zFac;
    unsigned int maxD;
};


massMatrix::massMatrix(daType da) {
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

bool massMatrix::initStencils() {
  typedef int* int1Ptr;
  typedef int** int2Ptr;
 
  if (m_daType == PETSC) {
    int Bjk[8][8] =    {
      { 64, 32, 32, 16, 32, 16, 16,  8},
      { 32, 64, 16, 32, 16, 32,  8, 16},
      { 32, 16, 64, 32, 16,  8, 32, 16},
      { 16, 32, 32, 64,  8, 16, 16, 32},
      { 32, 16, 16,  8, 64, 32, 32, 16},
      { 16, 32,  8, 16, 32, 64, 16, 32},
      { 16,  8, 32, 16, 32, 16, 64, 32},
      {  8, 16, 16, 32, 16, 32, 32, 64}
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
      //Type-0:No Hanging
      {{ 64, 32, 32, 16, 32, 16, 16,  8},
        { 32, 64, 16, 32, 16, 32,  8, 16},
        { 32, 16, 64, 32, 16,  8, 32, 16},
        { 16, 32, 32, 64,  8, 16, 16, 32},
        { 32, 16, 16,  8, 64, 32, 32, 16},
        { 16, 32,  8, 16, 32, 64, 16, 32},
        { 16,  8, 32, 16, 32, 16, 64, 32},
        {  8, 16, 16, 32, 16, 32, 32, 64}},

      //Type-1: Y Hanging
      {{ 112,  40,  32,  32,  40,  20,  32,  16},
        {  40,  64,   8,  32,  16,  32,   8,  16},
        {  32,   8,  16,  16,   8,   4,  16,   8},
        {  32,  32,  16,  64,   8,  16,  16,  32},
        {  40,  16,   8,   8,  64,  32,  32,  16},
        {  20,  32,   4,  16,  32,  64,  16,  32},
        {  32,   8,  16,  16,  32,  16,  64,  32},
        {  16,  16,   8,  32,  16,  32,  32,  64}},

      //Type-2: X and Y Hanging
      {{ 168,  36,  36,  48,  48,  36,  36,  24},
        {  36,  16,   4,  16,   8,  16,   4,   8},
        {  36,   4,  16,  16,   8,   4,  16,   8},
        {  48,  16,  16,  64,   8,  16,  16,  32},
        {  48,   8,   8,   8,  64,  32,  32,  16},
        {  36,  16,   4,  16,  32,  64,  16,  32},
        {  36,   4,  16,  16,  32,  16,  64,  32},
        {  24,   8,   8,  32,  16,  32,  32,  64}},

      //Type-3: X and Y and Z Hanging
      {{ 232,  40,  40,  52,  40,  52,  52,  32},
        {  40,  16,   4,  16,   4,  16,   4,   8},
        {  40,   4,  16,  16,   4,   4,  16,   8},
        {  52,  16,  16,  64,   4,  16,  16,  32},
        {  40,   4,   4,   4,  16,  16,  16,   8},
        {  52,  16,   4,  16,  16,  64,  16,  32},
        {  52,   4,  16,  16,  16,  16,  64,  32},
        {  32,   8,   8,  32,   8,  32,  32,  64}},

      //Type-4:XY and X and Y Hanging
      {{ 196,  56,  56,  16,  50,  40,  40,  32},
        {  56,  28,  16,   8,  10,  20,   8,  16},
        {  56,  16,  28,   8,  10,   8,  20,  16},
        {  16,   8,   8,   4,   2,   4,   4,   8},
        {  50,  10,  10,   2,  64,  32,  32,  16},
        {  40,  20,   8,   4,  32,  64,  16,  32},
        {  40,   8,  20,   4,  32,  16,  64,  32},
        {  32,  16,  16,   8,  16,  32,  32,  64}},

      //Type-5:XY and X and Y and Z Hanging
      {{ 262,  61,  61,  17,  41,  56,  56,  40},
        {  61,  28,  16,   8,   5,  20,   8,  16},
        {  61,  16,  28,   8,   5,   8,  20,  16},
        {  17,   8,   8,   4,   1,   4,   4,   8},
        {  41,   5,   5,   1,  16,  16,  16,   8},
        {  56,  20,   8,   4,  16,  64,  16,  32},
        {  56,   8,  20,   4,  16,  16,  64,  32},
        {  40,  16,  16,   8,   8,  32,  32,  64}},

      //Type-6:XY and YZ and X and Y and Z Hanging
      {{ 294,  63,  84,  18,  63,  60,  18,  48},
        {  63,  28,  18,   8,   7,  20,   2,  16},
        {  84,  18,  42,   9,  18,  12,   9,  24},
        {  18,   8,   9,   4,   2,   4,   1,   8},
        {  63,   7,  18,   2,  28,  20,   8,  16},
        {  60,  20,  12,   4,  20,  64,   4,  32},
        {  18,   2,   9,   1,   8,   4,   4,   8},
        {  48,  16,  24,   8,  16,  32,   8,  64}},

      //Type-7: All 6 Hanging
      {{ 328,  87,  87,  19,  87,  19,  19,  56},
        {  87,  42,  21,   9,  21,   9,   3,  24},
        {  87,  21,  42,   9,  21,   3,   9,  24},
        {  19,   9,   9,   4,   3,   1,   1,   8},
        {  87,  21,  21,   3,  42,   9,   9,  24},
        {  19,   9,   3,   1,   9,   4,   1,   8},
        {  19,   3,   9,   1,   9,   1,   4,   8},
        {  56,  24,  24,   8,  24,   8,   8,  64}}
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

bool massMatrix::preMatVec() {
  if (m_daType == PETSC) {
    // compute Hx
    int ierr;
    PetscInt mx,my,mz;
    ierr = DAGetInfo(m_DA,0, &mx, &my, &mz, 0,0,0,0,0,0,0); CHKERRQ(ierr); 

    m_dHx = m_dLx/(mx -1);
    m_dHx = m_dHx*m_dHx*m_dHx;
    m_dHx /= 1728.0;
    CHKERRQ(ierr);
  } else {
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

bool massMatrix::ElementalMatVec(unsigned int i, PetscScalar *in, PetscScalar *out, double scale) {
  unsigned int lev = m_octDA->getLevel(i);
  double hx = xFac*(1u << (maxD - lev));
  double hy = yFac*(1u << (maxD - lev));
  double hz = zFac*(1u << (maxD - lev));

  double fac = scale*hx*hy*hz/1728.0;
  
  stdElemType elemType;
  unsigned int idx[8];
  
  int ***Aijk = (int ***)m_stencil;
  
  alignElementAndVertices(m_octDA, elemType, idx);       
  

  for (int k = 0;k < 8;k++) {
    for (int j=0;j<8;j++) {
      out[m_uiDof*idx[k]] += fac*(Aijk[elemType][k][j])*in[m_uiDof*idx[j]];
		//		std::cout << Aijk[elemType][k][j] << ";" << fac << " ; " << in[m_uiDof*idx[j]] << std::endl;
		//std::cout << xFac << "; " << hx << std::endl;
    }//end for j
  }//end for k
  
  return true;
}
bool massMatrix::ElementalMatVec(int i, int j, int k, PetscScalar ***in, PetscScalar ***out, double scale){
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

  double stencilScale =  m_dHx*scale;
  int **Ajk = (int **)m_stencil;
  
  for (int q = 0; q < 8; q++) {
    for (int r = 0; r < 8; r++) {
      out[idx[q][0]][idx[q][1]][idx[q][2]] += stencilScale*Ajk[q][r]*in[idx[r][0]][idx[r][1]][idx[r][2]];
		if(dof > 1)  out[idx[q][0]][idx[q][1]][idx[q][2]+1] += 0.0;
    }
  }
  return true;
  // std::cout << "Stencil scale is " << stencilScale << std::endl;
}

bool massMatrix::postMatVec() {

  return true;
}

bool massMatrix::ElementalMatGetDiagonal(unsigned int i, PetscScalar *diag, double scale) {
  unsigned int lev = m_octDA->getLevel(i);
  double hx = xFac*(1u << (maxD - lev));
  double hy = yFac*(1u << (maxD - lev));
  double hz = zFac*(1u << (maxD - lev));

  double fac = scale*hx*hy*hz/1728.0;
  
  stdElemType elemType;
  unsigned int idx[8];
  
  int ***Aijk = (int ***)m_stencil;
  
  alignElementAndVertices(m_octDA, elemType, idx);       
  
  for (int k = 0;k < 8;k++) {
      diag[m_uiDof*idx[k]] += fac*(Aijk[elemType][k][k]);
	}//end for k
  
  return true;
}

bool massMatrix::ElementalMatGetDiagonal(int i, int j, int k, PetscScalar ***diag, double scale){
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

  double stencilScale =  m_dHx*scale;
  int **Ajk = (int **)m_stencil;
  
  for (int q = 0; q < 8; q++) {
    diag[idx[q][0]][idx[q][1]][idx[q][2]] += stencilScale*Ajk[q][q];
		if(dof > 1)  diag[idx[q][0]][idx[q][1]][idx[q][2]+1] += 0.0;
  }
  return true;
}


#endif /*_MASSMATRIX_H_*/

