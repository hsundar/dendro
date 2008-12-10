
/** 
 * @file omg.h
 * @brief 		Defines the basic MG class, which provides an octree
 * 			based Multigrid compatible with PETSc.  
 * @author		Rahul S Sampath, rahul.sampath@gmail.com 
 * @date 		2007-04-22
 *
 * Defines the basic MG class, which provides an octree based Multigrid
 * compatible with PETSc.
 **/ 

#ifndef __OCT_MG_H__
#define __OCT_MG_H__

#include "petscpc.h"
#include "petscksp.h"
#include <vector>

#ifdef PETSC_USE_LOG

#include "petscsys.h"

namespace ot {
  extern int setDaEvent;
  extern int setUpEvent;
  extern int createRp1Event;
  extern int createRp2Event;
  extern int setKspEvent;
  extern int restrictEvent;
  extern int dummyRestrictEvent;
  extern int prolongEvent;
  extern int scatterEvent;
  extern int damgInitEvent;
  extern int damgFinalEvent;
  
  extern int setDAstage1Event;
  extern int setDAstage2Event;
  extern int setDAstage3Event;
  extern int setDAstage4Event;
  extern int setDAstage5Event;
  extern int setDAstage6Event;
}

#define PROF_SET_DA_STAGE1_BEGIN PetscLogEventBegin(setDAstage1Event,0,0,0,0);
#define PROF_SET_DA_STAGE1_END PetscLogEventEnd(setDAstage1Event,0,0,0,0); 

#define PROF_SET_DA_STAGE2_BEGIN PetscLogEventBegin(setDAstage2Event,0,0,0,0);
#define PROF_SET_DA_STAGE2_END PetscLogEventEnd(setDAstage2Event,0,0,0,0); 

#define PROF_SET_DA_STAGE3_BEGIN PetscLogEventBegin(setDAstage3Event,0,0,0,0);
#define PROF_SET_DA_STAGE3_END PetscLogEventEnd(setDAstage3Event,0,0,0,0); 

#define PROF_SET_DA_STAGE4_BEGIN PetscLogEventBegin(setDAstage4Event,0,0,0,0);
#define PROF_SET_DA_STAGE4_END PetscLogEventEnd(setDAstage4Event,0,0,0,0); 

#define PROF_SET_DA_STAGE5_BEGIN PetscLogEventBegin(setDAstage5Event,0,0,0,0);
#define PROF_SET_DA_STAGE5_END PetscLogEventEnd(setDAstage5Event,0,0,0,0); 

#define PROF_SET_DA_STAGE6_BEGIN PetscLogEventBegin(setDAstage6Event,0,0,0,0);
#define PROF_SET_DA_STAGE6_END PetscLogEventEnd(setDAstage6Event,0,0,0,0); 

#define PROF_MG_INIT_BEGIN  PetscFunctionBegin; \
  PetscLogEventBegin(damgInitEvent,0,0,0,0);

#define PROF_MG_INIT_END PetscLogEventEnd(damgInitEvent,0,0,0,0); \
  PetscFunctionReturn(0);

#define PROF_MG_FINAL_BEGIN  PetscFunctionBegin; \
  PetscLogEventBegin(damgFinalEvent,0,0,0,0);

#define PROF_MG_FINAL_END PetscLogEventEnd(damgFinalEvent,0,0,0,0); \
  PetscFunctionReturn(0);

#define PROF_MG_SET_KSP_BEGIN  PetscFunctionBegin; \
  PetscLogEventBegin(setKspEvent,0,0,0,0);

#define PROF_MG_SET_KSP_END PetscLogEventEnd(setKspEvent,0,0,0,0); \
  PetscFunctionReturn(0);

#define PROF_MG_SET_DA_BEGIN  PetscFunctionBegin; \
  PetscLogEventBegin(setDaEvent,0,0,0,0);

#define PROF_MG_SET_DA_END PetscLogEventEnd(setDaEvent,0,0,0,0); \
  PetscFunctionReturn(0);

#define PROF_MG_SETUP_BEGIN  PetscFunctionBegin; \
  PetscLogEventBegin(setUpEvent,0,0,0,0);

#define PROF_MG_SETUP_END PetscLogEventEnd(setUpEvent,0,0,0,0); \
  PetscFunctionReturn(0);

#define PROF_MG_CREATE_RP1_BEGIN  PetscFunctionBegin; \
  PetscLogEventBegin(createRp1Event,0,0,0,0);

#define PROF_MG_CREATE_RP1_END PetscLogEventEnd(createRp1Event,0,0,0,0); \
  PetscFunctionReturn(0);

#define PROF_MG_CREATE_RP2_BEGIN  PetscFunctionBegin; \
  PetscLogEventBegin(createRp2Event,0,0,0,0);

#define PROF_MG_CREATE_RP2_END PetscLogEventEnd(createRp2Event,0,0,0,0); \
  PetscFunctionReturn(0);

#define PROF_MG_RESTRICT_BEGIN PetscFunctionBegin; \
  PetscLogEventBegin(restrictEvent,0,0,0,0);

#define PROF_MG_RESTRICT_END PetscLogEventEnd(restrictEvent,0,0,0,0); \
  PetscFunctionReturn(0);

#define PROF_MG_RESTRICT_DUMMY_BEGIN PetscFunctionBegin; \
  PetscLogEventBegin(dummyRestrictEvent,0,0,0,0);

#define PROF_MG_RESTRICT_DUMMY_END PetscLogEventEnd(dummyRestrictEvent,0,0,0,0); \
  PetscFunctionReturn(0);

#define PROF_MG_PROLONG_BEGIN  PetscFunctionBegin; \
  PetscLogEventBegin(prolongEvent,0,0,0,0);

#define PROF_MG_PROLONG_END PetscLogEventEnd(prolongEvent,0,0,0,0); \
  PetscFunctionReturn(0);

#define PROF_MG_SCATTER_BEGIN PetscLogEventBegin(scatterEvent,0,0,0,0);

#define PROF_MG_SCATTER_END PetscLogEventEnd(scatterEvent,0,0,0,0); \
  return 0;

#else

#define PROF_SET_DA_STAGE1_BEGIN 
#define PROF_SET_DA_STAGE1_END 

#define PROF_SET_DA_STAGE2_BEGIN 
#define PROF_SET_DA_STAGE2_END 

#define PROF_SET_DA_STAGE3_BEGIN 
#define PROF_SET_DA_STAGE3_END 

#define PROF_SET_DA_STAGE4_BEGIN 
#define PROF_SET_DA_STAGE4_END 

#define PROF_SET_DA_STAGE5_BEGIN 
#define PROF_SET_DA_STAGE5_END 

#define PROF_SET_DA_STAGE6_BEGIN 
#define PROF_SET_DA_STAGE6_END 

#define PROF_MG_INIT_BEGIN PetscFunctionBegin; 

#define PROF_MG_INIT_END PetscFunctionReturn(0);

#define PROF_MG_FINAL_BEGIN PetscFunctionBegin; 

#define PROF_MG_FINAL_END PetscFunctionReturn(0);

#define PROF_MG_SET_KSP_BEGIN PetscFunctionBegin; 

#define PROF_MG_SET_KSP_END PetscFunctionReturn(0);

#define PROF_MG_SET_DA_BEGIN PetscFunctionBegin; 

#define PROF_MG_SET_DA_END PetscFunctionReturn(0);

#define PROF_MG_SETUP_BEGIN PetscFunctionBegin; 

#define PROF_MG_SETUP_END PetscFunctionReturn(0);

#define PROF_MG_CREATE_RP1_BEGIN PetscFunctionBegin; 

#define PROF_MG_CREATE_RP1_END PetscFunctionReturn(0);

#define PROF_MG_CREATE_RP2_BEGIN PetscFunctionBegin; 

#define PROF_MG_CREATE_RP2_END PetscFunctionReturn(0);

#define PROF_MG_SCATTER_BEGIN 

#define PROF_MG_SCATTER_END return 0;

#define PROF_MG_RESTRICT_BEGIN PetscFunctionBegin; 

#define PROF_MG_RESTRICT_END PetscFunctionReturn(0);

#define PROF_MG_RESTRICT_DUMMY_BEGIN PetscFunctionBegin; 

#define PROF_MG_RESTRICT_DUMMY_END PetscFunctionReturn(0);

#define PROF_MG_PROLONG_BEGIN PetscFunctionBegin; 

#define PROF_MG_PROLONG_END PetscFunctionReturn(0);

#endif

namespace ot {

  /**
    @struct	PC_KSP_Shell
    @brief A private preconditioner object used within DAMG at the coarsest grid
    when not all processors are active on the coarsest grid.
    @author Rahul Sampath
    This must be used with KSPPREONLY only 
    */
  typedef struct {
    bool iAmActive; /**< True if the calling processor is active. */
    MPI_Comm commActive; /**< The active processors */
    PC pc; /**< The PCSHELL itself */
    KSP ksp_private; /**< Internal KSP */
    Vec rhs_private; /**< Internal rhs vector for ksp_private */ 
    Vec sol_private; /**< Internal rhs vector for sol_private */
  } PC_KSP_Shell;

  PetscErrorCode PC_KSP_Shell_SetUp(void* ctx);

  PetscErrorCode PC_KSP_Shell_Apply(void* ctx, Vec in, Vec out);

  PetscErrorCode PC_KSP_Shell_Destroy(void* ctx);

  /**
    @class	FineTouchedStatus
    @brief A private class used inside Restriction and Prolongation.
    @author Rahul Sampath

    The masks for the inter-grid transfer operators are stored in and object of this class.
    */
  class FineTouchedStatus {
    public:
      unsigned char flags[8];

      /** @name Constructors */
      //@{
      FineTouchedStatus() {
	for(int i=0;i<8;i++) {
	  flags[i] = 0;
	}
      }    

      FineTouchedStatus(const FineTouchedStatus & other) {
	for(int i=0; i<8; i++) {
	  this->flags[i] = other.flags[i];
	}
      }
      //@}

      /** @name Overloaded Operators */
      //@{
      FineTouchedStatus & operator = (FineTouchedStatus const  & other) {
	if(this == (&other)) {return *this;}	
	for(int i=0; i<8; i++) {
	  this->flags[i] = other.flags[i];
	}
	return *this;
      }//end fn.

      FineTouchedStatus & operator +=(FineTouchedStatus const & other) {
	//Self-Assignment not a problem!
	for(int i=0; i<8; i++) {
	  this->flags[i] |= other.flags[i];
	}
	return *this;
      }//end fn.

      bool  operator == ( FineTouchedStatus const  &other) const {
	for(int i = 0; i < 8; i++ ) {
	  if(this->flags[i] != other.flags[i]) {
	    return false;
	  }
	}
	return true;
      }

      bool  operator != ( FineTouchedStatus const  &other) const {
	return (!((*this) == other));
      }

      const FineTouchedStatus operator + (FineTouchedStatus const &other) const {
	return ((FineTouchedStatus(*this)) += other);
      }
      //@}

  };

  /**
    @class	FineTouchedDummyStatus
    @brief A class used within the dummy restriction operation
    @author Rahul Sampath
    */
  class FineTouchedDummyStatus {
    public:
      unsigned char flags[16];

      /** @name Constructors */
      //@{
      FineTouchedDummyStatus() {
	for(int i=0;i<16;i++) {
	  flags[i] = 0;
	}
      }    

      FineTouchedDummyStatus(const FineTouchedDummyStatus & other) {
	for(int i=0; i<16; i++) {
	  this->flags[i] = other.flags[i];
	}
      }
      //@}

      /** @name Overloaded operators */
      //@{
      FineTouchedDummyStatus & operator = (FineTouchedDummyStatus const  & other) {
	if(this == (&other)) {return *this;}	
	for(int i=0; i<16; i++) {
	  this->flags[i] = other.flags[i];
	}
	return *this;
      }//end fn.

      FineTouchedDummyStatus & operator +=(FineTouchedDummyStatus const & other) {
	//Self-Assignment not a problem!
	for(int i=0; i<16; i++) {
	  this->flags[i] |= other.flags[i];
	}
	return *this;
      }//end fn.

      bool  operator == ( FineTouchedDummyStatus const  &other) const {
	for(int i = 0; i < 16; i++ ) {
	  if(this->flags[i] != other.flags[i]) {
	    return false;
	  }
	}
	return true;
      }

      bool  operator != ( FineTouchedDummyStatus const  &other) const {
	return (!((*this) == other));
      }

      const FineTouchedDummyStatus operator + (FineTouchedDummyStatus const &other) const {
	return ((FineTouchedDummyStatus(*this)) += other);
      }
      //@}

  };

  //Forward Declarations
  class DA;
  class TreeNode;

  /**
    @struct	TransferOpData
    @brief Context used for storing information required by the inter-grid transfer operators
    @author Rahul Sampath
    */
  typedef struct {
    ot::DA* dac; /**< Coarse mesh */
    ot::DA* daf; /**< Fine mesh */
    unsigned int minIndependentSize; /**< Min (across processors) number of independent
				       elements on the coarse grid. This us used for
				       overlapping communication and computation inside
				       the R and P matvecs.  */
    unsigned char* suppressedDOFc; /**< Dirichlet nodes on the coarse mesh */
    unsigned char* suppressedDOFf; /**< Dirichlet nodes on the fine mesh */
    std::vector<ot::FineTouchedStatus >* fineTouchedFlags; /**< The masks used for Restriction/Prolongation */
    unsigned int dof; /**< The number of degrees of freedom per node */
    MPI_Comm comm;
    Vec tmp; //For R/P-type2 scatter
    Vec addRtmp;
    Vec addPtmp;
    /** @name Communication Stuff for Scatters */
    //@{
    int * sendSzP; 
    int * sendOffP; 
    int * recvSzP; 
    int * recvOffP;
    int * sendSzR; 
    int * sendOffR; 
    int * recvSzR; 
    int * recvOffR;
    //@}
  } TransferOpData;

  /**
    @struct	_p_DAMG
    @brief The Octree-Multigrid Object
    @author Rahul Sampath
    */
  struct _p_DAMG {

    ot::DA* da; /**< octree mesh used for smoothing at this level */
    ot::DA* da_aux;  /**< Pseudo-mesh used in inter-grid transfers (Scatters). This is not used for smoothing. 
		       This is NULL, if the immediate coarser grid shares the same partition. */

    unsigned int dof; /**< the number of degrees of freedom per node */

    unsigned char* suppressedDOF; /**< A list of Dirichlet boundary nodes on 'da'. Must be the size of (dof*localbufferSize) */
    unsigned char* suppressedDOFaux; /**< A list of Dirichlet boundary nodes on 'da_aux'. Must be the size of (dof*localbufferSize) */ 

    Vec            x,b,r;                /**< global vectors used in multigrid preconditioner for this level*/
    Mat            J;                    /**< matrix on this level */
    Mat            B;			/**< preconditioning matrix for this level */	
    Mat            R;                    /**< interpolation to next finer level */
    int       nlevels;              /**< number of levels above this one (total number of levels on level 0)*/
    int       totalLevels; /**< The total number of levels. */
    MPI_Comm       comm; /**< The communicator */

    /**
      The pointer to the solve function is just taken from the PETSc implementation.
      In PETSc, it is used to choose between a KSP solve and SNES solve or FAS solve.
      However, OTK only supports KSP solve. This can be extended later if required. 
      */
    PetscErrorCode (*solve)(_p_DAMG**, int);

    void           *user;        /**< User context */ 

    KSP            ksp;  /**< The solver */           
    PetscErrorCode (*initialguess)(_p_DAMG*, Vec); /**< Function handle to compute the initial guess vector */
    PetscErrorCode (*rhs)(_p_DAMG*,Vec); /**< Function handle to compute the RHS vector */
    PetscTruth     matricesset;  /**< true if user had called stsDMMGSetKSP() and the matrices have been computed */
  };//end struct definition

  /** The multigrid object */
  typedef struct _p_DAMG* DAMG; 

  //Public Functions

  /** 
    @brief Constructs the Multigrid object.
    @param comm The communicator
    @param nlevels maximum number of multigrid levels. This may be reset to something smaller within the function.
    @param user User context
    @param damg The multigrid object
    @param finestOctree A 2:1 balanced octree for the finest level. This will be destroyed within the function
    @param dof Degrees of freedom per node
    @param loadFac Tolerance on load imbalance at each level
    @param compressLut Option to compress the element-to-node mappings at each level
    @param incCorner Option to balance across corners as well
    @return error flag
    */
  PetscErrorCode DAMGCreateAndSetDA(MPI_Comm comm, int& nlevels,
      void* user, DAMG** damg, std::vector<ot::TreeNode>& finestOctree,
      unsigned int dof = 1, double loadFac = 1.5,
      bool compressLut = true, bool incCorner = true);

  /**
    @brief This is the function used to create the J matrix. This must be called before calling DAMGSetKSP if you want to set
    J and B to be different. If this is not called, by default J and B will be equal.
    */
  PetscErrorCode DAMGCreateJMatrix(DAMG damg, PetscErrorCode (*crJ)(DAMG, Mat*J));

  /**
    @brief Set function handles to create and set values in the matrices and compute the RHS 
    @param damg the multigrid object
    @param crjac function handle to the function that allocates memory for the precondition matrix at each level
    @param computeJac function handle to the function that sets matrix entries/ setups the matrices at each level.
    @param rhs function handle to compute the RHS vector
    */
  PetscErrorCode DAMGSetKSP(DAMG* damg, PetscErrorCode (*crjac)(DAMG, Mat* B),
      PetscErrorCode (*computeJac)(DAMG, Mat J, Mat B), PetscErrorCode (*rhs)(DAMG, Vec));

  /**@brief Solves the problem */
  PetscErrorCode DAMGSolve(DAMG* damg);

  /**@brief Set a function handle that will be used to generate the initial guess. */
  PetscErrorCode DAMGSetInitialGuess(DAMG*, PetscErrorCode (*)(DAMG, Vec));

  /**@brief Use the current vector as the initial guess. */
  PetscErrorCode DAMGInitialGuessCurrent(DAMG, Vec);

  /**@brief Prints detailed information about the meshes for each level */
  void PrintDAMG(DAMG*);

  /**@brief Call this function to allocate memory for vectors used to mark the dirichlet nodes.*/
  int DAMGCreateSuppressedDOFs(DAMG* damg);

  /**@brief Loads the stencils used in R and P. Processor 0 reads the files and
  ** broadcasts to others.
  */
  void DAMG_InitPrivateType1(MPI_Comm comm);

  /**@brief Loads the stencils used in R and P. The comm is split into many
  ** groups (each group containing a maximum of 1000 processors) and the
  first processor in each group reads the file and broadcasts
  ** the stencils to the other processors in the group*/
  void DAMG_InitPrivateType2(MPI_Comm comm);

  /**@brief Loads the stencils used in R and P. No communication. Each
  ** processor opens a unique file.  */
  void DAMG_InitPrivateType3(MPI_Comm comm);

  /**@brief Initializes the stencils used in R and P */
  PetscErrorCode DAMG_Initialize(MPI_Comm comm);

  /**@brief Destroys the stencils used in R and P */
  PetscErrorCode DAMG_Finalize();

  /**@brief Destroy the DAMG object */
  PetscErrorCode DAMGDestroy(DAMG* damg);

  /**
    @def DAMGGetRHS
    @brief The RHS vector at the finest level 
    */
#define DAMGGetRHS(ctx)              (ctx)[(ctx)[0]->nlevels-1]->b

  /**
    @def DAMGGetr
    @brief The residual vector at the finest level
    */
#define DAMGGetr(ctx)              (ctx)[(ctx)[0]->nlevels-1]->r

  /**
    @def DAMGGetx
    @brief The finest level Vec object containing the solution
    */
#define DAMGGetx(ctx)              (ctx)[(ctx)[0]->nlevels-1]->x

  /**
    @def DAMGGetComm
    @brief Returns the MPI communicator 
    */
#define DAMGGetComm(ctx)           (ctx)[(ctx)[0]->nlevels-1]->comm

  /**
    @def DAMGGetJ
    @brief The finest level matrix
    */
#define DAMGGetJ(ctx)              (ctx)[(ctx)[0]->nlevels-1]->J

  /** 
    @def DAMGGetB
    @brief The finest level preconditioning matrix 
    */
#define DAMGGetB(ctx)              (ctx)[(ctx)[0]->nlevels-1]->B

  /**
    @def DAMGGetFine
    @brief The finest DAMG object 
    */
#define DAMGGetFine(ctx)           (ctx)[(ctx)[0]->nlevels-1]

  /**
    @def DAMGGetKSP
    @brief The finest grid KSP object 
    */
#define DAMGGetKSP(ctx)            (ctx)[(ctx)[0]->nlevels-1]->ksp

  /** 
    @def DAMGGetDA
    @brief The finest grid ot::DA* object 
    */
#define DAMGGetDA(ctx)             ((ctx)[(ctx)[0]->nlevels-1]->da)

  /** 
    @def DAMGGetUser
    @brief The user context at this level 
    */
#define DAMGGetUser(ctx,level)     ((ctx)[level]->user)

  /**
    @def DAMGSetUser
    @brief Set the user context for this level 
    */
#define DAMGSetUser(ctx,level,usr) ((ctx)[level]->user = usr,0)

  /** 
    @def DAMGGetLevels
    @brief The total number of levels 
    */
#define DAMGGetLevels(ctx)         (ctx)[0]->nlevels

  /**
    @def DAMGGetDAMG
    @brief The finest DAMG. 
    */
#define DAMGGetDAMG(ctx)              (ctx)[(ctx)[0]->nlevels-1]

  PetscErrorCode DAMGSetNullSpace(DAMG* damg, PetscTruth, int, PetscErrorCode (*)(DAMG,Vec[]));

  /** @name Private Functions */
  //@{
  /**
    @author Rahul Sampath	
    @brief Redistributes the vector, preserving the relative ordering of the elements
    @param in the input vector
    @param out the output vector
    @param inSz local size of in
    @param outSz local size of out
    Memory for out will not be allocated within the function.
    */
  int scatterValues(Vec in, Vec out, PetscInt inSz, PetscInt outSz,
      int *& sendSz, int *& sendOff, int *& recvSz, int *& recvOff, MPI_Comm comm);

  PetscErrorCode DAMGSetUp(DAMG*);

  PetscErrorCode DAMGSolveKSP(DAMG *damg, int level);

  PetscErrorCode DAMGSetUpLevel(DAMG* damg, KSP ksp, int nlevels);

  /*Matrix-Free Intergrid Transfer Operators */
  int destroyRmatType1Stencil(double *****&lut);
  int destroyRmatType2Stencil(double ****&lut);
  int destroyVtxMaps(unsigned short ****&map1, unsigned short *****&map2,
      unsigned short *****&map3, unsigned short ******&map4);

  int readRmatType1Stencil(double *****&lut);
  int readRmatType2Stencil(double ****&lut);
  int readVtxMaps(unsigned short ****&map1, unsigned short *****&map2,
      unsigned short *****&map3, unsigned short ******&map4);

  int IreadRmatType1Stencil(double *****&lut, int rank);
  int IreadRmatType2Stencil(double ****&lut, int rank);
  int IreadVtxMaps(unsigned short ****&map1, unsigned short *****&map2,
      unsigned short *****&map3, unsigned short ******&map4, int rank);


  //Coarse and fine are aligned. No need to use da_aux
  PetscErrorCode createInterpolationType1(DAMG , DAMG, Mat *);
  PetscErrorCode prolongMatVecType1(Mat, Vec, Vec);
  PetscErrorCode restrictMatVecType1(Mat, Vec, Vec);
  PetscErrorCode dummyRestrictMatVecType1(TransferOpData* data);

  //Coarse and fine are NOT aligned. Need to use da_aux.
  PetscErrorCode  createInterpolationType2(DAMG , DAMG, Mat *);
  PetscErrorCode prolongMatVecType2(Mat, Vec, Vec);
  PetscErrorCode   restrictMatVecType2(Mat, Vec, Vec);

  PetscErrorCode addProlongMatVec(Mat, Vec, Vec, Vec);
  PetscErrorCode  addRestrictMatVec(Mat, Vec, Vec, Vec);

  PetscErrorCode  rpMatDestroy(Mat);

  //@}

}//end namespace


namespace par {

  //Forward Declaration
  template <typename T>
    class Mpi_datatype;

  //Used only with ReadFromGhosts
  template <>
    class Mpi_datatype< ot::FineTouchedStatus > {
      public: 
	static MPI_Datatype value()
	{
	  static bool  first = true;
	  static MPI_Datatype datatype;

	  if (first)
	  {
	    first = false;
	    MPI_Type_contiguous(sizeof(ot::FineTouchedStatus), MPI_BYTE, &datatype);
	    MPI_Type_commit(&datatype);
	  }

	  return datatype;
	}
    };

  //Used only with WriteToGhosts
  template <>
    class Mpi_datatype< ot::FineTouchedDummyStatus > {
      public: 
	static MPI_Datatype value()
	{
	  static bool  first = true;
	  static MPI_Datatype datatype;

	  if (first)
	  {
	    first = false;
	    MPI_Type_contiguous(sizeof(ot::FineTouchedDummyStatus), MPI_BYTE, &datatype);
	    MPI_Type_commit(&datatype);
	  }

	  return datatype;
	}
    };

}//end namespace par


#endif

