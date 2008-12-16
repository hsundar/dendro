
/**
 * @file odaUtils.h
 * @author		Rahul S. Sampath, rahul.sampath@gmail.com 
 * @brief A list of non-member functions for the ot::DA class.
**/ 

#ifndef __ODAUTILS_H__
#define __ODAUTILS_H__

#include "mpi.h"
#include <vector>

#ifdef __DEBUG__
#ifndef __DEBUG_DA__
#define __DEBUG_DA__
#endif
#endif

#ifdef PETSC_USE_LOG

#include "petscsys.h"

namespace ot {
  extern int pickGhostsEvent;
  extern int addBdySiblingsEvent;
  extern int DAaprioriCommEvent;

  extern int daInitEvent;
  extern int daFinalEvent;
  extern int DAbPart1Event;
  extern int DAbPart2Event;
  extern int DAbPart3Event;
}

#define PROF_DA_INIT_BEGIN PetscLogEventBegin(daInitEvent,0,0,0,0);
#define PROF_DA_INIT_END PetscLogEventEnd(daInitEvent,0,0,0,0); \
  return ;

#define PROF_DA_FINAL_BEGIN PetscLogEventBegin(daFinalEvent,0,0,0,0);
#define PROF_DA_FINAL_END PetscLogEventEnd(daFinalEvent,0,0,0,0); \
  return ;

#define PROF_PICK_GHOSTS_BEGIN PetscLogEventBegin(pickGhostsEvent,0,0,0,0);
#define PROF_PICK_GHOSTS_END PetscLogEventEnd(pickGhostsEvent,0,0,0,0); \
  return ;

#define PROF_ADD_BDY_SIBLINGS_BEGIN PetscLogEventBegin(addBdySiblingsEvent,0,0,0,0);
#define PROF_ADD_BDY_SIBLINGS_END PetscLogEventEnd(addBdySiblingsEvent,0,0,0,0); \
  return ;

#define PROF_DA_APRIORI_COMM_BEGIN PetscLogEventBegin(DAaprioriCommEvent,0,0,0,0);
#define PROF_DA_APRIORI_COMM_END PetscLogEventEnd(DAaprioriCommEvent,0,0,0,0); \
  return ;

#define PROF_DA_BPART1_BEGIN PetscLogEventBegin(DAbPart1Event,0,0,0,0);
#define PROF_DA_BPART1_END PetscLogEventEnd(DAbPart1Event,0,0,0,0); 

#define PROF_DA_BPART2_BEGIN PetscLogEventBegin(DAbPart2Event,0,0,0,0);
#define PROF_DA_BPART2_END PetscLogEventEnd(DAbPart2Event,0,0,0,0); 

#define PROF_DA_BPART3_BEGIN PetscLogEventBegin(DAbPart3Event,0,0,0,0);
#define PROF_DA_BPART3_END PetscLogEventEnd(DAbPart3Event,0,0,0,0); 

#else

#define PROF_DA_INIT_BEGIN 
#define PROF_DA_INIT_END return ;

#define PROF_DA_FINAL_BEGIN 
#define PROF_DA_FINAL_END return ;

#define PROF_ADD_BDY_SIBLINGS_BEGIN 
#define PROF_ADD_BDY_SIBLINGS_END return ;

#define PROF_PICK_GHOSTS_BEGIN 
#define PROF_PICK_GHOSTS_END return;

#define PROF_DA_APRIORI_COMM_BEGIN
#define PROF_DA_APRIORI_COMM_END return;

#define PROF_DA_BPART1_BEGIN 
#define PROF_DA_BPART1_END  

#define PROF_DA_BPART2_BEGIN 
#define PROF_DA_BPART2_END  

#define PROF_DA_BPART3_BEGIN 
#define PROF_DA_BPART3_END  

#endif

namespace ot { 

  namespace DA_FLAGS {

    /**
      @brief A flag to identify the relative order of the vertices of an octant
      @author Hari Sundar
      */	
    enum SortOrderType {
      XYZ, YXZ, XZY, ZXY, YZX, ZYX
    };

    /**           
      @brief A flag to identify the type of octant
      @author Hari Sundar
      */
    enum ElementFlagType {
      DEP_ELEM=32, FOREIGN=255
    };

    /**
      @author Rahul Sampath
      @brief A flag to determine the type of iterator used.

      The set of elements a processor owns will be referred to as local elements. 
      ot::DA allows each processor to loop over its local elements and also the 
      set of ghost octants recieved from processors with a lower rank than itself (pre-ghosts).
      Since, the octants are globally sorted to begin with, the pre-ghosts will
      have a lower Morton id compared to the first local element. 
      The elements (both ghosts and local) will be traversed in the Morton order.

      Each element has 8 node indices stored in it. Some of them could be pointing to ghosts. 

      use ALL to loop over all the local elements and pre-ghost elements as well.

      use WRITABLE to loop over all the elements the current processor owns.

      use DEPENDENT to loop over pre-ghosts as well as local elements, which have
      atleast one node index pointing to some ghost octant. 

      use INDEPENDENT to loop over local elements, which do not 
      have any node index pointing to a ghost octant.

      use W_DEPENDENT to loop over local elements, which have atleast one
      node index pointing to some ghost octant.

      use FROM_STORED to initialize a loop at the location of
      the last call to the function currWithInfo(). Note unlike,
      the other loopTypes FROM_STORED must be used only in the init portion. 

      A sample usage would look like:
     * @code 
     * for ( init<loopType> ; curr() < end<loopType>(); next<loopType>() ) { 
     *     // process current element ...
     * } 
     * @endcode 
     */
    enum loopType {
      ALL, WRITABLE, DEPENDENT, INDEPENDENT, W_DEPENDENT, FROM_STORED
    };

  }//end namespace

  //Forward Declarations
  class TreeNode;
  class DA;

  /**
    @brief A function that maps a given child number and hanging mask to one of the 18 standard hanging configurations.
    @author Rahul Sampath
    @param hnMask The hanging mask returned using the function getHangingNodeIndex()
    @see getHangingNodeIndex
    @see getChildNumber
    */
  template<unsigned char cNum>
    unsigned char getElemType(unsigned char hnMask);

  /**
    @author Rahul Sampath
    @author Hari Sundar
    @brief A function to determine the sort order, i.e.
    the relative orders of the
    indices of the 8 vertices of the current octant
    */
  unsigned int getSortOrder(unsigned int x, unsigned int y, 
      unsigned int z, unsigned int sz);

  /**
    @author Rahul Sampath
    @param curr The current pre-ghost octant
    @param next The next pre-ghost octant
    @brief A helper function required for compressing/uncompressing
    Pre-ghost offsets
    @return a flag describing how the two octants touch each 
    other (if at all they touch)
    */
  unsigned char getTouchConfig(const ot::TreeNode& curr,
      const ot::TreeNode& next, unsigned int maxDepth);

  /**
    @author Rahul Sampath
    @brief Creates a nodal, non-ghosted, single dof vector and set the boundary flag for each node. 
    The flag is one of the enumerations in BoundaryType3
    @see TreeNode::BoundaryType3        
    */
  void assignBoundaryFlags(ot::DA* da, 
      std::vector<unsigned char> & bdyFlagVec);

  /**
    @author Rahul Sampath	
    @return 'true' if the octree is a regular grid
    */
  bool isRegularGrid(ot::DA* da);

  /**
    @author Rahul Sampath
    @brief Interpolates the function and (optionally) its gradient at the specified points
    @param da The octree mesh
    @param in input values 
    @param out output values
    @param gradOut Pass NULL if gradient is not required, otherwise
    pass the address of the desired Vec object. gradOut will be stored
    such that all the components of a node are contiguous. gradOut will be multidimensional
    with dof_of_gradOut = dof_of_in*3. The first 3 components of gradOut correspond to
    the gradient of the first component of in, the next 3 correspond to
    the gradient of the next component of in and so on.
    @param dof Degrees of freedom per node
    @param pts Points to evaluate the field 
    */
  void interpolateData(ot::DA* da, const std::vector<double>& in,
      std::vector<double>& out, std::vector<double>* gradOut,
      unsigned int dof, std::vector<double>& pts);

  /**
    @author Rahul Sampath
    @return the minimum level in the octree
    */
  unsigned int getGlobalMinLevel(ot::DA* da);

  /**
    @author Rahul Sampath
    @return the maximum level in the octree 
    */
  unsigned int getGlobalMaxLevel(ot::DA* da);	  

  /**
    @author Rahul Sampath
    @brief saves the partition in VTK format
    */
  void writePartitionVTK(ot::DA* da, const char* outFilename);

  //@deprecated
  void pickGhostCandidates(const std::vector<ot::TreeNode> & blocks,
      const std::vector<ot::TreeNode> &nodes, std::vector<ot::TreeNode>& res,
      unsigned int dim, unsigned int maxDepth) ;

  void includeSiblingsOfBoundary(std::vector<ot::TreeNode>& allBoundaryLeaves, 
      const ot::TreeNode& myFirstOctant, const ot::TreeNode& myLastOctant);

  void prepareAprioriCommMessagesInDAtype1(const std::vector<ot::TreeNode>& in,
      std::vector<ot::TreeNode>& allBoundaryLeaves, std::vector<ot::TreeNode>& blocks,
      const std::vector<ot::TreeNode>& allBlocks, int myRank, int npes, int* sendCnt,
      std::vector<std::vector<unsigned int> >& sendNodes);

  void prepareAprioriCommMessagesInDAtype2(const std::vector<ot::TreeNode>& in,
      std::vector<ot::TreeNode>& allBoundaryLeaves, std::vector<ot::TreeNode>& blocks,
      const std::vector<ot::TreeNode>& minsOfBlocks, int myRank, int npes, int* sendCnt,
      std::vector<std::vector<unsigned int> >& sendNodes);


  //@deprecated
  int addSecondRing(const std::vector<ot::TreeNode> & nodes,
      const std::vector<unsigned int> & firstLayer,
      const std::vector<ot::TreeNode> & blocks,
      std::vector<unsigned int> & secondRing);

  int DA_blockPartStage2(std::vector<TreeNode> &nodes,
      std::vector<TreeNode> &globalCoarse,
      unsigned int dim, unsigned int maxDepth, MPI_Comm commActive);

  int DA_blockPartStage3(std::vector<TreeNode> &nodes,
      std::vector<TreeNode>& globalCoarse, std::vector<ot::TreeNode>& minsAllBlocks,
      unsigned int dim, unsigned int maxDepth, MPI_Comm commActive);

  /**@brief Initializes the stencils used in the oda module */
  void DA_Initialize(MPI_Comm comm);

  /**@brief Destroys the stencils used in R and P */
  void DA_Finalize();

  int createShapeFnCoeffs_Type1(MPI_Comm comm);
  int createShapeFnCoeffs_Type2(MPI_Comm comm);
  int createShapeFnCoeffs_Type3(MPI_Comm comm);

}//end namespace

/**  
  @def GET_ETYPE_BLOCK
  @author Rahul Sampath	
  @brief A macro used to map the childnumber and hanging mask
  of an octant to one of the 18 hanging configurations.

  @see getChildNumber()
  @see getHangingNodeIndex()
  */
#define GET_ETYPE_BLOCK(type,hnMask,cNum) {\
  if(hnMask) {\
    switch(cNum) {\
      case 0:{\
               type = ot::getElemType<0>(hnMask);\
               break;\
             }\
      case 1:{\
               type = ot::getElemType<1>(hnMask);\
               break;\
             }\
      case 2:{\
               type = ot::getElemType<2>(hnMask);\
               break;\
             }\
      case 3:{\
               type = ot::getElemType<3>(hnMask);\
               break;\
             }\
      case 4:{\
               type = ot::getElemType<4>(hnMask);\
               break;\
             }\
      case 5:{\
               type = ot::getElemType<5>(hnMask);\
               break;\
             }\
      case 6:{\
               type = ot::getElemType<6>(hnMask);\
               break;\
             }\
      case 7:{\
               type = ot::getElemType<7>(hnMask);\
               break;\
             }\
      default:{\
                assert(false);\
              }\
    }\
  }\
}

#include "odaUtils.txx"

#endif /*ODAUTILS_H_*/

