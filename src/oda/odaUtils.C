
/**
  @file odaUtils.C
  @brief A collection of simple functions for supporting octree-mesh related operations.  
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  */

#include "odaUtils.h"
#include "TreeNode.h"
#include "oda.h"
#include "parUtils.h"
#include "seqUtils.h"

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

namespace ot {

  void writePartitionVTK(ot::DA* da, const char* outFileName) {
    int rank = da->getRankAll();
    //Only processor writes
    if(!rank) {
      std::vector<ot::TreeNode> minBlocks = da->getMinAllBlocks();
      unsigned int maxDepth = da->getMaxDepth();
      ot::TreeNode root(3, maxDepth);
      std::vector<ot::TreeNode> allBlocks;
      ot::completeSubtree(root, minBlocks, allBlocks, 3, maxDepth, true, true);

      FILE* outfile = fopen(outFileName,"w");

      //Set the weights of allBlocks to be the processor ids. 
      for(unsigned int i = 0; i < allBlocks.size(); i++) {
        unsigned int pId;
        bool found = seq::maxLowerBound<ot::TreeNode>(minBlocks, allBlocks[i], pId, NULL, NULL);
        assert(found);
        allBlocks[i].setWeight(pId);
      }

      unsigned int numNode = static_cast<unsigned int>(allBlocks.size());

      float coord[8][3] = {
        {0.0,0.0,0.0},
        {1.0,0.0,0.0},
        {0.0,1.0,0.0},
        {1.0,1.0,0.0},
        {0.0,0.0,1.0},
        {1.0,0.0,1.0},
        {0.0,1.0,1.0},
        {1.0,1.0,1.0}
      };

      fprintf(outfile,"# vtk DataFile Version 3.0\n");
      fprintf(outfile,"Octree field file\n");
      fprintf(outfile,"ASCII\n");
      fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");
      fprintf(outfile,"POINTS %d float\n",(numNode*8));

      for (unsigned int i = 0; i < numNode; i++) {
        unsigned int x = allBlocks[i].getX();
        unsigned int y = allBlocks[i].getY(); 
        unsigned int z = allBlocks[i].getZ();
        unsigned int d = allBlocks[i].getLevel();
        float fx, fy, fz,hx;
        fx = ((float)x)/((float)(1u<<maxDepth));
        fy = ((float)y)/((float)(1u<<maxDepth));
        fz = ((float)z)/((float)(1u<<maxDepth));
        hx = ((float)(1u<<(maxDepth-d)))/((float)(1u<<maxDepth));

        for(int j = 0; j < 8; j++)
        {
          float fxn,fyn,fzn;
          fxn = fx + coord[j][0]*hx;
          fyn = fy + coord[j][1]*hx;
          fzn = fz + coord[j][2]*hx;
          fprintf(outfile,"%f %f %f \n",fxn,fyn,fzn);
        }
      }

      fprintf(outfile,"\nCELLS %d %d\n",numNode,numNode*9);

      for(int i = 0; i < numNode; i++)
      {
        fprintf(outfile,"8 ");

        for(int j = 0; j < 8; j++)
        {
          int index = (8*i)+j;
          fprintf(outfile,"%d ",index);
        }
        fprintf(outfile,"\n");
      }

      fprintf(outfile,"\nCELL_TYPES %d\n",numNode);

      for(int i = 0; i < numNode; i++)
      {
        fprintf(outfile,"11 \n");
      }

      fprintf(outfile,"\nCELL_DATA %d\n",numNode);
      fprintf(outfile,"SCALARS scalars unsigned_int\n");
      fprintf(outfile,"LOOKUP_TABLE default\n");

      for (unsigned int i =0; i< numNode; i++) {
        unsigned int v = allBlocks[i].getWeight();
        fprintf(outfile,"%u \n", v);
      }

      fclose(outfile);

    }//end if p0
  }//end function

  unsigned int getGlobalMinLevel(ot::DA* da) {

    unsigned int myMinLev = ot::TreeNode::MAX_LEVEL;
    unsigned int globalMinLev;

    //It is sufficient to loop over the elements, since boundaries were added at the same level as some element. 
    if(da->iAmActive()) {
      for(da->init<ot::DA_FLAGS::WRITABLE>();
          da->curr() < da->end<ot::DA_FLAGS::WRITABLE>();
          da->next<ot::DA_FLAGS::WRITABLE>())  
      {
        unsigned int currLevel = da->getLevel(da->curr());
        if(currLevel < myMinLev) {
          myMinLev = currLevel;
        }
      }      
    }

    par::Mpi_Allreduce<unsigned int>(&myMinLev, &globalMinLev, 1, MPI_MIN, da->getComm() );

    //The initial octree is not root. So the min lev in the initial octree is atleast 1. So in the embedded octree the minlev is atleast 2. 
    assert(globalMinLev > 1);

    //Return the result in the original octree configuration
    return (globalMinLev -1);
  }

  unsigned int getGlobalMaxLevel(ot::DA* da) {

    unsigned int myMaxLev = 0;
    unsigned int globalMaxLev;

    //It is sufficient to loop over the elements, since boundaries were added at the same level as some element. 
    if( da->iAmActive() ) {
      for(da->init<ot::DA_FLAGS::WRITABLE>(); 
          da->curr() < da->end<ot::DA_FLAGS::WRITABLE>();
          da->next<ot::DA_FLAGS::WRITABLE>())  
      {
        unsigned int currLevel = da->getLevel(da->curr());
        if(currLevel > myMaxLev) {
          myMaxLev = currLevel;
        }
      }      
    }

    par::Mpi_Allreduce<unsigned int>(&myMaxLev, &globalMaxLev, 1, MPI_MAX, da->getComm() );

    //m_uiMaxDepth will be 0 on inActive processors.
    if(da->iAmActive()) {
      assert(globalMaxLev <= da->getMaxDepth() );
    }else {
      assert(globalMaxLev <= ot::TreeNode::MAX_LEVEL);
    }

    //Return the result in the original octree configuration
    return (globalMaxLev -1);	
  }

  unsigned int getSortOrder(unsigned int x, unsigned int y, 
      unsigned int z, unsigned int sz) {

    // first compare the x, y, and z to determine which one dominates ...
    unsigned int _x = x^(x+sz);
    unsigned int _y = y^(y+sz);
    unsigned int _z = z^(z+sz);

    // compute order ...
    if (_x > _y) {
      if ( _y > _z) {
        return ot::DA_FLAGS::ZYX;
      } else if ( _x > _z ) {
        return ot::DA_FLAGS::YZX;
      } else {
        return ot::DA_FLAGS::YXZ;
      }
    } else {
      if ( _y > _z) {
        return ot::DA_FLAGS::ZXY;
      } else if ( _x > _z ) {
        return ot::DA_FLAGS::ZXY;
      } else {
        return ot::DA_FLAGS::XYZ;
      }
    }
  }

  unsigned char getTouchConfig(const ot::TreeNode& curr, 
      const ot::TreeNode& next, unsigned int maxDepth) {
    unsigned char c = 0;
    unsigned int cx = curr.minX();
    unsigned int cy = curr.minY();
    unsigned int cz = curr.minZ();
    unsigned int nx = next.minX();
    unsigned int ny = next.minY();
    unsigned int nz = next.minZ();
    unsigned int cd = curr.getLevel();
    unsigned int nd = next.getLevel();
    unsigned int cs = (1u<<(maxDepth - cd));
    unsigned int ns = (1u<<(maxDepth - nd));

    //_zzyyxxT  
    bool isTouchingX = false;
    bool isTouchingY = false;
    bool isTouchingZ = false;

    if (cx == nx) {
      isTouchingX = true;
    }

    if (cx == (nx + ns) ) {
      isTouchingX = true;
      c += (1 << 1);
    }

    if ( (cx + cs) == nx ) {
      isTouchingX = true;
      c += (2 << 1);
    }

    if ( (cx + cs) == (nx + ns) ) {
      isTouchingX = true;
      c += (3 << 1);
    }

    if (cy == ny) {
      isTouchingY = true;
    }

    if (cy == (ny + ns) ) {
      isTouchingY = true;
      c += (1 << 3);
    }

    if ( (cy + cs) == ny ) {
      isTouchingY = true;
      c += (2 << 3);
    }

    if ( (cy + cs) == (ny + ns) ) {
      isTouchingY = true;
      c += (3 << 3);
    }

    if (cz == nz) {
      isTouchingZ = true;
    }

    if (cz == (nz + ns) ) {
      isTouchingZ = true;
      c += (1 << 5);
    }

    if ( (cz + cs) == nz ) {
      isTouchingZ = true;
      c += (2<<5);
    }

    if ( (cz + cs) == (nz + ns) ) {
      isTouchingZ = true;
      c += (3<<5);
    }

    if ( isTouchingX && isTouchingY && isTouchingZ ) {
      c += 1;
    } else {
      c = 0;
    }

    return c;
  }

  bool isRegularGrid(ot::DA* da) {
    int iHaveHanging = 0;
    if(da->iAmActive()) {
      for(da->init<ot::DA_FLAGS::WRITABLE>(); 
          da->curr() < da->end<ot::DA_FLAGS::WRITABLE>();
          da->next<ot::DA_FLAGS::WRITABLE>()) { 
        if(da->isHanging(da->curr())) {
          iHaveHanging = 1;
          std::cout<<(da->getRankActive())<<
            " found a hanging node for element with id: "
            <<(da->curr())<<std::endl;

          break;
        }
      }//end for writable
    }//end if active

    int anyOneHasHanging;
    par::Mpi_Allreduce<int>(&iHaveHanging, &anyOneHasHanging, 1, MPI_SUM, da->getComm());

    return (!anyOneHasHanging);
  }//end function

  void assignBoundaryFlags(ot::DA* da, 
      std::vector<unsigned char> & bdyFlagVec) {

    da->createVector(bdyFlagVec, false, false, 1);//Nodal, Non-ghosted, single dof

    for(int i = 0; i < bdyFlagVec.size(); i++) {
      bdyFlagVec[i] = 0;
    }//initialization loop

    unsigned char *bdyFlagArr = NULL;
    da->vecGetBuffer<unsigned char>(bdyFlagVec,bdyFlagArr,
        false,false,false,1);

    if(da->iAmActive()) {
      //We can only loop over the elements, hence the positive boundary elements
      //will add the flags for the external positive boundary nodes.
      for(da->init<ot::DA_FLAGS::ALL>(); 
          da->curr() < da->end<ot::DA_FLAGS::ALL>();
          da->next<ot::DA_FLAGS::ALL>()) { 
        unsigned char currentFlags;
        bool calledGetNodeIndices = false;
        if(da->isBoundaryOctant(&currentFlags)) {
          //The status of the anchor of any real octant is determined by the
          //negative face boundaries only
          int xNegBdy = (currentFlags & ot::TreeNode::X_NEG_BDY);
          int yNegBdy = (currentFlags & ot::TreeNode::Y_NEG_BDY);
          int zNegBdy = (currentFlags & ot::TreeNode::Z_NEG_BDY);

          unsigned char hnMask = da->getHangingNodeIndex(da->curr());
          if(!(hnMask & 1)) {
            //Anchor is not hanging           	    
            if(xNegBdy) {
              if(yNegBdy && zNegBdy){
                bdyFlagArr[da->curr()] = ot::TreeNode::CORNER_BDY;
              }else if(yNegBdy || zNegBdy) {
                bdyFlagArr[da->curr()] = ot::TreeNode::EDGE_BDY;
              }else {
                bdyFlagArr[da->curr()] = ot::TreeNode::FACE_BDY;
              }
            }else if(yNegBdy) {
              if(zNegBdy) {
                bdyFlagArr[da->curr()] = ot::TreeNode::EDGE_BDY;
              }else {
                bdyFlagArr[da->curr()] = ot::TreeNode::FACE_BDY;
              }
            }else if(zNegBdy) {
              bdyFlagArr[da->curr()] = ot::TreeNode::FACE_BDY;
            }
          }//end if anchor hanging

          if(currentFlags > ot::TreeNode::NEG_POS_DEMARCATION) {
            //Has at least one positive boundary
            //May have negative boundaries as well
            unsigned int indices[8];
            da->getNodeIndices(indices); 
            calledGetNodeIndices = true;
            int xPosBdy = (currentFlags & ot::TreeNode::X_POS_BDY);
            int yPosBdy = (currentFlags & ot::TreeNode::Y_POS_BDY);
            int zPosBdy = (currentFlags & ot::TreeNode::Z_POS_BDY);

            if(!(hnMask & (1 << 1))) {
              bool xBdy = xPosBdy;
              bool yBdy = yNegBdy;
              bool zBdy = zNegBdy;
              if(xBdy) {
                if(yBdy && zBdy){
                  bdyFlagArr[indices[1]] = ot::TreeNode::CORNER_BDY;
                }else if(yBdy || zBdy) {
                  bdyFlagArr[indices[1]] = ot::TreeNode::EDGE_BDY;
                }else {
                  bdyFlagArr[indices[1]] = ot::TreeNode::FACE_BDY;
                }
              }else if(yBdy) {
                if(zBdy) {
                  bdyFlagArr[indices[1]] = ot::TreeNode::EDGE_BDY;
                }else {
                  bdyFlagArr[indices[1]] = ot::TreeNode::FACE_BDY;
                }
              }else if(zBdy) {
                bdyFlagArr[indices[1]] = ot::TreeNode::FACE_BDY;
              }
            }

            if(!(hnMask & (1 << 2))) {
              bool xBdy = xNegBdy;
              bool yBdy = yPosBdy;
              bool zBdy = zNegBdy;
              if(xBdy) {
                if(yBdy && zBdy){
                  bdyFlagArr[indices[2]] = ot::TreeNode::CORNER_BDY;
                }else if(yBdy || zBdy) {
                  bdyFlagArr[indices[2]] = ot::TreeNode::EDGE_BDY;
                }else {
                  bdyFlagArr[indices[2]] = ot::TreeNode::FACE_BDY;
                }
              }else if(yBdy) {
                if(zBdy) {
                  bdyFlagArr[indices[2]] = ot::TreeNode::EDGE_BDY;
                }else {
                  bdyFlagArr[indices[2]] = ot::TreeNode::FACE_BDY;
                }
              }else if(zBdy) {
                bdyFlagArr[indices[2]] = ot::TreeNode::FACE_BDY;
              }
            }

            if(!(hnMask & (1 << 3))) {
              bool xBdy = xPosBdy;
              bool yBdy = yPosBdy;
              bool zBdy = zNegBdy;
              if(xBdy) {
                if(yBdy && zBdy){
                  bdyFlagArr[indices[3]] = ot::TreeNode::CORNER_BDY;
                }else if(yBdy || zBdy) {
                  bdyFlagArr[indices[3]] = ot::TreeNode::EDGE_BDY;
                }else {
                  bdyFlagArr[indices[3]] = ot::TreeNode::FACE_BDY;
                }
              }else if(yBdy) {
                if(zBdy) {
                  bdyFlagArr[indices[3]] = ot::TreeNode::EDGE_BDY;
                }else {
                  bdyFlagArr[indices[3]] = ot::TreeNode::FACE_BDY;
                }
              }else if(zBdy) {
                bdyFlagArr[indices[3]] = ot::TreeNode::FACE_BDY;
              }
            }

            if(!(hnMask & (1 << 4))) {
              bool xBdy = xNegBdy;
              bool yBdy = yNegBdy;
              bool zBdy = zPosBdy;
              if(xBdy) {
                if(yBdy && zBdy){
                  bdyFlagArr[indices[4]] = ot::TreeNode::CORNER_BDY;
                }else if(yBdy || zBdy) {
                  bdyFlagArr[indices[4]] = ot::TreeNode::EDGE_BDY;
                }else {
                  bdyFlagArr[indices[4]] = ot::TreeNode::FACE_BDY;
                }
              }else if(yBdy) {
                if(zBdy) {
                  bdyFlagArr[indices[4]] = ot::TreeNode::EDGE_BDY;
                }else {
                  bdyFlagArr[indices[4]] = ot::TreeNode::FACE_BDY;
                }
              }else if(zBdy) {
                bdyFlagArr[indices[4]] = ot::TreeNode::FACE_BDY;
              }
            }

            if(!(hnMask & (1 << 5))) {
              bool xBdy = xPosBdy;
              bool yBdy = yNegBdy;
              bool zBdy = zPosBdy;
              if(xBdy) {
                if(yBdy && zBdy){
                  bdyFlagArr[indices[5]] = ot::TreeNode::CORNER_BDY;
                }else if(yBdy || zBdy) {
                  bdyFlagArr[indices[5]] = ot::TreeNode::EDGE_BDY;
                }else {
                  bdyFlagArr[indices[5]] = ot::TreeNode::FACE_BDY;
                }
              }else if(yBdy) {
                if(zBdy) {
                  bdyFlagArr[indices[5]] = ot::TreeNode::EDGE_BDY;
                }else {
                  bdyFlagArr[indices[5]] = ot::TreeNode::FACE_BDY;
                }
              }else if(zBdy) {
                bdyFlagArr[indices[5]] = ot::TreeNode::FACE_BDY;
              }
            }

            if(!(hnMask & (1 << 6))) {
              bool xBdy = xNegBdy;
              bool yBdy = yPosBdy;
              bool zBdy = zPosBdy;
              if(xBdy) {
                if(yBdy && zBdy){
                  bdyFlagArr[indices[6]] = ot::TreeNode::CORNER_BDY;
                }else if(yBdy || zBdy) {
                  bdyFlagArr[indices[6]] = ot::TreeNode::EDGE_BDY;
                }else {
                  bdyFlagArr[indices[6]] = ot::TreeNode::FACE_BDY;
                }
              }else if(yBdy) {
                if(zBdy) {
                  bdyFlagArr[indices[6]] = ot::TreeNode::EDGE_BDY;
                }else {
                  bdyFlagArr[indices[6]] = ot::TreeNode::FACE_BDY;
                }
              }else if(zBdy) {
                bdyFlagArr[indices[6]] = ot::TreeNode::FACE_BDY;
              }
            }

            if(!(hnMask & (1 << 7))) {
              bool xBdy = xPosBdy;
              bool yBdy = yPosBdy;
              bool zBdy = zPosBdy;
              if(xBdy) {
                if(yBdy && zBdy){
                  bdyFlagArr[indices[7]] = ot::TreeNode::CORNER_BDY;
                }else if(yBdy || zBdy) {
                  bdyFlagArr[indices[7]] = ot::TreeNode::EDGE_BDY;
                }else {
                  bdyFlagArr[indices[7]] = ot::TreeNode::FACE_BDY;
                }
              }else if(yBdy) {
                if(zBdy) {
                  bdyFlagArr[indices[7]] = ot::TreeNode::EDGE_BDY;
                }else {
                  bdyFlagArr[indices[7]] = ot::TreeNode::FACE_BDY;
                }
              }else if(zBdy) {
                bdyFlagArr[indices[7]] = ot::TreeNode::FACE_BDY;
              }
            }
          }//end if-else has positive boundaries
        }//end if boundary
        if( (!calledGetNodeIndices) && (da->isLUTcompressed()) ) {
          da->updateQuotientCounter();
        }
      }//end for all own elements
    }//end if active

    da->vecRestoreBuffer<unsigned char>(bdyFlagVec,bdyFlagArr,false,false,false,1);
  }//end function

  void includeSiblingsOfBoundary(std::vector<ot::TreeNode>& allBoundaryLeaves, 
      const ot::TreeNode& myFirstOctant, const ot::TreeNode& myLastOctant) {
    PROF_ADD_BDY_SIBLINGS_BEGIN

      std::vector<ot::TreeNode> tmpAllBoundaryLeaves;

    for(unsigned int i = 0; i < allBoundaryLeaves.size(); i++) {
      ot::TreeNode thisOct = allBoundaryLeaves[i];
      ot::TreeNode thisOctParent = thisOct.getParent();
      unsigned int thisCnum = ((unsigned int)(thisOct.getChildNumber()));
      std::vector<ot::TreeNode> siblingsToAdd;
      thisOctParent.addChildren(siblingsToAdd);
      for(unsigned int j=0; j < 8; j++) {
        if( (j != thisCnum) && (j != (7-thisCnum)) ) {
          if( (siblingsToAdd[j] >= myFirstOctant) && 
              (siblingsToAdd[j] <= myLastOctant) ) {
            tmpAllBoundaryLeaves.push_back(siblingsToAdd[j]);
          }
        }
      }   
    }//end for i

    seq::makeVectorUnique<ot::TreeNode>(tmpAllBoundaryLeaves,false);

    std::vector<ot::TreeNode> tmp2Vec;

    unsigned int tmpCnt = 0;
    unsigned int bdyCnt = 0;

    //The two lists are independently sorted and unique, Now we do a linear
    //pass and merge them so that the result is also sorted and unique.

    while ( (tmpCnt < tmpAllBoundaryLeaves.size()) &&
        (bdyCnt < allBoundaryLeaves.size()) ) {
      if ( tmpAllBoundaryLeaves[tmpCnt] < allBoundaryLeaves[bdyCnt] ) {
        tmp2Vec.push_back(tmpAllBoundaryLeaves[tmpCnt]);
        tmpCnt++;
      }else if( tmpAllBoundaryLeaves[tmpCnt] > allBoundaryLeaves[bdyCnt] ) {
        tmp2Vec.push_back(allBoundaryLeaves[bdyCnt]);
        bdyCnt++;
      }else {
        tmp2Vec.push_back(allBoundaryLeaves[bdyCnt]);
        bdyCnt++;
        tmpCnt++; //tmpCnt must also be incremented to preserve uniqueness.
      }
    }

    while (bdyCnt < allBoundaryLeaves.size()) {
      tmp2Vec.push_back(allBoundaryLeaves[bdyCnt]);
      bdyCnt++;
    }

    while (tmpCnt < tmpAllBoundaryLeaves.size()) {
      tmp2Vec.push_back(tmpAllBoundaryLeaves[tmpCnt]);
      tmpCnt++;
    }

    allBoundaryLeaves = tmp2Vec;

    tmp2Vec.clear();
    tmpAllBoundaryLeaves.clear();

    PROF_ADD_BDY_SIBLINGS_END
  }//end function

  void prepareAprioriCommMessagesInDAtype1(const std::vector<ot::TreeNode>& in,
      std::vector<ot::TreeNode>& allBoundaryLeaves, std::vector<ot::TreeNode>& blocks,
      const std::vector<ot::TreeNode>& allBlocks, int myRank, int npes, int* sendCnt,
      std::vector<std::vector<unsigned int> >& sendNodes) {
    PROF_DA_APRIORI_COMM_BEGIN

      std::vector<unsigned int> bdy2elem;

    unsigned int bdyCnt = 0;
    unsigned int allCnt = 0;

    while (bdyCnt < allBoundaryLeaves.size()) {
      if ( allBoundaryLeaves[bdyCnt] < in[allCnt]) {
        bdyCnt++;
      } else if ( allBoundaryLeaves[bdyCnt] > in[allCnt]) {
        allCnt++;
      } else {
        //Both are equal.		
        bdy2elem.push_back(allCnt);
        bdyCnt++;
      }
    }

    //This step is necessary because some elements of allBoundaryLeaves were not
    //copied from the "in" vector, instead we generated them directly. So these
    //octants will not have correct flags set. So we find the corresponding copy
    //in the "in" vector.  
    allBoundaryLeaves.clear();
    for(unsigned int i = 0; i < bdy2elem.size(); i++) {
      allBoundaryLeaves.push_back(in[bdy2elem[i]]);
    }

    // 3. Reduce the list of global blocks to a smaller list that only has the
    // blocks which neighbour ones own blocks.

    //First mark your own blocks as being singular or not.
    unsigned int allBdyCnt = 0;
    for (unsigned int i = 0; i < blocks.size(); i++) {
      while( (allBdyCnt < allBoundaryLeaves.size()) &&
          (allBoundaryLeaves[allBdyCnt] < blocks[i]) ) {
        allBdyCnt++;
      }
      //Note, even if a block is Singular but is
      //not a ghost candidates (i.e., it is completely internal).
      // It will be treated as being NOT singular.
      bool isSingular = false;
      if( (allBdyCnt < allBoundaryLeaves.size()) &&
          (allBoundaryLeaves[allBdyCnt] == blocks[i]) ) {
        isSingular = true;
      }
      if(isSingular) {
        blocks[i].setWeight(1);
      }else {
        blocks[i].setWeight(0);
      }
    }

    std::vector<ot::TreeNode> myNhBlocks;
    for (int j = 0; j < blocks.size(); j++) {
      unsigned int myMaxX;
      unsigned int myMaxY;
      unsigned int myMaxZ;
      unsigned int myMinX;
      unsigned int myMinY;
      unsigned int myMinZ;
      if(blocks[j].getWeight()) {
        //Since the subsequent selection is done based on the octant's parent's
        //neighbours, Singular blocks must be handled differently.
        myMaxX = blocks[j].getParent().maxX();
        myMaxY = blocks[j].getParent().maxY();
        myMaxZ = blocks[j].getParent().maxZ();
        myMinX = blocks[j].getParent().minX();
        myMinY = blocks[j].getParent().minY();
        myMinZ = blocks[j].getParent().minZ();
      }else {
        myMaxX = blocks[j].maxX();
        myMaxY = blocks[j].maxY();
        myMaxZ = blocks[j].maxZ();
        myMinX = blocks[j].minX();
        myMinY = blocks[j].minY();
        myMinZ = blocks[j].minZ();
      }
      double myLenX = (double)(myMaxX-myMinX);
      double myLenY = (double)(myMaxY-myMinY);
      double myLenZ = (double)(myMaxZ-myMinZ);
      double myXc  = ((double)(myMinX + myMaxX))/2.0;
      double myYc  = ((double)(myMinY + myMaxY))/2.0;
      double myZc  = ((double)(myMinZ + myMaxZ))/2.0;
      for (int k = 0; k < allBlocks.size(); k++) {
        if ( (allBlocks[k] >= blocks[0]) && 
            (allBlocks[k] <= blocks[blocks.size()-1]) ) {
          //ignore my own blocks
          continue;
        }
        unsigned int hisMinX = allBlocks[k].minX();
        unsigned int hisMaxX = allBlocks[k].maxX();
        unsigned int hisMinY = allBlocks[k].minY();
        unsigned int hisMaxY = allBlocks[k].maxY();
        unsigned int hisMinZ = allBlocks[k].minZ();
        unsigned int hisMaxZ = allBlocks[k].maxZ();
        double hisLenX = (double)(hisMaxX-hisMinX);
        double hisLenY = (double)(hisMaxY-hisMinY);
        double hisLenZ = (double)(hisMaxZ-hisMinZ);
        double hisXc  = ((double)(hisMinX + hisMaxX))/2.0;
        double hisYc  = ((double)(hisMinY + hisMaxY))/2.0;
        double hisZc  = ((double)(hisMinZ + hisMaxZ))/2.0;
        double deltaX = ( (hisXc > myXc) ? (hisXc - myXc) : (myXc - hisXc ) );
        double deltaY = ( (hisYc > myYc) ? (hisYc - myYc) : (myYc - hisYc ) );
        double deltaZ = ( (hisZc > myZc) ? (hisZc - myZc) : (myZc - hisZc ) );

        //Note: This test will pass if the octants intersect (say one is an
        //ancestor of another) or they simply touch externally.

        if ((deltaX <= ((hisLenX+myLenX)/2.0)) &&
            (deltaY <= ((hisLenY+myLenY)/2.0)) &&
            (deltaZ <= ((hisLenZ+myLenZ)/2.0))) {
          //We touch
          myNhBlocks.push_back(allBlocks[k]);
        }//end if
      }//end for k
    }//end for j

    //This also sorts myNhBlocks
    seq::makeVectorUnique<ot::TreeNode>(myNhBlocks,false);

    sendNodes.resize(npes);
    for (int i=0; i < npes; i++) {
      sendNodes[i].clear();
      sendCnt[i] = 0;
    } 

    for (int j=0; j<allBoundaryLeaves.size(); j++) {
      //It is important to make the selection using the parent. This tackles
      //the cases where some nodes of an octant are hanging and so even if the octant itself
      //does not touch a processor, its parent does and so the replacement for
      //the hanging node will be mapped to that processor. This also handles
      //the case where two or more siblings belong to different processors. By
      //making the test using the parent they will be sent to each other.
      //Moreover, some of the nodes of these children could be hanging and the
      //replacement could belong to a third processor. Even in this case, the
      //octants will be sent to the correct processors. 

      ot::TreeNode parNode = allBoundaryLeaves[j].getParent();
      unsigned int myMaxX = parNode.maxX();
      unsigned int myMaxY = parNode.maxY();
      unsigned int myMaxZ = parNode.maxZ();
      unsigned int myMinX = parNode.minX();
      unsigned int myMinY = parNode.minY();
      unsigned int myMinZ = parNode.minZ();
      double myLenX = (double)(myMaxX-myMinX);
      double myLenY = (double)(myMaxY-myMinY);
      double myLenZ = (double)(myMaxZ-myMinZ);
      double myXc  = ((double)(myMinX + myMaxX))/2.0;
      double myYc  = ((double)(myMinY + myMaxY))/2.0;
      double myZc  = ((double)(myMinZ + myMaxZ))/2.0;
      unsigned int lastP = npes;
      for (int k = 0; k < myNhBlocks.size(); k++) {
        if (myNhBlocks[k].getWeight() == lastP) {
          continue;
        }
        unsigned int hisMinX = myNhBlocks[k].minX();
        unsigned int hisMaxX = myNhBlocks[k].maxX();
        unsigned int hisMinY = myNhBlocks[k].minY();
        unsigned int hisMaxY = myNhBlocks[k].maxY();
        unsigned int hisMinZ = myNhBlocks[k].minZ();
        unsigned int hisMaxZ = myNhBlocks[k].maxZ();
        double hisLenX = (double)(hisMaxX-hisMinX);
        double hisLenY = (double)(hisMaxY-hisMinY);
        double hisLenZ = (double)(hisMaxZ-hisMinZ);
        double hisXc  = ((double)(hisMinX + hisMaxX))/2.0;
        double hisYc  = ((double)(hisMinY + hisMaxY))/2.0;
        double hisZc  = ((double)(hisMinZ + hisMaxZ))/2.0;
        double deltaX = ( (hisXc > myXc) ? (hisXc - myXc) : (myXc - hisXc ) );
        double deltaY = ( (hisYc > myYc) ? (hisYc - myYc) : (myYc - hisYc ) );
        double deltaZ = ( (hisZc > myZc) ? (hisZc - myZc) : (myZc - hisZc ) );

        //Note: This test will pass if the octants intersect (say one is an
        //ancestor of another) or they simply touch externally.

        if ((deltaX <= ((hisLenX+myLenX)/2.0)) &&
            (deltaY <= ((hisLenY+myLenY)/2.0)) &&
            (deltaZ <= ((hisLenZ+myLenZ)/2.0))) {
          //We touch
          sendNodes[myNhBlocks[k].getWeight()].push_back(bdy2elem[j]); 
          sendCnt[myNhBlocks[k].getWeight()]++;
          lastP = myNhBlocks[k].getWeight();
        }//end if
      }//end for k
    }//end for j

    PROF_DA_APRIORI_COMM_END
  }//end function

  void prepareAprioriCommMessagesInDAtype2(const std::vector<ot::TreeNode>& in,
      std::vector<ot::TreeNode>& allBoundaryLeaves, std::vector<ot::TreeNode>& blocks,
      const std::vector<ot::TreeNode>& minsOfBlocks, int myRank, int npes, int* sendCnt,
      std::vector<std::vector<unsigned int> >& sendNodes) {
    PROF_DA_APRIORI_COMM_BEGIN

      std::vector<unsigned int> bdy2elem;

    unsigned int bdyCnt = 0;
    unsigned int allCnt = 0;
    unsigned int maxDepth = in[0].getMaxDepth();
    unsigned int dim = in[0].getDim(); 

    while (bdyCnt < allBoundaryLeaves.size()) {
      if ( allBoundaryLeaves[bdyCnt] < in[allCnt]) {
        bdyCnt++;
      } else if ( allBoundaryLeaves[bdyCnt] > in[allCnt]) {
        allCnt++;
      } else {
        //Both are equal.		
        bdy2elem.push_back(allCnt);
        bdyCnt++;
      }
    }

    //This step is necessary because some elements of allBoundaryLeaves were not
    //copied from the "in" vector, instead we generated them directly. So these
    //octants will not have correct flags set. So we find the corresponding copy
    //in the "in" vector.  
    allBoundaryLeaves.clear();
    for(unsigned int i = 0; i < bdy2elem.size(); i++) {
      allBoundaryLeaves.push_back(in[bdy2elem[i]]);
    }

    sendNodes.resize(npes);
    for (int i = 0; i < npes; i++) {
      sendNodes[i].clear();
      sendCnt[i] = 0;
    } 

    for (int j = 0; j < allBoundaryLeaves.size(); j++) {
      //It is important to make the selection using the parent. This tackles
      //the cases where some nodes of an octant are hanging and so even if the octant itself
      //does not touch a processor, its parent does and so the replacement for
      //the hanging node will be mapped to that processor. This also handles
      //the case where two or more siblings belong to different processors. By
      //making the test using the parent they will be sent to each other.
      //Moreover, some of the nodes of these children could be hanging and the
      //replacement could belong to a third processor. Even in this case, the
      //octants will be sent to the correct processors. 

      //1. We must not miss any pre-ghost elements that point to one of our own
      //octants, i.e. any pre-ghost element that we integrate over. This is
      //because we need to build the nlist for these octants as well.
      //2. We might get some extra ghosts they will be marked as FOREIGN and
      //will be ignored.
      //3. We must try to get as many post-ghost octants as possible to avoid
      //misses later and hence reduce subsequent communication. We should get
      //all direct post-ghost octants.
      //4. Post-ghosts are read-only and so it is sufficient to test for
      //post-ghost using its anchor
      std::vector<ot::TreeNode> myVertices;
      std::vector<ot::TreeNode> parVertices;
      std::vector<ot::TreeNode> anchorMirrors;

      unsigned int myX = allBoundaryLeaves[j].getX();
      unsigned int myY = allBoundaryLeaves[j].getY();
      unsigned int myZ = allBoundaryLeaves[j].getZ();

      //keys to check if you are a pre-ghost
      //Positive boundaries will not be pre-ghosts
      if(!(allBoundaryLeaves[j].getFlag() & ot::TreeNode::BOUNDARY)) {
        unsigned int myLev = allBoundaryLeaves[j].getLevel();
        unsigned int mySz = (1u<<(maxDepth - myLev));

        //All vertices except my anchor. Since my anchor belongs to my processor
        myVertices.push_back(ot::TreeNode((myX + mySz), myY, myZ, maxDepth, dim, maxDepth));
        myVertices.push_back(ot::TreeNode(myX, (myY + mySz), myZ, maxDepth, dim, maxDepth));
        myVertices.push_back(ot::TreeNode((myX + mySz), (myY + mySz), myZ, maxDepth, dim, maxDepth));
        myVertices.push_back(ot::TreeNode(myX, myY, (myZ + mySz), maxDepth, dim, maxDepth));
        myVertices.push_back(ot::TreeNode((myX + mySz), myY, (myZ + mySz), maxDepth, dim, maxDepth));
        myVertices.push_back(ot::TreeNode(myX, (myY + mySz), (myZ + mySz), maxDepth, dim, maxDepth));
        myVertices.push_back(ot::TreeNode((myX + mySz), (myY + mySz), (myZ + mySz), maxDepth, dim, maxDepth));

        ot::TreeNode parNode = allBoundaryLeaves[j].getParent();
        unsigned int myCnum = allBoundaryLeaves[j].getChildNumber();
        unsigned int parX = parNode.getX();
        unsigned int parY = parNode.getY();
        unsigned int parZ = parNode.getZ();
        unsigned int parLev = parNode.getLevel();
        unsigned int parSz = (1u<<(maxDepth - parLev));

        //vertices numbers myCnum and (7 - myCnum) can't be hanging and so they
        //will not be mapped to the corresponding vertices of my parent
        if( (myCnum != 0) && (myCnum != 7) ) {
          parVertices.push_back(ot::TreeNode(parX, parY, parZ, maxDepth, dim, maxDepth));
        }
        if( (myCnum != 1) && (myCnum != 6) ) {
          parVertices.push_back(ot::TreeNode((parX + parSz), parY, parZ, maxDepth, dim, maxDepth));
        }
        if( (myCnum != 2) && (myCnum != 5) ) {
          parVertices.push_back(ot::TreeNode(parX, (parY + parSz), parZ, maxDepth, dim, maxDepth));
        }
        if( (myCnum != 3) && (myCnum != 4) ) {
          parVertices.push_back(ot::TreeNode((parX + parSz), (parY + parSz), parZ, maxDepth, dim, maxDepth));
        }
        if( (myCnum != 4) && (myCnum != 3) ) {
          parVertices.push_back(ot::TreeNode(parX, parY, (parZ + parSz), maxDepth, dim, maxDepth));
        }
        if( (myCnum != 5) && (myCnum != 2) ) {
          parVertices.push_back(ot::TreeNode((parX + parSz), parY, (parZ + parSz), maxDepth, dim, maxDepth));
        }
        if( (myCnum != 6) && (myCnum != 1) ) {
          parVertices.push_back(ot::TreeNode(parX, (parY + parSz), (parZ + parSz), maxDepth, dim, maxDepth));
        }
        if( (myCnum != 7) && (myCnum != 0) ) {
          parVertices.push_back(ot::TreeNode((parX + parSz), (parY + parSz),
                (parZ + parSz), maxDepth, dim, maxDepth));
        }

      }//end if positive boundary

      //Keys to check if you are a post-ghost
      //If the anchor is hanging we need not send it. When some processor searches
      //for this node as its primary key and does not find it, it will be
      //understood that this is hanging  
      if(allBoundaryLeaves[j].getFlag() & ot::TreeNode::NODE) {
        //-x
        if( myX ) {
          anchorMirrors.push_back(ot::TreeNode((myX - 1), myY, myZ, 
                maxDepth, dim, maxDepth));
        }
        //-y
        if( myY ) {
          anchorMirrors.push_back(ot::TreeNode(myX, (myY - 1), myZ, 
                maxDepth, dim, maxDepth));
        }
        //-z
        if( myZ ) {
          anchorMirrors.push_back(ot::TreeNode(myX, myY, (myZ - 1), 
                maxDepth, dim, maxDepth));
        }
        //-xy
        if( myX && myY ) {
          anchorMirrors.push_back(ot::TreeNode((myX - 1), (myY - 1), myZ, 
                maxDepth, dim, maxDepth));
        }
        //-yz
        if( myY && myZ ) {
          anchorMirrors.push_back(ot::TreeNode(myX, (myY - 1), (myZ - 1),
                maxDepth, dim, maxDepth));
        }
        //-zx
        if( myZ && myX ) {
          anchorMirrors.push_back(ot::TreeNode((myX - 1), myY, (myZ - 1),
                maxDepth, dim, maxDepth));
        }
        //-xyz
        if( myX && myY && myZ ) {
          anchorMirrors.push_back(ot::TreeNode((myX - 1), (myY - 1), (myZ - 1),
                maxDepth, dim, maxDepth));
        }
      }//end if hanging anchor

      std::vector<unsigned int> pIds;
      for(int k = 0; k < parVertices.size(); k++) {
        unsigned int idx;
        seq::maxLowerBound<ot::TreeNode>(minsOfBlocks, parVertices[k], idx, NULL, NULL);
        pIds.push_back(idx);
      }

      for(int k = 0; k < anchorMirrors.size(); k++) {
        unsigned int idx;
        seq::maxLowerBound<ot::TreeNode>(minsOfBlocks, anchorMirrors[k], idx, NULL, NULL);
        pIds.push_back(idx);
      }

      for(int k = 0; k < myVertices.size(); k++) {
        unsigned int idx;
        seq::maxLowerBound<ot::TreeNode>(minsOfBlocks, myVertices[k], idx, NULL, NULL);
        pIds.push_back(idx);
      }//end for k

      //Do not send the same octant to the same processor twice 
      seq::makeVectorUnique<unsigned int>(pIds, false);

      for(int k = 0; k < pIds.size(); k++) {
        //Send to processor pIds[k] 
        //Only send to other processors  
        if(pIds[k] != myRank) {
          sendNodes[pIds[k]].push_back(bdy2elem[j]); 
          sendCnt[pIds[k]]++;
        }
      }//end for k
    }//end for j

    PROF_DA_APRIORI_COMM_END
  }//end function

}//end namespace


