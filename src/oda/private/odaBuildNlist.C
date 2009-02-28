
/**
  @file odaBuildNlist.C
  @brief Octree Meshing: Building Element-to-node-mappings
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  @author Hari Sundar, hsundar@gmail.com
  */

#include "oda.h"
#include "parUtils.h"
#include "seqUtils.h"
#include "colors.h"
#include "nodeAndRanks.h"
#include "testUtils.h"

#ifdef __DEBUG__
#ifndef __DEBUG_DA_NLIST__
#define __DEBUG_DA_NLIST__
#endif
#endif

#ifdef __DEBUG_DA__
#ifndef __DEBUG_DA_NLIST__
#define __DEBUG_DA_NLIST__
#endif
#endif

#ifdef __DEBUG_DA_NLIST__
#ifndef __MEASURE_BUILD_NLIST__
#define __MEASURE_BUILD_NLIST__
#endif
#endif

//Warning: DO NOT MEDDLE WITH THIS PIECE OF CODE!!!

namespace ot {

#define CHECK_FAST_MAX_LOWER_BOUND(arr,key,fastIdx,fastResult) { \
  unsigned int tmpMlbIdx;\
  bool tmpMlbResult;\
  tmpMlbResult = seq::maxLowerBound<ot::TreeNode> (arr, key, tmpMlbIdx,NULL,NULL);\
  assert(tmpMlbResult == fastResult);\
  assert(tmpMlbIdx == fastIdx);\
}

//Build Node list using 4-way searches...
void DA::buildNodeList(std::vector<ot::TreeNode> &in) {
#ifdef __PROF_WITH_BARRIER__
  MPI_Barrier(m_mpiCommActive);
#endif
  PROF_BUILD_NLIST_BEGIN
    // everybody except for the boundary and positive ghosts should be elements,
    // This means that anything that is not a boundary should be an element.
    // The only extra elements added are the ghost elements.
    // Initially will store all 8 indices. will not be compressed.

    // compute number of elements.
    unsigned int nelem = m_uiPreGhostElementSize + m_uiElementSize;

  std::vector<unsigned int> nlist;

  /*
     The first iteration is for pre-ghosts only.
     The second iteration is for own elements. 
     The first iteration does not use 4-way searches. The second iteration does.
     */
  for(unsigned int numFullLoopCtr = 0; numFullLoopCtr < 2; numFullLoopCtr++) {

    // Some storage ...
    std::vector<ot::TreeNode> primaryKeys;
    std::vector<ot::TreeNode> secondaryKeys;
    std::vector<ot::TreeNode> extraAtEnd;

#ifdef __DEBUG_DA_NLIST__
    MPI_Barrier(m_mpiCommActive);
    std::vector<ot:: TreeNode > checkSecondRing;
    std::vector<ot::TreeNode> chkMissedPrimary;
#endif

    unsigned int iLoopSt,iLoopEnd;
    if( numFullLoopCtr == 0) {
      //PreGhosts only...
      iLoopSt = 0;
      iLoopEnd = m_uiPreGhostElementSize;
    }else {
      //nelem and m_uiPreGhostElementSize would have been changed in the
      //first iteration. It's ok. Use the new values only.
      //Own elements only...
      iLoopSt = m_uiPreGhostElementSize;
      iLoopEnd = nelem;
    }

#ifdef __DEBUG_DA_NLIST__
    MPI_Barrier(m_mpiCommActive);
    if(!m_iRankActive) {
      std::cout<<"numFullLoopCtr: "<<numFullLoopCtr<<std::endl;
    }
    std::cout<<m_iRankActive<<": preElemSz: "<<m_uiPreGhostElementSize
      <<" elemBeg: "<<m_uiElementBegin<<" elemEnd: "<<m_uiElementEnd
      <<" postGhostBegin: "<<m_uiPostGhostBegin
      <<" locBufferSz: "<<m_uiLocalBufferSize
      <<" nelem: "<<nelem
      <<" iLoopSt: "<<iLoopSt
      <<" iLoopEnd: "<<iLoopEnd<<std::endl;
    assert(m_uiElementBegin < m_uiPostGhostBegin);
    std::cout<<m_iRankActive<<" my First Octant(elem/Bnd): "<<in[m_uiElementBegin]<<std::endl;
    MPI_Barrier(m_mpiCommActive);
#endif


    //Malloc...
    //In the first iteration, we need to only allocate for pre-ghosts.
    //The second iteration, it would already be resized to the correct size.
    //So there would be no change.
    if( numFullLoopCtr == 0) {
      nlist.resize(8*iLoopEnd);
      m_ucpLutMasks.resize(2*iLoopEnd);
    }

    //Loop through all the elements in this set and set LUTs.
    for (unsigned int i = iLoopSt; i < iLoopEnd; i++) {

      m_ucpLutMasks[2*i + 1] = 0;

      std::vector<ot::TreeNode> nodeLocations(8);
      std::vector<ot::TreeNode> parNodeLocations(8);
      // get basic info ...
      unsigned int d   = in[i].getLevel();
      unsigned int x = in[i].getX();
      unsigned int y = in[i].getY();
      unsigned int z = in[i].getZ();

      //Cryptic Implementation:
      unsigned int parX = ( ( x >> ( m_uiMaxDepth - d + 1 ) ) << ( m_uiMaxDepth - d + 1 ) );
      unsigned int parY = ( ( y >> ( m_uiMaxDepth - d + 1 ) ) << ( m_uiMaxDepth - d + 1 ) );
      unsigned int parZ = ( ( z >> ( m_uiMaxDepth - d + 1 ) ) << ( m_uiMaxDepth - d + 1 ) );

      unsigned int sz = 1u << (m_uiMaxDepth - d);
      unsigned int len_par = (unsigned int)(1u<<( m_uiMaxDepth  - d +1 ) );

      unsigned int a = x % len_par;
      unsigned int b = y % len_par;
      unsigned int c = z % len_par;

      a /= sz;
      b /= sz;
      c /= sz;

      unsigned int ch_num = (4*c + 2*b + a);

#ifdef __DEBUG_DA_NLIST__
      if ( !ch_num || (ch_num==7) ) {
        if ( !(in[i].getFlag() & ot::TreeNode::NODE) ) {
          std::cerr << RED"Nodes are marked wrongly "NRM << std::endl;
          assert(false);
        }
      }
#endif

      bool found[8];
      // haven't found anything yet. Set Default values.
      for (unsigned int k = 0; k < 8; k++) {
        nlist[8*i + k] = m_uiLocalBufferSize;
        found[k] = false;
      }//end for k

      //~~~~~~~~~~~~~~~~~~~~~~~NEGATIVE SEARCH~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(numFullLoopCtr == 1) {
        unsigned int idnx, idny, idnz;
        //Basic Idea of a negative search. 
        //1. Search and find elements on the negative faces. 
        //2. If the result is a brother and the Vtx in question is hanging, swap and copy.
        //3. In all other cases, simply copy.
        //4. Note, copy only if the nlist is pointing to a valid location

        bool foundNegX=false, foundNegY=false, foundNegZ=false;

        // first lets do the 3 negative searches and we'll decide how to use it later ...
        // search negative X
        if (x) {
          ot::TreeNode knx( x-1, y, z, m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
          foundNegX = seq::maxLowerBound<ot::TreeNode> (in, knx, idnx,NULL,&i);
#ifdef __DEBUG_DA_NLIST__
          CHECK_FAST_MAX_LOWER_BOUND(in,knx,idnx,foundNegX)
#endif

            if ( foundNegX && (!( (in[idnx].isAncestor(knx) ) || (in[idnx] == knx) )) ) {
              foundNegX=false;
            }
          if ( foundNegX && (m_ucpLutMasks[2*idnx+1] == ot::DA_FLAGS::FOREIGN) ) {
            foundNegX=false;
          }
        }

        // search negative Y
        if (y) {
          ot::TreeNode kny( x, y-1, z, m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
          foundNegY = seq::maxLowerBound<ot::TreeNode> (in, kny, idny,NULL,&i) ;
#ifdef __DEBUG_DA_NLIST__
          CHECK_FAST_MAX_LOWER_BOUND(in,kny,idny,foundNegY)
#endif
            if ( foundNegY && (!( (in[idny].isAncestor(kny) ) || (in[idny] == kny) )) ) {
              foundNegY=false;
            }
          if ( foundNegY && (m_ucpLutMasks[2*idny+1] == ot::DA_FLAGS::FOREIGN) ) {
            foundNegY=false;
          }
        }

        // search negative Z
        if (z) {
          ot::TreeNode knz( x, y, z-1, m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
          foundNegZ = seq::maxLowerBound<ot::TreeNode> (in, knz, idnz,NULL,&i) ;
#ifdef __DEBUG_DA_NLIST__
          CHECK_FAST_MAX_LOWER_BOUND(in,knz,idnz,foundNegZ)
#endif
            if ( foundNegZ && (!( (in[idnz].isAncestor(knz) ) || (in[idnz] == knz) )) ) {
              foundNegZ=false;
            }
          if ( foundNegZ && (m_ucpLutMasks[2*idnz+1] == ot::DA_FLAGS::FOREIGN) ) {
            foundNegZ=false;
          }
        }

        // How we use the results of the negative search depends on the child number ...
        // if here my is the current element and x, y, and z are the elements in the 
        // respective negative directions, then the default mapping sans corrections is
        // 
        // my[0] = x[1] = y[2] = z[4]
        // my[1] = y[3] = z[5]
        // my[2] = x[3] = z[6]
        // my[3] = z[7];
        // my[4] = x[5] = y[6]
        // my[5] = y[7];
        // my[6] = x[7];
        // 
        // The corrections will appear for all children except for child zero ... 
        // and will be explained when they are performed.

#define NEG_SEARCH_BLOCK1(idx,n1l,n1r,n2l,n2r,n3l,n3r,n4l,n4r) {\
  if(nlist[8*idx+n1r] < m_uiLocalBufferSize) {\
    nlist[8*i+n1l] = nlist[8*idx+n1r];\
    found[n1l] = true;\
  }\
  if(nlist[8*idx+n2r] < m_uiLocalBufferSize) {\
    nlist[8*i+n2l] = nlist[8*idx+n2r];\
    found[n2l] = true;\
  }\
  if(nlist[8*idx+n3r] < m_uiLocalBufferSize) {\
    nlist[8*i+n3l] = nlist[8*idx+n3r];\
    found[n3l] = true;\
  }\
  if(nlist[8*idx+n4r] < m_uiLocalBufferSize) {\
    nlist[8*i+n4l] = nlist[8*idx+n4r];\
    found[n4l] = true;\
  }\
}

#define NEG_SEARCH_BLOCK2(idx,n1l,n1r,n2l,n2r,n3l,n3r,n4l,n4r) {\
  unsigned char negMask;\
  if ( in[i].getLevel() != in[idx].getLevel()) {\
    negMask = 0;\
  } else {\
    negMask = m_ucpLutMasks[2*idx+1];\
  }\
  if (negMask != ot::DA_FLAGS::FOREIGN) {\
    if ( negMask & (1 << n1r)) {\
      if(nlist[8*idx+n1l] < m_uiLocalBufferSize) {\
        nlist[8*i+n1l] = nlist[8*idx+n1l];\
        found[n1l] = true;\
      }\
    } else {\
      if(nlist[8*idx+n1r] < m_uiLocalBufferSize) {\
        nlist[8*i+n1l] = nlist[8*idx+n1r];\
        found[n1l] = true;\
      }\
    }\
    if ( negMask & (1 << n2r)) {\
      if(nlist[8*idx+n2l] < m_uiLocalBufferSize) {\
        nlist[8*i+n2l] = nlist[8*idx+n2l];\
        found[n2l] = true;\
      }\
    } else {\
      if(nlist[8*idx+n2r] < m_uiLocalBufferSize) {\
        nlist[8*i+n2l] = nlist[8*idx+n2r];\
        found[n2l] = true;\
      }\
    }\
    if ( negMask & (1 << n3r)) {\
      if(nlist[8*idx+n3l] < m_uiLocalBufferSize) {\
        nlist[8*i+n3l] = nlist[8*idx+n3l];\
        found[n3l] = true;\
      }\
    } else {\
      if(nlist[8*idx+n3r] < m_uiLocalBufferSize) {\
        nlist[8*i+n3l] = nlist[8*idx+n3r];\
        found[n3l] = true;\
      }\
    }\
    if ( negMask & (1 << n4r)) {\
      if(nlist[8*idx+n4l] < m_uiLocalBufferSize) {\
        nlist[8*i+n4l] = nlist[8*idx+n4l];\
        found[n4l] = true;\
      }\
    } else {\
      if(nlist[8*idx+n4r] < m_uiLocalBufferSize) {\
        nlist[8*i+n4l] = nlist[8*idx+n4r];\
        found[n4l] = true;\
      }\
    }\
  }\
}

#define NEG_SEARCH_BLOCK1X NEG_SEARCH_BLOCK1(idnx,0,1,2,3,4,5,6,7)
#define NEG_SEARCH_BLOCK2X NEG_SEARCH_BLOCK2(idnx,0,1,2,3,4,5,6,7)

#define NEG_SEARCH_BLOCK1Y NEG_SEARCH_BLOCK1(idny,0,2,1,3,4,6,5,7)
#define NEG_SEARCH_BLOCK2Y NEG_SEARCH_BLOCK2(idny,0,2,1,3,4,6,5,7)

#define NEG_SEARCH_BLOCK1Z NEG_SEARCH_BLOCK1(idnz,0,4,1,5,2,6,3,7)
#define NEG_SEARCH_BLOCK2Z NEG_SEARCH_BLOCK2(idnz,0,4,1,5,2,6,3,7)

switch (ch_num) {
  case 0: {
            //No Negative brothers.
            if (foundNegX) {
              NEG_SEARCH_BLOCK1X
            }
            if (foundNegY) {
              NEG_SEARCH_BLOCK1Y
            }
            if (foundNegZ) {
              NEG_SEARCH_BLOCK1Z
            }
            break;
          }
  case 1: {
            //negX could be a brother. 
            if (foundNegX) {
              NEG_SEARCH_BLOCK2X
            }
            if (foundNegY) {
              NEG_SEARCH_BLOCK1Y
            }
            if (foundNegZ) {
              NEG_SEARCH_BLOCK1Z
            }
            break;
          }
  case 2: {
            //negY could be a brother.
            if (foundNegX) {
              NEG_SEARCH_BLOCK1X
            }
            if (foundNegY) {
              NEG_SEARCH_BLOCK2Y
            }
            if (foundNegZ) {
              NEG_SEARCH_BLOCK1Z
            }
            break;
          }
  case 3: {
            //negX and negY could be brothers.
            if (foundNegX) {
              NEG_SEARCH_BLOCK2X
            }
            if (foundNegY) {
              NEG_SEARCH_BLOCK2Y
            }
            if (foundNegZ) {
              NEG_SEARCH_BLOCK1Z
            }
            break;
          }
  case 4: {
            //negZ could be a brother.
            if (foundNegX) {
              NEG_SEARCH_BLOCK1X
            }
            if (foundNegY) { 
              NEG_SEARCH_BLOCK1Y
            }
            if (foundNegZ) {
              NEG_SEARCH_BLOCK2Z
            }
            break;
          }
  case 5: {
            //negX and negZ could be brothers.
            if (foundNegX) {
              NEG_SEARCH_BLOCK2X
            }
            if (foundNegY) {  
              NEG_SEARCH_BLOCK1Y
            }
            if (foundNegZ) {
              NEG_SEARCH_BLOCK2Z
            }
            break;
          }
  case 6: {
            //negY and negZ could be brothers.
            if (foundNegX) {
              NEG_SEARCH_BLOCK1X
            }
            if (foundNegY) {
              NEG_SEARCH_BLOCK2Y
            }
            if (foundNegZ) {
              NEG_SEARCH_BLOCK2Z
            }
            break;
          }
  case 7: {
            //negX negY and negZ could be brothers.
            if (foundNegX) {
              NEG_SEARCH_BLOCK2X
            }
            if (foundNegY) {
              NEG_SEARCH_BLOCK2Y
            }
            if (foundNegZ) {
              NEG_SEARCH_BLOCK2Z
            }
            break;
          }
  default: {
             std::cerr << "Wrong Child Number " << ch_num << std::endl;
             assert(false);
             break;
           }
}//end cases
}//end if for negative

#undef NEG_SEARCH_BLOCK1X
#undef NEG_SEARCH_BLOCK2X
#undef NEG_SEARCH_BLOCK1Y
#undef NEG_SEARCH_BLOCK2Y 
#undef NEG_SEARCH_BLOCK1Z 
#undef NEG_SEARCH_BLOCK2Z
#undef NEG_SEARCH_BLOCK1
#undef NEG_SEARCH_BLOCK2

//~~~~~~~~~~~~~~~~~~~~~POSITIVE SEARCH~~~~~~~~~~~~~~~~~~~~~~
#ifdef __DEBUG_DA_NLIST__
#define DEBUG_CHECK_FAST_MAX_LOWER_BOUND(arr,key,fastIdx,fastResult) CHECK_FAST_MAX_LOWER_BOUND(arr,key,fastIdx,fastResult)

#define POS_SEARCH_DEBUG_BLOCK1(debugV1,debugV2) {\
  if ( (ch_num == debugV1) || (ch_num == debugV2) ) {\
    std::cout << "Failing for index " << i << " with child num " << ch_num << std::endl;\
    assert(false);\
  }\
}

#define  POS_SEARCH_DEBUG_BLOCK2 {\
  assert(sKey > in[m_uiPostGhostBegin - 1]);\
  chkMissedPrimary.push_back(sKey);\
}

#define POS_SEARCH_DEBUG_BLOCK3 {\
  assert( (idx+k) < m_uiLocalBufferSize );\
  assert( in[idx+k].getParent() == newKey );\
  assert( !(in[idx+k].getFlag() & ot::TreeNode::NODE) );\
  assert( ((idx+k) < m_uiElementBegin) || (((idx+k) >= m_uiPostGhostBegin)) );\
}

#else
#define DEBUG_CHECK_FAST_MAX_LOWER_BOUND(arr,key,fastIdx,fastResult) 
#define POS_SEARCH_DEBUG_BLOCK1(debugV1,debugV2)
#define POS_SEARCH_DEBUG_BLOCK2
#define POS_SEARCH_DEBUG_BLOCK3
#endif

#define POS_SECONDARY_SEARCH_BLOCK(debugV1,debugV2) {\
  /* if this is not a node, i.e., it is hanging*/\
  POS_SEARCH_DEBUG_BLOCK1(debugV1,debugV2)\
  /* All other cNums are anchored at the parent+(2*sz,0,0).*/\
  sKey = parNodeLocations[j];\
  lastLevel = d-1;\
  foundKey = seq::maxLowerBound<ot::TreeNode>(in, sKey, idx,&i,NULL);\
  DEBUG_CHECK_FAST_MAX_LOWER_BOUND(in, sKey,idx, foundKey) \
  if ( foundKey && !((in[idx].getAnchor()==sKey.getAnchor()) && (in[idx].getFlag()&ot::TreeNode::NODE))) {\
    foundKey=false;\
  }\
  if ( foundKey ) {\
    nlist[8*i+j] = idx;\
  } else {\
    /* Can happen for own elements too. */\
    findGhost = true;\
  }\
}

#define POS_SEARCH_BLOCK(debugV1,debugV2) {\
  /* first search in default location */\
  sKey = nodeLocations[j];\
  lastLevel = d;\
  foundKey = seq::maxLowerBound<ot::TreeNode>(in, sKey, idx,&i,NULL);\
  DEBUG_CHECK_FAST_MAX_LOWER_BOUND(in, sKey,idx, foundKey) \
  if(foundKey && !( (in[idx].isAncestor(sKey) ) || (in[idx]==sKey))) {\
    foundKey=false;\
  }\
  if(foundKey) {\
    /*found somebody*/\
    if((in[idx].getAnchor()==sKey.getAnchor())&&(in[idx].getFlag() & ot::TreeNode::NODE)) {\
      /*found a node, so set it.*/\
      nlist[8*i+j] = idx;\
    } else {\
      POS_SECONDARY_SEARCH_BLOCK(debugV1,debugV2)\
    }\
  } else {\
    if(i >= m_uiElementBegin) {\
      /*The primary search for some vertex of my own element failed*/\
      /*This should only happen if the node we are looking for is*/\
      /*a post-ghost and it is hanging. When we do the check later*/\
      /*we must also test that this node is a real anchor. */\
      POS_SEARCH_DEBUG_BLOCK2\
      /*Treat this case just as if the search successfully returned*/\
      /*a node, but it turned out to be hanging*/\
      POS_SECONDARY_SEARCH_BLOCK(debugV1,debugV2)\
    } else {\
      /* this is a pre-ghost so it is normal*/\
      /*to miss some primary searches */\
      findGhost=true;\
    }\
  }\
  if ( findGhost ) {\
    findGhost=false;\
    /* need to find the ghost. */\
    /* check if idx+k is valid */\
    ot::TreeNode newKey(sKey.getX(),sKey.getY(),sKey.getZ(),\
        lastLevel,m_uiDimension,m_uiMaxDepth);\
    unsigned int k=1;\
    while ( (idx+k) < in.size() ) {\
      if ( in[idx+k].getParent() == newKey ) {\
        if ( in[idx+k].getFlag() & ot::TreeNode::NODE ) {\
          k++;\
          continue;\
        } else {\
          /* found the correct node. (hanging)*/\
          findGhost=true;\
          break;\
        }\
      }\
      if ( in[idx+k] > newKey.getDLD() ) {\
        findGhost = false;\
        break;\
      } else {\
        k++;\
      }\
    }\
    if (findGhost) {\
      nlist[8*i+j] = idx+k;\
      POS_SEARCH_DEBUG_BLOCK3\
    } else {\
      nlist[8*i+j] = m_uiLocalBufferSize;\
    }\
  }\
}

// find the eight vertices ...
//Only 0 is a special case.
for (unsigned int j = 0; j < 8; j++) {
  if(found[j]) {
    continue;
  }
  bool foundKey=false;
  bool findGhost=false;
  unsigned int idx;
  unsigned int lastLevel;
  ot::TreeNode sKey(m_uiDimension, m_uiMaxDepth);

  switch (j) {
    case 0: {
              nodeLocations[j] =
                ot::TreeNode(x, y, z, m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
              parNodeLocations[j] = 
                ot::TreeNode(parX, parY, parZ, m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
              // first search is not required since we are searching for i
              if ( !(in[i].getFlag() & ot::TreeNode::NODE ) ) {
                // if this is not a node, i.e., it is hanging
#ifdef __DEBUG_DA_NLIST__
                if ( (ch_num == 0) || (ch_num == 7) ) {
                  std::cout << "Failing for index " << i << " with child num " << ch_num << std::endl;
                  assert(false);
                }
#endif
                // All other child numbers are anchored at the parent.
                sKey = parNodeLocations[j];
                lastLevel = d-1;
                foundKey = seq::maxLowerBound<ot::TreeNode> (in, sKey, idx,NULL,&i);
#ifdef __DEBUG_DA_NLIST__
                CHECK_FAST_MAX_LOWER_BOUND(in, sKey,idx, foundKey) 
#endif
                  if ( foundKey && 
                      (!( (in[idx].getAnchor() == sKey.getAnchor()) && (in[idx].getFlag() & ot::TreeNode::NODE) )) ) {
                    foundKey=false;
                  }
                if ( !foundKey ) {
#ifdef __DEBUG_DA_NLIST__
                  // should only happen for ghosts ...
                  if(i >= m_uiElementBegin) {
                    std::cout<<m_iRankActive<<" i = "<<i<<" preGhostElemEnd "<<m_uiPreGhostElementSize
                      <<" elemBeg: "<<m_uiElementBegin<<" elemEnd: "<<m_uiElementEnd<<std::endl;
                  }
                  assert (i < m_uiElementBegin);
#endif
                  // for node zero, simply default to i.
                  nlist[8*i+j] = i;
                } else {
                  nlist[8*i+j] = idx;
                }
              } else {
                nlist[8*i+j] = i;
              }
              break;
            }
    case 1: {
              nodeLocations[j] = 
                ot::TreeNode(x+sz, y, z, m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
              parNodeLocations[j] = 
                ot::TreeNode(parX+(sz<<1u), parY, parZ, m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
              POS_SEARCH_BLOCK(1,6) 
                break;
            }
    case 2: {
              nodeLocations[j] = 
                ot::TreeNode(x, y+sz, z, m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
              parNodeLocations[j] = 
                ot::TreeNode(parX, parY+(sz<<1u), parZ, m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
              POS_SEARCH_BLOCK(2,5) 
                break;
            }
    case 3: {
              nodeLocations[j] = 
                ot::TreeNode(x+sz, y+sz, z, m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
              parNodeLocations[j] = 
                ot::TreeNode(parX+(sz<<1u), parY+(sz<<1u), parZ, m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
              POS_SEARCH_BLOCK(3,4) 
                break;
            }
    case 4: {
              nodeLocations[j] = 
                ot::TreeNode(x, y, z+sz, m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
              parNodeLocations[j] = 
                ot::TreeNode(parX, parY, parZ+(sz<<1u), m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
              POS_SEARCH_BLOCK(4,3) 
                break;
            }
    case 5: {
              nodeLocations[j] = 
                ot::TreeNode(x+sz, y, z+sz, m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
              parNodeLocations[j] = 
                ot::TreeNode(parX+(sz<<1u), parY, parZ+(sz<<1u), m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
              POS_SEARCH_BLOCK(5,2) 
                break;
            }
    case 6: {
              nodeLocations[j] = 
                ot::TreeNode(x, y+sz, z+sz, m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
              parNodeLocations[j] = 
                ot::TreeNode(parX, parY+(sz<<1u), parZ+(sz<<1u), m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
              POS_SEARCH_BLOCK(6,1) 
                break;
            }
    case 7: {
              nodeLocations[j] = 
                ot::TreeNode(x+sz, y+sz, z+sz, m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
              parNodeLocations[j] = 
                ot::TreeNode(parX+(sz<<1u), parY+(sz<<1u), parZ+(sz<<1u), m_uiMaxDepth, m_uiDimension, m_uiMaxDepth);
              POS_SEARCH_BLOCK(7,0) 
                break;
            }
    default: {
               std::cerr << RED<<"Wrong node number in " << __func__ <<NRM<< std::endl;
               assert(false);
               break;
             }
  } // end switch
} // end loop over the 8 vertices of this element...

#undef POS_SEARCH_BLOCK
#undef POS_SECONDARY_SEARCH_BLOCK
#undef POS_SEARCH_DEBUG_BLOCK1
#undef POS_SEARCH_DEBUG_BLOCK2
#undef POS_SEARCH_DEBUG_BLOCK3

// FINISHED Searching for Node Indices ...

//Ensure that the anchor of the local element is not pointing to ghost.
//This can happen only if the octant in question is a singular block
//and its anchor is hanging and the 0th child of its parent is sitting on a different processor.
//BlockPart should have detected this case and prevented this.
#ifdef __DEBUG_DA_NLIST__
if ( (i >= m_uiElementBegin) && (i < m_uiElementEnd) &&
    ( (nlist[8*i] < m_uiElementBegin) || (nlist[8*i] >= m_uiPostGhostBegin) ) ) {
  std::cout << "At index " << i << " anchor is at " << nlist[8*i] <<
    " where elemBegin is at " << m_uiElementBegin << " and postGhBegin is "
    << m_uiPostGhostBegin << std::endl;
  std::cout << "RANK is " << m_iRankActive << std::endl;
  std::cout << "NList is ";
  for (unsigned int j=0; j<8; j++) {
    std::cout << nlist[8*i+j] << ", ";
  }
  std::cout << std::endl;
  std::cout<<m_iRankActive<<" failingOct: "<<in[i]<<std::endl<<" it's parent: "<<in[i].getParent()<<std::endl;
  if( nlist[8*i] < in.size() ) {
    std::cout<<m_iRankActive<<" failingOct's anchor is actually mapped to: "<<in[nlist[8*i]]<<std::endl;  
  }
  assert(false);
}
#endif

#ifdef __DEBUG_DA_NLIST__
//Check if you sent yourself apriori (Second Ring). 
if( (i >= m_uiElementBegin) && (i < m_uiPostGhostBegin) ) {
  for(unsigned int j=0; j< 8; j++) {
    if(nlist[8*i + j] < m_uiElementBegin) {
      std::cout<<m_iRankActive<<" Trying to send yourself as a Post Ghost ELEMENT for  i = "
        <<i<<" j = "<<j<<std::endl;
      assert(false);
    }
    if( (nlist[8*i+j] >= m_uiPostGhostBegin) && (nlist[8*i+j] < m_uiLocalBufferSize) ) {
      TreeNode tmpToSend = in[nlist[8*i+j]];
      tmpToSend.setWeight(i);
      checkSecondRing.push_back(tmpToSend);
    }
  }//end for j
}
#endif

// compute the hanging node mask for this element ...
unsigned char _mask=0;
unsigned int numOutOfBounds = 0;
for ( unsigned int j=0; j<8; j++) {
  if ( ( nlist[8*i+j] < m_uiElementBegin ) || (nlist[8*i+j] >= m_uiPostGhostBegin) ) {
    numOutOfBounds++;
  }

  if(nlist[8*i+j] >= m_uiLocalBufferSize) {
    //Skip setting mask for this vtx.
    continue;
  }

  // check if any of the nodes is hanging ...
  unsigned int _x,_y,_z, _d;
  _x = in[nlist[8*i+j]].getX(); 
  _y = in[nlist[8*i+j]].getY(); 
  _z = in[nlist[8*i+j]].getZ(); 
  _d = in[nlist[8*i+j]].getLevel();
  if ( !(in[nlist[8*i+j]].getFlag() & ot::TreeNode::NODE ) ) {
    // std::cout << "For i=" << i << " looking at parent" << std::endl;
    _x  = ( ( _x >> ( m_uiMaxDepth - _d + 1 ) ) << ( m_uiMaxDepth - _d + 1 ) ); 
    _y  = ( ( _y >> ( m_uiMaxDepth - _d + 1 ) ) << ( m_uiMaxDepth - _d + 1 ) );
    _z  = ( ( _z >> ( m_uiMaxDepth - _d + 1 ) ) << ( m_uiMaxDepth - _d + 1 ) );
  }

  switch (j) {
    case 0:
      if ( !(in[i].getFlag() & ot::TreeNode::NODE ) ) {
        _mask |= (1 << j);
      }
      break;
    case 1:
      // look at the anchor of the +x neighbor
      // if ( (x+sz) != _x ) {
      if ( ( (x+sz) != _x ) || ( y != _y ) || ( z != _z ) ) {
        _mask |= (1 << j);
      }
      break;
    case 2:
      // look at the anchor of the +y neighbor
      // if ( (y+sz) != _y ) {
      if ( ( x != _x ) || ( (y+sz) != _y ) || ( z != _z ) ) {
        _mask |= (1 << j);
      }
      break;
    case 3:
      // look at both x and y anchors ...
      // if ( ( (x+sz) != _x ) || ( (y+sz) != _y )  ) {
      if ( ( (x+sz) != _x ) || ( (y+sz) != _y ) || ( z != _z ) ) {
        _mask |= (1 << j);
      }
      break;
    case 4:
      // look at z anchor
      // if ( (z+sz) != _z ) {
      if ( ( x != _x ) || ( y != _y ) || ( (z+sz) != _z ) ) {
        _mask |= (1 << j);
      }
      break;
    case 5:
      // look at +z,+x
      // if ( ( (x+sz) != _x ) || ( (z+sz) != _z )  ) {
      if ( ( (x+sz) != _x ) || ( y != _y ) || ( (z+sz) != _z ) ) {
        _mask |= (1 << j);
      }
      break;
    case 6:
      // look at +z, +y
      // if ( ( (z+sz) != _z ) || ( (y+sz) != _y )  ) {
      if ( ( x != _x ) || ( (y+sz) != _y ) || ( (z+sz) != _z ) ) {
        _mask |= (1 << j);
      }
      break;
    case 7:
      if ( ( (x+sz) != _x ) || ( (y+sz) != _y ) || ( (z+sz) != _z ) ) {
        _mask |= (1 << j);
      }
      break;
  }//end switch-case
  }//end for j

  // store the mask ...
  if (numOutOfBounds == 8) {
    //This does not even have one writable node and hence this is not an element.
    _mask = ot::DA_FLAGS::FOREIGN;
  } else if (numOutOfBounds) {
    //A Dependent Element, is one which has atleast one writable node and atleast one ghosted node.
    in[i].orFlag(ot::DA_FLAGS::DEP_ELEM);        
  }

  //Prepare for the Ugly portion...
  //Sometimes, we might find a primary key but it could turn out to be hanging
  //and then we may not find the secondary key. In such cases, we still
  //generate both the priimary and secondary keys and search for both in the
  //following second ring correction phase. This might seem like an overkill
  //but this situation occurs only for pre-ghosts and singular blocks. Suppose,
  //the octant A searches for a secondary key B and does not find it and
  //suppose A is not a pre-ghost (i.e. this processor owns A). Since B is one
  //of the vertices of A's parent this means the sibling of A that shares the
  //vertex B with A's parent is not on the same processor as A (since all
  //direct vertices which are not hanging are communicated apriori and B can't
  //be hanging). Hence, A is singular
  for (int j=0; j<8; j++) {
    if (nlist[8*i+j] >= m_uiLocalBufferSize) {
      if (_mask != ot::DA_FLAGS::FOREIGN) {
#ifdef __DEBUG_DA_NLIST__
        assert(j);
#endif
        // add to list ...
        nodeLocations[j].setWeight(m_iRankActive);
        parNodeLocations[j].setWeight(m_iRankActive);
        primaryKeys.push_back(nodeLocations[j]);
        secondaryKeys.push_back(parNodeLocations[j]);
        // correct lookUp Table ...
        nlist[8*i+j] = m_uiLocalBufferSize +
          static_cast<unsigned int>(extraAtEnd.size());
        nodeLocations[j].setWeight(i);
        extraAtEnd.push_back(nodeLocations[j]);
      }//end if foreign
    }//end if invalid
  }//end for j

  // Store the hanging node mask
  m_ucpLutMasks[2*i+1] = _mask;
  } // end for i: All elements in this set 

#ifdef __DEBUG_DA_NLIST__
  MPI_Barrier(m_mpiCommActive);
  if(!m_iRankActive) {
    std::cout<<std::endl;
    std::cout<<"Finished Elemental Loop for Set# "<<numFullLoopCtr<<std::endl;
    std::cout<<std::endl;
  }
  MPI_Barrier(m_mpiCommActive);
#endif

#ifdef __DEBUG_DA_NLIST__
  MPI_Barrier(m_mpiCommActive);
  seq::makeVectorUnique<ot::TreeNode>(chkMissedPrimary, false);
  int* chkMissedPrimarySendCounts = new int[m_iNpesActive];
  for(int i = 0; i < m_iNpesActive; i++) {
    chkMissedPrimarySendCounts[i] = 0; 
  }//end for i

  for(unsigned int i = 0; i < chkMissedPrimary.size(); i++) {
    //maxLB returns the last index in a sorted array such that a[ind] <= key and  a[index +1] > key
    unsigned int idx;
    bool found = seq::maxLowerBound<TreeNode >(m_tnMinAllBlocks, chkMissedPrimary[i], idx, NULL, NULL);
    assert(found);
    //missed keys must be post-ghosts
    assert(idx > m_iRankActive);
    chkMissedPrimarySendCounts[idx]++; 
  }//end for i

  int* chkMissedPrimaryRecvCounts = new int[m_iNpesActive];
  par::Mpi_Alltoall<int>(chkMissedPrimarySendCounts, chkMissedPrimaryRecvCounts, 1, m_mpiCommActive);

  int* chkMissedPrimarySendOffsets = new int[m_iNpesActive];
  int* chkMissedPrimaryRecvOffsets = new int[m_iNpesActive];
  chkMissedPrimarySendOffsets[0] = 0;
  chkMissedPrimaryRecvOffsets[0] = 0;
  for(int i = 1; i < m_iNpesActive; i++) {
    chkMissedPrimarySendOffsets[i] = chkMissedPrimarySendOffsets[i-1] + chkMissedPrimarySendCounts[i-1] ;
    chkMissedPrimaryRecvOffsets[i] = chkMissedPrimaryRecvOffsets[i-1] + chkMissedPrimaryRecvCounts[i-1] ;
  }//end for i
  unsigned int chkMissedPrimaryRecvSz = chkMissedPrimaryRecvOffsets[m_iNpesActive - 1] +
    chkMissedPrimaryRecvCounts[m_iNpesActive - 1];
  std::vector<ot::TreeNode> chkMissedPrimaryRecvBuffer(chkMissedPrimaryRecvSz);
  ot::TreeNode* chkMissedPrimarySendPtr = NULL;
  ot::TreeNode* chkMissedPrimaryRecvPtr = NULL;
  if(!(chkMissedPrimary.empty())) {
    chkMissedPrimarySendPtr = (&(*(chkMissedPrimary.begin())));
  }
  if(!(chkMissedPrimaryRecvBuffer.empty())) {
    chkMissedPrimaryRecvPtr = (&(*(chkMissedPrimaryRecvBuffer.begin())));
  }
  par::Mpi_Alltoallv_sparse<ot::TreeNode>( chkMissedPrimarySendPtr, chkMissedPrimarySendCounts,
      chkMissedPrimarySendOffsets, chkMissedPrimaryRecvPtr,
      chkMissedPrimaryRecvCounts, chkMissedPrimaryRecvOffsets, m_mpiCommActive);

  for(unsigned int i = 0; i < chkMissedPrimaryRecvBuffer.size(); i++) {
    unsigned int idx;
    bool found = seq::maxLowerBound<TreeNode >(in, chkMissedPrimaryRecvBuffer[i], idx, NULL, NULL);
    assert(found);
    assert( (idx >= m_uiElementBegin) && (idx < m_uiPostGhostBegin) );
    //Although this is a hanging node this is also some anchor
    assert(in[idx].getAnchor() == chkMissedPrimaryRecvBuffer[i].getAnchor());
    assert( !(in[idx].getFlag() & ot::TreeNode::NODE) );
  }//end for i

  delete [] chkMissedPrimarySendCounts;
  delete [] chkMissedPrimaryRecvCounts;
  delete [] chkMissedPrimarySendOffsets;
  delete [] chkMissedPrimaryRecvOffsets;
  chkMissedPrimaryRecvBuffer.clear();
  MPI_Barrier(m_mpiCommActive);
  if(!m_iRankActive) {
    std::cout<<std::endl;
    std::cout<<"Passed Test for Missed Primary "<<numFullLoopCtr<<std::endl;
    std::cout<<std::endl;
  }
  MPI_Barrier(m_mpiCommActive);
#endif

  /*
     NOW, The Ugly Parallel Book-keeping and Corrections for the Missing Entries in the Previous Step...
     */

  // ~~~~~~~~~~~~~~~~~~ SECONDARY ~~~~~~~~~~~~~~~~~~~~~~~
  std::vector<unsigned int>               ScndScatterMap;

  std::vector<unsigned int>               ScndSendProcs;
  std::vector<unsigned int>               ScndSendCounts;

  std::vector<unsigned int>               ScndRecvProcs;
  std::vector<unsigned int>               ScndRecvCounts;

  //~~~~~~~~~~~~~~~~~~~~~~~~ Search for Failed Keys ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  /*
     1. Each octant which has a missing entry in its node-list would have generated a primary and secondary key for that entry.
     This is the last for-j loop within the main i-loop.
     2. Multiple octants might want the same entry. Hence, we make these keys unique before sending these requests
     to the respective processors.
     3. However, we need to map the results back to the octants which generated these keys.
     So we keep a back-up of the original list of keys. This list will have duplicates.
     4. Send the unique set of Primary and Secondary keys computed above to the processors which control those domains.
     5. Each processor recieves the keys that lie within its domain.
     Now, multiple processors might have requested for the same key. So, we make the list of keys recieved unique.
     We need to map them back to the processors that requested them. Hence, we store the pair of keys and 
     a list of processors that requested that key. This is done using the NodeAndRanks class. The sort order
     for this class is defined only on TreeNode.
     6. Each processor performs a local search with this set of unique keys. Note only non-hanging nodes are considered as matches.
     7. The results are returned to the processors that requested the respective keys.
     8. The results are matched with the octants that generated the keys. 
     9. If a primary key returned a positive result, then it is used else the secondary key is used. 
     10. There is a special for the secondary key for the nlist of a pre-ghost element. It is possible that while the primary key was NOT recieved during the a-priori communication, but the secondary key was recieved. In such situations if the secondary key is to be selected, then we must select the copy that was recieved during the a-priori comm and not the one got from this second-ring correction. This must be done in order to prevent having duplicate elements in our local buffer.
     11. Each of the octants that were selected is given an unique id. This will be used to fix the nlist.
     12. The hanging masks need to be corrected as well.
     13. Since, the primary key could be owned by one processor and the secondary key by another. Only the processor which owns the octant that generated
     the keys can decide what key is actually picked. The decision is then communicated to the processors that own the keys.
     14. A processor might be sending some primary results and some secondary results to the same processor. So both the sets must be merged and
     a single scattermap is built. This is the secondary scattermap. This is from the primary scattermap, which is built in the constructor 
     before entering this function (BuildNodeList).
     15. Finally, the same processor might have sent some of its elements in the first ring (apriori comm inside the constructor) and it
     might send some of its elements in the second ring (below). When the actual data is sent, we want the two sets to come together.
     We want both the first set of octants and the second set of octants to be merged and to be sorted.
     16. The first set of octants and the second set of octants have to be merged and sorted. All the second set of octants are marked as FOREIGNs. So, we don't loop over them in the MatVec and we don't build their LUTs. They are simply place holders for ghost values.
     17. The old LUT has to be re-mapped to the new LUT, since the indices will change after step 15.
     18. Similarly, the scattermaps must be merged and corrected to point to the new indices.
     19. We can have missing nodes both for pre-ghosts as well as own elements. Hence, this second ring correction is done inside the main outer loop (numFullLoopCtr).
     */

  //~~~~~~~~~~~Correction for Second Ring Begins~~~~~~~~~~~~~~~~~~//
  //First Get the min and max from each processor.

  ot::TreeNode rootNode(m_uiDimension,m_uiMaxDepth);

  assert(!in.empty());
  assert(m_uiElementBegin < in.size());
  assert(m_uiPostGhostBegin >= 1);
  assert((m_uiPostGhostBegin-1) < in.size());

  std::vector<ot::TreeNode> failedPrimaryKeys = primaryKeys;
  std::vector<ot::TreeNode> failedSecondaryKeys = secondaryKeys;

  //Sort and Make Unique
  seq::makeVectorUnique<ot::TreeNode>(failedPrimaryKeys, false);
  seq::makeVectorUnique<ot::TreeNode>(failedSecondaryKeys, false);

#ifdef __DEBUG_DA_NLIST__
  MPI_Barrier(m_mpiCommActive);

  if(!m_iRankActive) {
    std::cout<<std::endl;
    std::cout<<"Made Failed Keys Unique. Finding Partition Next."<<std::endl;
    std::cout<<std::endl;
  }
  MPI_Barrier(m_mpiCommActive);
#endif

  unsigned int PrimaryKeysSz = static_cast<unsigned int>(failedPrimaryKeys.size());
  unsigned int SecondaryKeysSz = static_cast<unsigned int>(failedSecondaryKeys.size());

  //Now determine the processors which own these keys.
  unsigned int *partPrimary = NULL;
  if(PrimaryKeysSz) {
    partPrimary = new unsigned int[PrimaryKeysSz];    
  }
  unsigned int *partSecondary = NULL;
  if(SecondaryKeysSz) {
    partSecondary = new unsigned int[SecondaryKeysSz];    
  }

  for (unsigned int i=0; i<PrimaryKeysSz; i++) {
    unsigned int idx;
    //maxLB returns the last index in a sorted array such that a[ind] <= key and  a[index +1] > key
    bool found = seq::maxLowerBound<TreeNode >(m_tnMinAllBlocks, failedPrimaryKeys[i], idx, NULL, NULL);
    assert(found);
    partPrimary[i] = idx;
  }//end for i

  for (unsigned int i=0; i<SecondaryKeysSz; i++) {
    unsigned int idx;
    //maxLB returns the last index in a sorted array such that a[ind] <= key and  a[index +1] > key
    bool found = seq::maxLowerBound<TreeNode >(m_tnMinAllBlocks, failedSecondaryKeys[i], idx, NULL, NULL);
    assert(found);
    partSecondary[i] = idx;
  }//end for i

  int *numKeysSendP = new int[m_iNpesActive];
  int *numKeysSendS = new int[m_iNpesActive];
  int *numKeysRecvP = new int[m_iNpesActive];    
  int *numKeysRecvS = new int[m_iNpesActive];

  for (int i=0; i<m_iNpesActive; i++) {
    numKeysSendP[i] = 0;
    numKeysSendS[i] = 0;
  }
  // calculate the number of keys to send ...
  for (unsigned int i=0; i<PrimaryKeysSz; i++) {
    numKeysSendP[partPrimary[i]]++;      
  }
  for (unsigned int i=0; i<SecondaryKeysSz; i++) {
    numKeysSendS[partSecondary[i]]++;
  }    

#ifdef __DEBUG_DA_NLIST__
  MPI_Barrier(m_mpiCommActive);
  if(!m_iRankActive) {
    std::cout<<std::endl;
    std::cout<<"First ALL2ALL for Second Ring..."<<std::endl;
    std::cout<<std::endl;
  }
  MPI_Barrier(m_mpiCommActive);
#endif


  // Now do an All2All to get inumKeysRecv
  PROF_BUILD_NLIST_COMM_BEGIN

    par::Mpi_Alltoall<int>(numKeysSendP, numKeysRecvP, 1, m_mpiCommActive);
  par::Mpi_Alltoall<int>(numKeysSendS, numKeysRecvS, 1, m_mpiCommActive);

  PROF_BUILD_NLIST_COMM_END

    // Now create sendK
    int *sendOffsetsP = new int[m_iNpesActive]; sendOffsetsP[0] = 0;
  int *recvOffsetsP = new int[m_iNpesActive]; recvOffsetsP[0] = 0;
  int *numKeysTmpP = new int[m_iNpesActive]; numKeysTmpP[0] = 0; 

  int *sendOffsetsS = new int[m_iNpesActive]; sendOffsetsS[0] = 0;
  int *recvOffsetsS = new int[m_iNpesActive]; recvOffsetsS[0] = 0;
  int *numKeysTmpS = new int[m_iNpesActive]; numKeysTmpS[0] = 0; 

  // compute offsets ...
  for (int i = 1; i < m_iNpesActive; i++) {
    sendOffsetsP[i] = sendOffsetsP[i-1] + numKeysSendP[i-1];
    recvOffsetsP[i] = recvOffsetsP[i-1] + numKeysRecvP[i-1];
    numKeysTmpP[i] = 0; 

    sendOffsetsS[i] = sendOffsetsS[i-1] + numKeysSendS[i-1];
    recvOffsetsS[i] = recvOffsetsS[i-1] + numKeysRecvS[i-1];
    numKeysTmpS[i] = 0; 
  }

  // create the send and recv buffers ...
  std::vector<ot::TreeNode> sendKp (PrimaryKeysSz);
  std::vector<ot::TreeNode> recvKp (recvOffsetsP[m_iNpesActive-1] + numKeysRecvP[m_iNpesActive-1]);

  std::vector<ot::TreeNode> sendKs (SecondaryKeysSz);
  std::vector<ot::TreeNode> recvKs (recvOffsetsS[m_iNpesActive-1] + numKeysRecvS[m_iNpesActive-1]);

  for (unsigned int i=0; i< PrimaryKeysSz; i++) {
    unsigned int ni = numKeysTmpP[partPrimary[i]];
    numKeysTmpP[partPrimary[i]]++;
    // set entry ...
    sendKp[sendOffsetsP[partPrimary[i]] + ni] = failedPrimaryKeys[i];      
  } 

  for (unsigned int i=0; i< SecondaryKeysSz; i++) {
    unsigned int ni = numKeysTmpS[partSecondary[i]];
    numKeysTmpS[partSecondary[i]]++;
    // set entry ...
    sendKs[sendOffsetsS[partSecondary[i]] + ni] = failedSecondaryKeys[i];       
  } 

  failedPrimaryKeys.clear();
  failedSecondaryKeys.clear();

  if(partPrimary) {
    delete [] partPrimary;
    partPrimary = NULL;
  }

  if(partSecondary) {
    delete [] partSecondary;
    partSecondary = NULL;
  }

  if(numKeysTmpP) {
    delete [] numKeysTmpP;    
    numKeysTmpP = NULL;
  }

  if(numKeysTmpS) {
    delete [] numKeysTmpS;
    numKeysTmpS = NULL;
  }

  ot::TreeNode* sendKpPtr = NULL;
  ot::TreeNode* recvKpPtr = NULL;
  ot::TreeNode* sendKsPtr = NULL;
  ot::TreeNode* recvKsPtr = NULL;
  if(!sendKp.empty()) {
    sendKpPtr = &(*(sendKp.begin()));
  }
  if(!recvKp.empty()) {
    recvKpPtr = &(*(recvKp.begin()));
  }
  if(!sendKs.empty()) {
    sendKsPtr = &(*(sendKs.begin()));
  }
  if(!recvKs.empty()) {
    recvKsPtr = &(*(recvKs.begin()));
  }

  PROF_BUILD_NLIST_COMM_BEGIN

    par::Mpi_Alltoallv_sparse<ot::TreeNode>( sendKpPtr, numKeysSendP, sendOffsetsP,
        recvKpPtr, numKeysRecvP, recvOffsetsP, m_mpiCommActive);

  par::Mpi_Alltoallv_sparse<ot::TreeNode>( sendKsPtr, numKeysSendS, sendOffsetsS,
      recvKsPtr, numKeysRecvS, recvOffsetsS, m_mpiCommActive);

  PROF_BUILD_NLIST_COMM_END

    sendKp.clear();
  sendKs.clear();

  delete [] sendOffsetsP;
  sendOffsetsP = NULL;

  delete [] recvOffsetsP;
  recvOffsetsP = NULL;

  delete [] numKeysSendP;
  numKeysSendP = NULL;

  delete [] numKeysRecvP;
  numKeysRecvP = NULL;

  delete [] sendOffsetsS;
  sendOffsetsS = NULL;

  delete [] recvOffsetsS;
  recvOffsetsS = NULL;

  delete [] numKeysSendS;
  numKeysSendS = NULL;

  delete [] numKeysRecvS;
  numKeysRecvS = NULL;

  std::vector<ot::NodeAndRanks> recvK2P;
  std::vector<ot::NodeAndRanks> recvK2S;
  //recvKp and recvKs are NOT sorted and NOT unique.
  //First merge recvK into recvK2.
  for (int i = 0; i < recvKp.size(); i++) {
    ot::NodeAndRanks tmp;
    tmp.node = recvKp[i];
    tmp.ranks.push_back(recvKp[i].getWeight());
    recvK2P.push_back(tmp);
  }

  for (int i=0;i<recvKs.size();i++) {
    ot::NodeAndRanks tmp;
    tmp.node = recvKs[i];
    tmp.ranks.push_back(recvKs[i].getWeight());
    recvK2S.push_back(tmp);
  }

  recvKp.clear();  
  recvKs.clear();  

  std::sort(recvK2P.begin(),recvK2P.end());
  std::sort(recvK2S.begin(),recvK2S.end());

  //Make recvK2P Unique and concatenate ranks.
  if (recvK2P.size() >= 2) {
    std::vector<ot::NodeAndRanks> tmp(recvK2P.size());
    tmp[0] = recvK2P[0];
    unsigned int tmpSize=1;
    for (unsigned int i=1;i<recvK2P.size();i++) {
      if (tmp[tmpSize-1] != recvK2P[i]) {
        //new entry
        tmp[tmpSize] = recvK2P[i];
        tmpSize++;
      } else {
        tmp[tmpSize-1].ranks.push_back(recvK2P[i].ranks[0]);
      }
    }//end for
    tmp.resize(tmpSize);
    for (unsigned int i=0;i<tmpSize;i++) {
      seq::makeVectorUnique<int>(tmp[i].ranks,false);
    }
    recvK2P = tmp; 
    tmp.clear();
  }

  //Make recvK2S Unique and concatenate ranks.
  if (recvK2S.size() >= 2) {
    std::vector<ot::NodeAndRanks> tmp(recvK2S.size());
    tmp[0] = recvK2S[0];
    unsigned int tmpSize=1;
    for (unsigned int i=1;i<recvK2S.size();i++) {
      if (tmp[tmpSize-1] != recvK2S[i]) {
        //new entry
        tmp[tmpSize] = recvK2S[i];
        tmpSize++;
      } else {
        tmp[tmpSize-1].ranks.push_back(recvK2S[i].ranks[0]);
      }
    }//end for
    tmp.resize(tmpSize);
    for (unsigned int i=0;i<tmpSize;i++) {
      seq::makeVectorUnique<int>(tmp[i].ranks,false);
    }
    recvK2S = tmp; 
    tmp.clear();
  }

  //Local Search and update sendNodes and sendCnt.		
  std::vector<std::vector<ot::TreeNode> > sendNodesP(m_iNpesActive);
  std::vector<std::vector<ot::TreeNode> > sendNodesS(m_iNpesActive);
  //int is necessary here. Set to -1 later.
  std::vector<std::vector<int> > idxP(m_iNpesActive);
  std::vector<std::vector<int> > idxS(m_iNpesActive);

#ifdef __DEBUG_DA_NLIST__
  MPI_Barrier(m_mpiCommActive);
  if(!m_iRankActive) {
    std::cout<<std::endl;
    std::cout<<"Starting Local Search for Second Ring..."<<std::endl;
    std::cout<<std::endl;
  }
  MPI_Barrier(m_mpiCommActive);
#endif

  //in is sorted and unique and linear and recvK2 is sorted and unique and linear.    
  for (unsigned int i=0; i<recvK2P.size();i++) {
    unsigned int idx;    
    bool found = seq::maxLowerBound<ot::TreeNode >(in, recvK2P[i].node, idx,NULL,NULL);

    if (found) {
#ifdef __DEBUG_DA_NLIST__
      assert( (in[idx].isAncestor(recvK2P[i].node)) || (in[idx] == (recvK2P[i].node)) );
      assert( (idx >= m_uiElementBegin) && (idx < m_uiPostGhostBegin) );
#endif

      if ( (in[idx].getAnchor() != recvK2P[i].node.getAnchor()) || 
          (!(in[idx].getFlag() & ot::TreeNode::NODE)) ) {
        found = false;
      }
    }

    //Send the result to all the processors that want it.
    for (int j=0; j < recvK2P[i].ranks.size(); j++) {
      if ( found ) {
        //Send in[idx];
        sendNodesP[recvK2P[i].ranks[j]].push_back(in[idx]);          
        idxP[recvK2P[i].ranks[j]].push_back(idx);
      } else {
        //Send rootNode;  
        sendNodesP[recvK2P[i].ranks[j]].push_back(rootNode);          
        idxP[recvK2P[i].ranks[j]].push_back(-1);
      }
      sendNodesP[recvK2P[i].ranks[j]][sendNodesP[recvK2P[i].ranks[j]].size()-1].setWeight(m_iRankActive);
    }//end for j
  }//end for i		

  for (unsigned int i=0; i<recvK2S.size();i++) {
    unsigned int idx;    
    bool found = seq::maxLowerBound<ot::TreeNode >(in, recvK2S[i].node, idx, NULL, NULL);
    if (found) {
#ifdef __DEBUG_DA_NLIST__
      assert( (in[idx].isAncestor(recvK2S[i].node)) || (in[idx] == (recvK2S[i].node)) );
      assert( (idx >= m_uiElementBegin) && (idx < m_uiPostGhostBegin) );
#endif
      if ( (in[idx].getAnchor() != recvK2S[i].node.getAnchor()) || (!(in[idx].getFlag() & ot::TreeNode::NODE)) ) {
        found = false;
      }
    }
    //Send the result to all the processors that want it.
    for (int j=0; j < recvK2S[i].ranks.size(); j++) {
      if ( found ) {
        //Send in[idx];
        sendNodesS[recvK2S[i].ranks[j]].push_back(in[idx]);
        idxS[recvK2S[i].ranks[j]].push_back(idx);
      } else {
        //Send rootNode;  
        sendNodesS[recvK2S[i].ranks[j]].push_back(rootNode);
        idxS[recvK2S[i].ranks[j]].push_back(-1);
      }
      sendNodesS[recvK2S[i].ranks[j]][sendNodesS[recvK2S[i].ranks[j]].size()-1].setWeight(m_iRankActive);
    }
  }//end for i		

  recvK2P.clear(); 
  recvK2S.clear(); 

  for (unsigned int i = 0; i < m_iNpesActive; i++) {
    seq::makeVectorUnique<TreeNode>(sendNodesP[i],false);
    seq::makeVectorUnique<TreeNode>(sendNodesS[i],false);
    seq::makeVectorUnique<int>(idxP[i],false);
    seq::makeVectorUnique<int>(idxS[i],false);
  }

  numKeysSendP = new int[m_iNpesActive];
  numKeysSendS = new int[m_iNpesActive];
  numKeysRecvP = new int[m_iNpesActive];    
  numKeysRecvS = new int[m_iNpesActive];

  for (int i=0; i<m_iNpesActive; i++) {
    numKeysSendP[i] = static_cast<int>(sendNodesP[i].size());
    numKeysSendS[i] = static_cast<int>(sendNodesS[i].size());    
  }        

#ifdef __DEBUG_DA_NLIST__
  MPI_Barrier(m_mpiCommActive);
  if(!m_iRankActive) {
    std::cout<<std::endl;
    std::cout<<"Sending Results after Local Search for Second Ring..."<<std::endl;
    std::cout<<std::endl;
  }
  MPI_Barrier(m_mpiCommActive);
#endif

  // Now do an All2All to get inumKeysRecv
  PROF_BUILD_NLIST_COMM_BEGIN

    par::Mpi_Alltoall<int>(numKeysSendP, numKeysRecvP, 1, m_mpiCommActive);
  par::Mpi_Alltoall<int>(numKeysSendS, numKeysRecvS, 1, m_mpiCommActive);    

  PROF_BUILD_NLIST_COMM_END

    // Now create sendK
    sendOffsetsP = new int[m_iNpesActive]; sendOffsetsP[0] = 0;
  sendOffsetsS = new int[m_iNpesActive]; sendOffsetsS[0] = 0;
  recvOffsetsP = new int[m_iNpesActive]; recvOffsetsP[0] = 0;   
  recvOffsetsS = new int[m_iNpesActive]; recvOffsetsS[0] = 0;

  // compute offsets ...
  for (int i=1; i<m_iNpesActive; i++) {
    sendOffsetsP[i] = sendOffsetsP[i-1] + numKeysSendP[i-1];
    recvOffsetsP[i] = recvOffsetsP[i-1] + numKeysRecvP[i-1];

    sendOffsetsS[i] = sendOffsetsS[i-1] + numKeysSendS[i-1];
    recvOffsetsS[i] = recvOffsetsS[i-1] + numKeysRecvS[i-1];    
  }

  // create the send and recv buffers ...
  sendKp.resize(sendOffsetsP[m_iNpesActive-1] + numKeysSendP[m_iNpesActive-1]);
  recvKp.resize(recvOffsetsP[m_iNpesActive-1] + numKeysRecvP[m_iNpesActive-1]);

  sendKs.resize(sendOffsetsS[m_iNpesActive-1] + numKeysSendS[m_iNpesActive-1]);
  recvKs.resize(recvOffsetsS[m_iNpesActive-1] + numKeysRecvS[m_iNpesActive-1]);

  std::vector<int> idxSendKp(sendOffsetsP[m_iNpesActive-1] + numKeysSendP[m_iNpesActive-1]);
  std::vector<int> idxSendKs(sendOffsetsS[m_iNpesActive-1] + numKeysSendS[m_iNpesActive-1]);

  for (unsigned int i = 0; i < m_iNpesActive; i++) {
    for (unsigned int j = 0; j < numKeysSendP[i]; j++) {
      // set entry ...
      sendKp[sendOffsetsP[i] + j] = sendNodesP[i][j];       
      idxSendKp[sendOffsetsP[i] + j] = idxP[i][j];       

#ifdef __DEBUG_DA_NLIST__
      assert( ((sendKp[sendOffsetsP[i] + j] == rootNode) && (idxSendKp[sendOffsetsP[i] + j] == -1)) ||
          (in[idxSendKp[sendOffsetsP[i] + j]] == sendKp[sendOffsetsP[i] + j]) );
#endif
    }
    for (unsigned int j = 0; j < numKeysSendS[i]; j++) {
      // set entry ...
      sendKs[sendOffsetsS[i] + j] = sendNodesS[i][j];    
      idxSendKs[sendOffsetsS[i] + j] = idxS[i][j];       

#ifdef __DEBUG_DA_NLIST__
      assert( ((sendKs[sendOffsetsS[i] + j] == rootNode) && (idxSendKs[sendOffsetsS[i] + j] == -1)) ||
          (in[idxSendKs[sendOffsetsS[i] + j]] == sendKs[sendOffsetsS[i] + j]) );
#endif
    }
  } 

  sendNodesP.clear();
  sendNodesS.clear();

  idxP.clear();
  idxS.clear();

  sendKpPtr = NULL;
  recvKpPtr = NULL;
  sendKsPtr = NULL;
  recvKsPtr = NULL;
  if(!sendKp.empty()) {
    sendKpPtr = &(*(sendKp.begin()));
  }
  if(!sendKs.empty()) {
    sendKsPtr = &(*(sendKs.begin()));
  }
  if(!recvKp.empty()) {
    recvKpPtr = &(*(recvKp.begin()));
  }
  if(!recvKs.empty()) {
    recvKsPtr = &(*(recvKs.begin()));
  }

  PROF_BUILD_NLIST_COMM_BEGIN

    par::Mpi_Alltoallv_sparse<ot::TreeNode>( sendKpPtr, numKeysSendP, sendOffsetsP,
        recvKpPtr, numKeysRecvP, recvOffsetsP, m_mpiCommActive);

  par::Mpi_Alltoallv_sparse<ot::TreeNode>( sendKsPtr, numKeysSendS, sendOffsetsS,
      recvKsPtr, numKeysRecvS, recvOffsetsS, m_mpiCommActive);

  PROF_BUILD_NLIST_COMM_END

    //Store the sizes, will need it later
    unsigned int actualRecvKpSz = static_cast<unsigned int>(recvKp.size());
  unsigned int actualRecvKsSz = static_cast<unsigned int>(recvKs.size());

  //The result is sorted except for the fact that there might be multiple
  //roots in the middle. So remove all roots.
  std::vector<TreeNode> tmpRecvKp;
  std::vector<TreeNode> tmpRecvKs;

  std::vector<unsigned int> kp2ActualKp;
  std::vector<unsigned int> ks2ActualKs;

  tmpRecvKp.push_back(rootNode);
  kp2ActualKp.push_back(recvKp.size());
  for(unsigned int i = 0; i < recvKp.size(); i++) {
    if(recvKp[i] > rootNode) {      
      tmpRecvKp.push_back(recvKp[i]);
      kp2ActualKp.push_back(i);
    }
  }

  tmpRecvKs.push_back(rootNode);
  ks2ActualKs.push_back(recvKs.size());
  for(unsigned int i = 0; i < recvKs.size(); i++) {
    if(recvKs[i] > rootNode) {
      tmpRecvKs.push_back(recvKs[i]);
      ks2ActualKs.push_back(i);
    }
  }

  recvKp = tmpRecvKp;
  recvKs = tmpRecvKs;

  tmpRecvKp.clear();
  tmpRecvKs.clear();

#ifdef __DEBUG_DA_NLIST__
  MPI_Barrier(m_mpiCommActive);
  assert(seq::test::isUniqueAndSorted(recvKp));
  assert(seq::test::isUniqueAndSorted(recvKs));
  if(!m_iRankActive) {
    std::cout<<std::endl;
    std::cout<<"Processing Results in Second Ring..."<<std::endl;
    std::cout<<std::endl;
  }
  MPI_Barrier(m_mpiCommActive);
#endif

  unsigned int *oldToNewIdx = NULL;
  bool *fromPrimary = NULL;
  bool *fromSecondary = NULL;
  if(!extraAtEnd.empty()) {
    oldToNewIdx = new unsigned int[extraAtEnd.size()];
    fromPrimary = new bool[extraAtEnd.size()];
    fromSecondary = new bool[extraAtEnd.size()];
  }
  std::vector<seq::IndexHolder<ot::TreeNode> > extraIndices (extraAtEnd.size());

  for (unsigned int i = 0; i < extraAtEnd.size(); i++) {
    unsigned int idx;    
    bool found = seq::maxLowerBound<ot::TreeNode >(recvKp, primaryKeys[i], idx, NULL, NULL);

    if (found) {

      if (recvKp[idx] == rootNode) {
        found = false;
        fromPrimary[i] = false;
      } else {
#ifdef __DEBUG_DA_NLIST__
        assert( recvKp[idx].getFlag() & ot::TreeNode::NODE );
#endif
        if( recvKp[idx].getAnchor() != primaryKeys[i].getAnchor() ) {
          //Treat just like it was root.
          found = false;
          fromPrimary[i] = false;
        }else {
          oldToNewIdx[i] = idx; 
          fromPrimary[i] = true;
          fromSecondary[i] = false;
          TreeNode tmp = recvKp[idx];
          tmp.setWeight(extraAtEnd[i].getWeight());
          extraAtEnd[i] = tmp;
          extraIndices[i].value = &(extraAtEnd[i]);
          extraIndices[i].index = i;
        }
      }        
    }else {
      assert(false);
    }

    if (!found) {
      unsigned int idx;    
      found = seq::maxLowerBound<ot::TreeNode >(recvKs, secondaryKeys[i], idx, NULL, NULL);

#ifdef __DEBUG_DA_NLIST__
      assert( found );
      assert( recvKs[idx] != rootNode );
      assert( recvKs[idx].getFlag() & ot::TreeNode::NODE );
      assert( recvKs[idx].getAnchor() == secondaryKeys[i].getAnchor() );
#endif
      /*
         Handle Special Case: When a pre-ghost element did not find its primary key in the local search, the secondary key is never searched for. In such situations, when the primary result is not usable and we must pick the secondary result got from the explicit parallel search, we must first check if this secondary result is already present in the local buffer (either own or elements got through a-priori comm.) If so, we must use the one that is already present and discard the one got through the parallel search. This must be done in order to avoid duplicate octants in the local buffer.
         */
      unsigned int idxTmp;
      ot::TreeNode* inPtr = NULL;
      if(!in.empty()) {
        inPtr = &(*(in.begin()));
      }
      found = seq::BinarySearch<ot::TreeNode>(inPtr,
          static_cast<unsigned int>(in.size()), recvKs[idx], &idxTmp);

      if(found) {
        fromSecondary[i] = false;
        oldToNewIdx[i] = idxTmp;
        TreeNode tmp = in[idxTmp];
        tmp.setWeight(extraAtEnd[i].getWeight());
        extraAtEnd[i] = tmp;
        extraIndices[i].value = &(extraAtEnd[i]);
        extraIndices[i].index = i;
      }else {    
        fromSecondary[i] = true;
        oldToNewIdx[i] = idx; 
        TreeNode tmp = recvKs[idx];
        tmp.setWeight(extraAtEnd[i].getWeight());
        extraAtEnd[i] = tmp;
        extraIndices[i].value = &(extraAtEnd[i]);
        extraIndices[i].index = i;
      }
    }
  }//end for i

  primaryKeys.clear();
  secondaryKeys.clear();
  extraAtEnd.clear();

  std::sort(extraIndices.begin(),extraIndices.end());

  std::vector<ot::TreeNode> sortedUniqueExtras(extraIndices.size());    
  std::vector<unsigned int> chosenIndexIntoOld(extraIndices.size());

  //Old refers to the unsorted and non-unique extraAtEnd list 
  std::vector<std::vector< unsigned int> > indicesIntoOld(extraIndices.size());

  //LUT refers to localOcts (in)
  std::vector<std::vector< unsigned int> > indicesIntoLUT(extraIndices.size());

  if (!extraIndices.empty()) {
    sortedUniqueExtras[0] = *(extraIndices[0].value);
    chosenIndexIntoOld[0] = extraIndices[0].index;

    indicesIntoOld[0].push_back(extraIndices[0].index);
    indicesIntoLUT[0].push_back(extraIndices[0].value->getWeight());
  }
  //Make  Unique and concatenate ranks.
  if (extraIndices.size() >= 2) {
    unsigned int tmpSize=1;
    for (unsigned int i=1;i<extraIndices.size();i++) {
      if (sortedUniqueExtras[tmpSize-1] != *(extraIndices[i].value)) {
        //New entry
        sortedUniqueExtras[tmpSize] = *(extraIndices[i].value);
        chosenIndexIntoOld[tmpSize] =  extraIndices[i].index;

        indicesIntoOld[tmpSize].push_back(extraIndices[i].index);
        indicesIntoLUT[tmpSize].push_back(extraIndices[i].value->getWeight());
        tmpSize++;
      } else {
        indicesIntoOld[tmpSize-1].push_back(extraIndices[i].index);
        indicesIntoLUT[tmpSize-1].push_back(extraIndices[i].value->getWeight());
      }
    }//end for

    sortedUniqueExtras.resize(tmpSize);

#ifdef __DEBUG_DA_NLIST__
    assert(seq::test::isUniqueAndSorted<TreeNode >(sortedUniqueExtras));
#endif

    chosenIndexIntoOld.resize(tmpSize);

    indicesIntoOld.resize(tmpSize);
    indicesIntoLUT.resize(tmpSize);      

    for (unsigned int i=0;i<tmpSize;i++) {
      std::sort(indicesIntoOld[i].begin(),indicesIntoOld[i].end());
      std::sort(indicesIntoLUT[i].begin(),indicesIntoLUT[i].end());        
#ifdef __DEBUG_DA_NLIST__
      //Ensuring that all the 8 vertices of an element map to different indices. 
      assert(seq::test::isUniqueAndSorted<unsigned int>(indicesIntoOld[i]));
      assert(seq::test::isUniqueAndSorted<unsigned int>(indicesIntoLUT[i]));
#endif
    }//end for i      

  }

  extraIndices.clear();

  /*
     std::vector<bool> pickPrimary(sortedUniqueExtras.size());
     std::vector<unsigned int> oldIndex(sortedUniqueExtras.size());

  //The same extra key can come from primary or secondary and from multiple
  //pre-ghosts. So pick one.
  //Some arbit default rule is used to make this decision...
  //Here. The rule is simply if this key was some element's primary key use
  //that. In this case, the smallest element is chosen. Else, the last
  //element is chosen (secondary key).
  //Note: indicesIntoOld[i] is sorted and unique.
  for (unsigned int i=0; i < sortedUniqueExtras.size(); i++) {
  pickPrimary[i] = false;
  for (unsigned int j = 0; j < indicesIntoOld[i].size(); j++) {
  oldIndex[i] = indicesIntoOld[i][j];
  if (fromPrimary[indicesIntoOld[i][j]]) {
  pickPrimary[i] = true;         
  break;       
  }
  }//end for j
  }//end for i
  */

  unsigned int *recvCntsExtras = new unsigned int[m_iNpesActive];        

  std::vector<char> pickingP(actualRecvKpSz);
  std::vector<char> pickingS(actualRecvKsSz);

  std::vector<char> pickedP(sendKp.size());
  std::vector<char> pickedS(sendKs.size());

  for (unsigned int i = 0; i < actualRecvKpSz; i++) {
    pickingP[i] = 0;
  }

  for (unsigned int i = 0; i < actualRecvKsSz; i++) {
    pickingS[i] = 0;
  }

  for (unsigned int i = 0; i < m_iNpesActive; i++) {
    recvCntsExtras[i] = 0;
  }

  // std::cout << m_iRankActive << " Correcting LuTs" << std::endl;

  //Do not do an in-place update of nlist. Instead store the corrections. 
  std::vector< std::vector<unsigned int> > secondRingCorrections;

  for (unsigned int i=0;i<sortedUniqueExtras.size();i++) {

    /*
       if (pickPrimary[i]) {
#ifdef __DEBUG_DA_NLIST__
assert( recvKp[oldToNewIdx[oldIndex[i]]] == sortedUniqueExtras[i] );
assert( recvKp[oldToNewIdx[oldIndex[i]]] != rootNode );
#endif
in.push_back(recvKp[oldToNewIdx[oldIndex[i]]]);        
recvCntsExtras[recvKp[oldToNewIdx[oldIndex[i]]].getWeight()]++;
pickingP[kp2ActualKp[oldToNewIdx[oldIndex[i]]]] = 1;
} else {
#ifdef __DEBUG_DA_NLIST__
assert( recvKs[oldToNewIdx[oldIndex[i]]] == sortedUniqueExtras[i] );
assert( recvKs[oldToNewIdx[oldIndex[i]]] != rootNode );
#endif
in.push_back(recvKs[oldToNewIdx[oldIndex[i]]]);        
recvCntsExtras[recvKs[oldToNewIdx[oldIndex[i]]].getWeight()]++;
pickingS[ks2ActualKs[oldToNewIdx[oldIndex[i]]]] = 1;
}

for (unsigned int j =0;j < indicesIntoLUT[i].size(); j++) {
for (unsigned int k=0; k < 8; k++) {
for (unsigned int l = 0;l < indicesIntoOld[i].size(); l++) {
if (nlist[8*(indicesIntoLUT[i][j])+k] == (m_uiLocalBufferSize + indicesIntoOld[i][l])) {
std::vector<unsigned int> ringCorrectionsTuple(3);
ringCorrectionsTuple[0] = indicesIntoLUT[i][j];
ringCorrectionsTuple[1] = k;
ringCorrectionsTuple[2] = (in.size()-1);
secondRingCorrections.push_back(ringCorrectionsTuple);
break;
} // if 
} // for l                   
}  // for k     
} // for j
*/

  if (fromPrimary[chosenIndexIntoOld[i]]) {
#ifdef __DEBUG_DA_NLIST__
    assert( recvKp[oldToNewIdx[chosenIndexIntoOld[i]]] == sortedUniqueExtras[i] );
    assert( recvKp[oldToNewIdx[chosenIndexIntoOld[i]]] != rootNode );
#endif
    in.push_back(recvKp[oldToNewIdx[chosenIndexIntoOld[i]]]);        
    recvCntsExtras[recvKp[oldToNewIdx[chosenIndexIntoOld[i]]].getWeight()]++;
    pickingP[kp2ActualKp[oldToNewIdx[chosenIndexIntoOld[i]]]] = 1;
  } else if(fromSecondary[chosenIndexIntoOld[i]]) {
#ifdef __DEBUG_DA_NLIST__
    assert( recvKs[oldToNewIdx[chosenIndexIntoOld[i]]] == sortedUniqueExtras[i] );
    assert( recvKs[oldToNewIdx[chosenIndexIntoOld[i]]] != rootNode );
#endif
    in.push_back(recvKs[oldToNewIdx[chosenIndexIntoOld[i]]]);        
    recvCntsExtras[recvKs[oldToNewIdx[chosenIndexIntoOld[i]]].getWeight()]++;
    pickingS[ks2ActualKs[oldToNewIdx[chosenIndexIntoOld[i]]]] = 1;
  }

//If the element is neither in primary nor in secondary,
//then 'in' already has it. So no need to append it into 'in' again.
//Note, that we will never select a copy of sortedUniqueExtras[i] from primary or secondary if the element already exists within in.

for (unsigned int j =0;j < indicesIntoLUT[i].size(); j++) {
  for (unsigned int k=0; k < 8; k++) {
    for (unsigned int l = 0;l < indicesIntoOld[i].size(); l++) {
      if (nlist[8*(indicesIntoLUT[i][j])+k] == (m_uiLocalBufferSize + indicesIntoOld[i][l])) {
        std::vector<unsigned int> ringCorrectionsTuple(3);
        ringCorrectionsTuple[0] = indicesIntoLUT[i][j];
        ringCorrectionsTuple[1] = k;

        if( fromPrimary[chosenIndexIntoOld[i]] || fromSecondary[chosenIndexIntoOld[i]] ) {
          ringCorrectionsTuple[2] = (static_cast<unsigned int>(in.size()) - 1);
        }else {
          ringCorrectionsTuple[2] = oldToNewIdx[chosenIndexIntoOld[i]];
        }

        secondRingCorrections.push_back(ringCorrectionsTuple);
        break;
      } // if 
    } // for l                   
  }  // for k     
} // for j

} // for i

for(unsigned int i = 0; i < secondRingCorrections.size(); i++) {
  //Actual Correction here...
  unsigned int elemId = secondRingCorrections[i][0];
  unsigned int vtxId = secondRingCorrections[i][1];
  unsigned int lutVal = secondRingCorrections[i][2];
  nlist[8*elemId+vtxId] = lutVal;

#ifdef __DEBUG_DA_NLIST__
  assert( lutVal < in.size() );
#endif

  // correct hanging node mask ...
  // check if any of the nodes is hanging ...
  unsigned int _x,_y,_z;

  unsigned int x = in[elemId].getX();
  unsigned int y = in[elemId].getY();
  unsigned int z = in[elemId].getZ();
  unsigned int d   = in[elemId].getLevel();
  unsigned int sz = (1u << (m_uiMaxDepth - d));

#ifdef __DEBUG_DA_NLIST__
  //You might have to find the actual node via final corrections
  //(explicit parallel search) although you sent yourself apriori
  //(Second Ring). 
  if( (elemId >= m_uiElementBegin) && (elemId < m_uiPostGhostBegin) ) {
    if(in[lutVal] < in[m_uiElementBegin]) {
      std::cout<<m_iRankActive<<" Trying to send yourself as a Post Ghost ELEMENT for  elemId = "
        <<elemId<<" vtxId = "<<vtxId<<std::endl;
      assert(false);
    }
    if( (m_uiPostGhostBegin < in.size()) && (in[lutVal] >= in[m_uiPostGhostBegin]) ) {
      TreeNode tmpToSend = in[lutVal];
      tmpToSend.setWeight(elemId);
      checkSecondRing.push_back(tmpToSend);
    }
  }
#endif

  _x = in[lutVal].getX(); 
  _y = in[lutVal].getY(); 
  _z = in[lutVal].getZ(); 

#ifdef __DEBUG_DA_NLIST__
  if ( !(in[lutVal].getFlag() & ot::TreeNode::NODE ) ) {
    //Second ring must find a node only 
    assert(false);
  }
#endif

  switch (vtxId) {
    case 0:
      assert(false);
      break;
    case 1:
      if ( ( (x+sz) != _x ) || ( y != _y ) || ( z != _z ) ) {
        m_ucpLutMasks[2*elemId+1] |= (1 << vtxId);
      }
      break;
    case 2:
      if ( ( x != _x ) || ( (y+sz) != _y ) || ( z != _z ) ) {
        m_ucpLutMasks[2*elemId+1] |= (1 << vtxId);
      }
      break;
    case 3:
      if ( ( (x+sz) != _x ) || ( (y+sz) != _y ) || ( z != _z ) ) {
        m_ucpLutMasks[2*elemId+1] |= (1 << vtxId);
      }
      break;
    case 4:
      if ( ( x != _x ) || ( y != _y ) || ( (z+sz) != _z ) ) {
        m_ucpLutMasks[2*elemId+1] |= (1 << vtxId);
      }
      break;
    case 5:
      if ( ( (x+sz) != _x ) || ( y != _y ) || ( (z+sz) != _z ) ) {
        m_ucpLutMasks[2*elemId+1] |= (1 << vtxId);
      }
      break;
    case 6:
      if ( ( x != _x ) || ( (y+sz) != _y ) || ( (z+sz) != _z ) ) {
        m_ucpLutMasks[2*elemId+1] |= (1 << vtxId);
      }
      break;
    case 7:
      if ( ( (x+sz) != _x ) || ( (y+sz) != _y ) || ( (z+sz) != _z ) ) {
        m_ucpLutMasks[2*elemId+1] |= (1 << vtxId);
      }
      break;
  } // switch (vtxId)

}//end for i

secondRingCorrections.clear();

//pickPrimary.clear();
//oldIndex.clear();

sortedUniqueExtras.clear();
chosenIndexIntoOld.clear();

indicesIntoLUT.clear();
indicesIntoOld.clear();

kp2ActualKp.clear();
ks2ActualKs.clear();

sendKp.clear();
sendKs.clear();
recvKp.clear(); 
recvKs.clear();  

if(oldToNewIdx) {
  delete [] oldToNewIdx;
  oldToNewIdx = NULL;
}

if(fromPrimary) {
  delete [] fromPrimary;    
  fromPrimary = NULL;
}

if(fromSecondary) {
  delete [] fromSecondary;
  fromSecondary = NULL;
}

// Set secondary comm variables ...
for (unsigned int i=0; i<m_iNpesActive; i++) {
  if (recvCntsExtras[i]) {
    ScndRecvProcs.push_back(i);
    ScndRecvCounts.push_back(recvCntsExtras[i]); 
#ifdef __DEBUG_DA_NLIST__
    std::cout<<m_iRankActive<<" Scnd Recv P : "<<i<<" Scnd Recv C : "<< recvCntsExtras[i]<<std::endl;
    assert(i != m_iRankActive);
#endif
  }
}

if(recvCntsExtras) {
  delete [] recvCntsExtras;    
  recvCntsExtras = NULL;
}

#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
if(!m_iRankActive) {
  std::cout<<std::endl;
  std::cout<<"Returning Selection in Second Ring..."<<std::endl;
  std::cout<<std::endl;
}
MPI_Barrier(m_mpiCommActive);
#endif

char* pickingPptr = NULL;
char* pickedPptr = NULL;
char* pickingSptr = NULL;
char* pickedSptr = NULL;

if(!pickingP.empty()) {
  pickingPptr = &(*(pickingP.begin()));
}
if(!pickingS.empty()) {
  pickingSptr = &(*(pickingS.begin()));
}
if(!pickedP.empty()) {
  pickedPptr = &(*(pickedP.begin()));
}
if(!pickedS.empty()) {
  pickedSptr = &(*(pickedS.begin()));
}

PROF_BUILD_NLIST_COMM_BEGIN

par::Mpi_Alltoallv_sparse<char>(pickingPptr ,numKeysRecvP, recvOffsetsP , 
    pickedPptr, numKeysSendP, sendOffsetsP, m_mpiCommActive);

par::Mpi_Alltoallv_sparse<char>(pickingSptr ,numKeysRecvS, recvOffsetsS, 
    pickedSptr, numKeysSendS, sendOffsetsS , m_mpiCommActive);

PROF_BUILD_NLIST_COMM_END

pickingP.clear();
pickingS.clear();

#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
if(!m_iRankActive) {
  std::cout<<std::endl;
  std::cout<<"Second ScatterMap..."<<std::endl;
  std::cout<<std::endl;
}
MPI_Barrier(m_mpiCommActive);
#endif

//Create ScatterMap here using sendKp, recvKp, pickedP...
unsigned int cnt=0, pickedCnt=0;
unsigned int pCnt=0, sCnt=0;
for (unsigned int i=0; i < static_cast<unsigned int>(m_iNpesActive); i++) {
  cnt=0;
  pickedCnt=0;
  // merge and pick entries from both primary and secondary lists being sent to procs ...
  while (cnt < (numKeysSendP[i] + numKeysSendS[i]) ) {
    // sOver := sCnt >= (sendOffsetsS[i] + numKeysSendS[i]);
    // pOver := pCnt >= (sendOffsetsP[i] + numKeysSendP[i]);
    // sRemains = sCnt < (sendOffsetsS[i] + numKeysSendS[i]);
    // pRemains = pCnt < (sendOffsetsP[i] + numKeysSendP[i]);
    if ( (sCnt >= (sendOffsetsS[i] + numKeysSendS[i])) ||
        ( ( pCnt < (sendOffsetsP[i] + numKeysSendP[i]) ) &&
          ( idxSendKp[pCnt] <= idxSendKs[sCnt]) ) ) {
#ifdef __DEBUG_DA_NLIST__
      if(idxSendKp[pCnt] == -1) {
        assert(!pickedP[pCnt]);
      }
#endif
      if (pickedP[pCnt]) {
        ScndScatterMap.push_back(idxSendKp[pCnt]);
        pickedCnt++;
      }
      pCnt++;
    } else if ( ( pCnt >=
          static_cast<unsigned int>(sendOffsetsP[i] + numKeysSendP[i]) ) ||
        ( sCnt < static_cast<unsigned int>(sendOffsetsS[i] + numKeysSendS[i]) ) ) {
#ifdef __DEBUG_DA_NLIST__
      if(idxSendKs[sCnt] == -1) {
        assert(!pickedS[sCnt]);
      }
#endif
      if (pickedS[sCnt]) {
        ScndScatterMap.push_back(idxSendKs[sCnt]);
        pickedCnt++;
      }
      sCnt++;
    } else {
      std::cout << "Scnd skipped both" << std::endl;
    }
    cnt++;
  }//end while
  if (pickedCnt) {
    ScndSendProcs.push_back(i);  
    ScndSendCounts.push_back(pickedCnt);
#ifdef __DEBUG_DA_NLIST__
    std::cout<<m_iRankActive<<" Scnd Send P : "<<i<<" Scnd Send C : "<< pickedCnt<<std::endl;
    assert(i != m_iRankActive);
#endif
  }
}//end for i
#ifdef __DEBUG_DA_NLIST__
for(unsigned int i = 0; i < m_iNpesActive; i++) {
  unsigned int numToSend_i = 0;
  for(unsigned int j = sendOffsetsP[i]; j < (sendOffsetsP[i] + numKeysSendP[i]); j++) {
    if(pickedP[j]) {
      numToSend_i++;
    }
  }//end for j
  for(unsigned int j = sendOffsetsS[i]; j < (sendOffsetsS[i] + numKeysSendS[i]); j++) {
    if(pickedS[j]) {
      numToSend_i++;
    }
  }//end for j
  if(numToSend_i) {
    std::cout<<m_iRankActive<<" Actually should send "<<numToSend_i<<" to "<<i<<std::endl;
  }
}//end for i
#endif

//Reset LocalBuffer Size
m_uiLocalBufferSize = static_cast<unsigned int>(in.size());

idxSendKp.clear();
idxSendKs.clear();

delete [] sendOffsetsP;
sendOffsetsP = NULL;

delete [] sendOffsetsS;
sendOffsetsS = NULL;

delete [] numKeysSendP;
numKeysSendP = NULL;

delete [] numKeysSendS;
numKeysSendS = NULL;

delete [] numKeysRecvP;
numKeysRecvP = NULL;

delete [] recvOffsetsP;    
recvOffsetsP = NULL;

delete [] recvOffsetsS;
recvOffsetsS = NULL;

delete [] numKeysRecvS;
numKeysRecvS = NULL;

pickedP.clear();
pickedS.clear();   

#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
if(!m_iRankActive) {
  std::cout<<std::endl;
  std::cout<<"Finished Correction for Second Ring..."<<std::endl;
  std::cout<<std::endl;
}
MPI_Barrier(m_mpiCommActive);
#endif

//~~~~~~~~~~~~~~~~Correction for Second Ring Ends~~~~~~~~~~~~~~//

#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
//CHECK if apriori comm for Second ring actually works.
//Check if primary scattermap contains the right element and if it is
//sent to the right processor
//checkSecondRing[i] is a post-ghost element on the local proc, to which one of
//my own elements are pointing to. The corresponding index of my own element is
//stored in the weight of checkSecondRing[i]. 
std::vector<std::vector<unsigned int> > secondRingExpectedScatterMap(m_iNpesActive);
for(unsigned int i = 0; i < checkSecondRing.size(); i++) {
  unsigned int bId;
  //Find the processor that owns the node.
  bool bucketFound = seq::maxLowerBound<ot::TreeNode>(m_tnMinAllBlocks, checkSecondRing[i], bId, NULL, NULL);
  assert(bucketFound);
  assert(bId < m_iNpesActive);
  assert( m_tnMinAllBlocks[bId] <= checkSecondRing[i] );
  secondRingExpectedScatterMap[bId].push_back(checkSecondRing[i].getWeight());
}

assert(seq::test::isUniqueAndSorted(m_uipSendProcs));

unsigned int toSendCnt = 0;
unsigned int currOffset = 0;
for(unsigned int i = 0; i < m_iNpesActive; i++) {
  seq::makeVectorUnique<unsigned int>(secondRingExpectedScatterMap[i], false) ;
  if( (toSendCnt < m_uipSendProcs.size()) && (m_uipSendProcs[toSendCnt] == i) ) {
    for(unsigned int j = currOffset; j < (currOffset + m_uipSendCounts[toSendCnt] - 1); j++ ) {
      assert(m_uipScatterMap[j] < m_uipScatterMap[j+1]);
    }
    unsigned int st = 0;    
    for(unsigned int j = 0; j < secondRingExpectedScatterMap[i].size(); j++) {
      bool foundIt = false;
      for(unsigned int k = st; k < m_uipSendCounts[toSendCnt]; k++) {
        assert( (currOffset + k) < m_uipScatterMap.size() );
        if( m_uipScatterMap[currOffset + k] == secondRingExpectedScatterMap[i][j] ) {
          foundIt = true;
          //secondRingExpectedScatterMap[i] is Sorted and Unique.
          //Primary ScatterMap is also sorted and Unique for each processor 
          st = (k+1); 
          break;
        }
      }//end for k
      if(!foundIt) {
        std::cout<<m_iRankActive<<" should have sent in["<<
          secondRingExpectedScatterMap[i][j]<<"]: "<<in[secondRingExpectedScatterMap[i][j]]
          <<"  to "<<i<<" in primary. But did not."<<std::endl; 
        assert(false);
      }
    }//end for j
    currOffset += m_uipSendCounts[toSendCnt];
    toSendCnt++;
  }else if( !(secondRingExpectedScatterMap[i].empty()) ) {
    std::cout<<m_iRankActive<<" should have sent some elements to "<<i<<" in primary."
      <<" But did not. toSendCnt: "<<toSendCnt<<" primarySz: "<<(m_uipSendProcs.size())<<std::endl; 
    std::cout<<m_iRankActive<<" For example, This element was not found in primary ScatterMap: "
      <<in[secondRingExpectedScatterMap[i][0]]<<std::endl;
    assert(false);
  }//end if ranks match
}//end for i

secondRingExpectedScatterMap.clear();
checkSecondRing.clear();

if(!m_iRankActive) {
  std::cout<<" Finished Checking Apriori send for Second Ring.."<<std::endl;
}
MPI_Barrier(m_mpiCommActive);
#endif

//~~~~~~~~~~~Sort In, Flag Corrections as FOREIGN, update Luts~~~~~~~~~~~~~~~~~

// Now resort in, and correct nlist and scattermap.
std::vector<seq::IndexHolder<ot::TreeNode> > inHolder(in.size());
for (unsigned int i=0; i < in.size(); i++) {
  in[i].setWeight(i);
  inHolder[i].value = &in[i];
  if( (i >= nelem) || ( (i < iLoopEnd) &&
        (m_ucpLutMasks[(2*i) + 1] == ot::DA_FLAGS::FOREIGN) ) ) {
    //Mark old Foreign, new octants, postghosts and boundaries
    inHolder[i].index = 0;
  } else {
    inHolder[i].index = 1;
  }
} // end i for resort set weights ...

#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
if(!m_iRankActive) {
  std::cout<<std::endl;
  std::cout<<"created inHolder..."<<std::endl;
  std::cout<<std::endl;
}
MPI_Barrier(m_mpiCommActive);
#endif

// Now resort ...
std::sort( inHolder.begin(), inHolder.end() );

assert(m_uiElementBegin < in.size());
TreeNode  myFirstOctant = in[m_uiElementBegin];
assert( (m_uiPostGhostBegin - 1) < in.size() );
TreeNode  myLastOctant = in[m_uiPostGhostBegin - 1];

//Correct nelem...
nelem = 0;
if(m_uiElementSize) {
  while( (nelem < in.size()) &&
      ((*(inHolder[nelem].value)) <= in[m_uiElementEnd - 1]) ) {
    nelem++; 
  }
} else {
  while( ( nelem < in.size() ) &&
      ( (*(inHolder[nelem].value)) < myFirstOctant ) &&
      ( !((inHolder[nelem].value)->getFlag() & ot::TreeNode::BOUNDARY) ) ) {
    nelem++; 
  }
}

#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
if(!m_iRankActive) {
  std::cout<<std::endl;
  std::cout<<"Corrected NELEM..."<<std::endl;
  std::cout<<std::endl;
}
MPI_Barrier(m_mpiCommActive);
#endif

// Wts of in gives us new->old mapping. However to correct
// nlist we need old->new mapping. So we need to generate this.
std::vector<unsigned int> oldToNew(in.size());

for (unsigned int i = 0; i < in.size(); i++) {
  oldToNew[inHolder[i].value->getWeight()] = i;
} // end i for remap indices

#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
if(!m_iRankActive) {
  std::cout<<std::endl;
  std::cout<<"Computed Old2New..."<<std::endl;
  std::cout<<std::endl;
}
MPI_Barrier(m_mpiCommActive);
#endif

//Re-shuffle in.
std::vector<ot::TreeNode> tmpIn(in.size());
for(unsigned int i = 0;  i < in.size(); i++) {
  tmpIn[i] = *(inHolder[i].value);
  tmpIn[i].setWeight(inHolder[i].index);
}//end for i

in = tmpIn;
tmpIn.clear();
inHolder.clear();

#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
if(!m_iRankActive) {
  std::cout<<std::endl;
  std::cout<<"Reshuffled In."<<std::endl;
  std::cout<<std::endl;
}
MPI_Barrier(m_mpiCommActive);
#endif

//Reset Counters...
// PreGhosts ....
m_uiPreGhostElementSize = 0;
while( (m_uiPreGhostElementSize < in.size()) && (in[m_uiPreGhostElementSize] < myFirstOctant) ) {
  if ( !(in[m_uiPreGhostElementSize].getFlag() & ot::TreeNode::BOUNDARY) ) {
    m_uiPreGhostElementSize++;
  }else {
    break;
  }
}//end while

if(nelem != (m_uiPreGhostElementSize + m_uiElementSize) ) {
  std::cout<<"Processor "<<m_iRankAll<<" failing: nelem = "<<nelem
    <<" pgElemSz = "<<m_uiPreGhostElementSize<<" elemSz = "<<m_uiElementSize<<std::endl;
}
assert( nelem == (m_uiPreGhostElementSize + m_uiElementSize) );
#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
if(!m_iRankActive) {
  std::cout<<std::endl;
  std::cout<<"PreGhost Counters Reset."<<std::endl;
  std::cout<<std::endl;
}
MPI_Barrier(m_mpiCommActive);
#endif

//Mine ...
m_uiElementBegin = m_uiPreGhostElementSize;
m_uiPostGhostBegin = m_uiPreGhostElementSize;
while( (m_uiPostGhostBegin < in.size()) && (in[m_uiPostGhostBegin] <= myLastOctant) ) {
  if(in[m_uiPostGhostBegin] < myFirstOctant) {
    m_uiElementBegin++;
  }
  m_uiPostGhostBegin++;
}//end while
m_uiElementEnd = m_uiElementBegin + m_uiElementSize;

#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
if(!m_iRankActive) {
  std::cout<<std::endl;
  std::cout<<"All Counters Reset."<<std::endl;
  std::cout<<std::endl;
}
MPI_Barrier(m_mpiCommActive);
#endif

std::vector<unsigned int> nlistNew(8*nelem);
std::vector<unsigned char> hnMaskNew(nelem);

#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
if(!m_iRankActive) {
  std::cout<<std::endl;
  std::cout<<"Allocated memory for the new nlist and masks."<<std::endl;
  std::cout<<std::endl;
}
MPI_Barrier(m_mpiCommActive);
#endif

// correct nLists ...
// The loop is in the order of the old Indices ...
//The first time only the nlist of preghosts are corrected.
//Since, the nlist of own elements are not built in the first
//outer-loop (numFullLoopCtr).
//The second time the preghosts will be visited once again.
//Skip correction of nlist for Foreigns 
for (unsigned int i = 0; i < iLoopEnd; i++) {
#ifdef __DEBUG_DA_NLIST__
  assert(i < oldToNew.size() );
  assert(oldToNew[i] < nelem);
  if( numFullLoopCtr == 0 ) {
    assert(oldToNew[i] < m_uiPreGhostElementSize);
  }
  assert( (2*i+1) < m_ucpLutMasks.size() );
#endif

  hnMaskNew[oldToNew[i]] = m_ucpLutMasks[2*i + 1];

  if(m_ucpLutMasks[2*i+1] == ot::DA_FLAGS::FOREIGN) {
    continue;
  }

  unsigned int iiOld = 8*i;
  unsigned int iiNew = 8*oldToNew[i];
  for (unsigned int j=0; j<8; j++) {
    // remap ...
#ifdef __DEBUG_DA_NLIST__
    assert( (iiNew + j) < nlistNew.size() );
    assert( (iiOld + j) < nlist.size()  );
    if( nlist[iiOld + j] >= oldToNew.size() ) {
      std::cout<<m_iRankActive<<": oldToNew.size(): "<<oldToNew.size()<<" nlist[8*"<<i<<"+"<<j<<"]: "<<nlist[iiOld+j]<<std::endl 
                                                       <<"oldToNew["<<i<<"] = "<<oldToNew[i]<<" elem: "<<in[oldToNew[i]]<<std::endl;                                        
      assert(false);
    }
#endif

    nlistNew[iiNew + j] = oldToNew[ nlist[iiOld + j] ];
  }//end for j
}//end for i

#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
if(!m_iRankActive) {
  std::cout<<std::endl;
  std::cout<<"Corrected Nlist."<<std::endl;
  std::cout<<std::endl;
}
MPI_Barrier(m_mpiCommActive);
#endif

//Replace the old lists by the new lists and clear the new
//lists...

m_ucpLutMasks.resize(2*nelem);
if(numFullLoopCtr) {
  for(unsigned int i = 0; i < nelem; i++) {
    if(in[i].getWeight() == 0) {
      m_ucpLutMasks[2*i + 1] = ot::DA_FLAGS::FOREIGN;
    }else {
      m_ucpLutMasks[2*i + 1] = hnMaskNew[i];
    }
  }//end for i
}else {
  for(unsigned int i = 0; i < m_uiPreGhostElementSize; i++) {
    if(in[i].getWeight() == 0) {
      m_ucpLutMasks[2*i + 1] = ot::DA_FLAGS::FOREIGN;
    }else {
      m_ucpLutMasks[2*i + 1] = hnMaskNew[i];
    }
  }//end for i
}
hnMaskNew.clear();

for(unsigned int i = 0; i < in.size(); i++) {
  in[i].setWeight(1);
}

nlist = nlistNew;
nlistNew.clear();

#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
if(!m_iRankActive) {
  std::cout<<std::endl;
  std::cout<<"Corrected HnMasks and Nlists."<<std::endl;
  std::cout<<std::endl;
}
MPI_Barrier(m_mpiCommActive);
#endif

// Merge scattermap and scnScatterMap
std::vector<unsigned int>               tmpScatterMap;
std::vector<unsigned int>               tmpSendProcs;
std::vector<unsigned int>               tmpSendCounts;
std::vector<unsigned int>               tmpRecvProcs;
std::vector<unsigned int>               tmpRecvCounts;

//Assumes primary and secondary lists to be sorted independently. the
//result will be sorted.
//There is a 1-1 mapping between counts and procs.
unsigned int primaryCnt = 0;
unsigned int secondaryCnt = 0;
while( (primaryCnt < m_uipSendProcs.size()) && (secondaryCnt < ScndSendProcs.size()) ) {
  if( m_uipSendProcs[primaryCnt] < ScndSendProcs[secondaryCnt] ) {
    tmpSendProcs.push_back(m_uipSendProcs[primaryCnt]);
    tmpSendCounts.push_back(m_uipSendCounts[primaryCnt]);
    primaryCnt++;
  }else if( m_uipSendProcs[primaryCnt] > ScndSendProcs[secondaryCnt] ) {
    tmpSendProcs.push_back(ScndSendProcs[secondaryCnt]);
    tmpSendCounts.push_back(ScndSendCounts[secondaryCnt]);
    secondaryCnt++;
  }else {
    //if both are equal select p only from one. Arbitrarily, the default
    //is picked as primary
    tmpSendProcs.push_back(m_uipSendProcs[primaryCnt]);
    //Sum the counts from the primary and secondary...
    tmpSendCounts.push_back(m_uipSendCounts[primaryCnt] + ScndSendCounts[secondaryCnt]);
    primaryCnt++;
    //skip secondary
    secondaryCnt++;
  }
}//end while

//only primary remains
while( primaryCnt < m_uipSendProcs.size() ) {
  tmpSendProcs.push_back(m_uipSendProcs[primaryCnt]);
  tmpSendCounts.push_back(m_uipSendCounts[primaryCnt]);
  primaryCnt++;
}

//only secondary remains
while( secondaryCnt < ScndSendProcs.size() ) {
  tmpSendProcs.push_back(ScndSendProcs[secondaryCnt]);
  tmpSendCounts.push_back(ScndSendCounts[secondaryCnt]);
  secondaryCnt++;
}

#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
assert( seq::test::isUniqueAndSorted<unsigned int>(m_uipSendProcs) );
assert( seq::test::isUniqueAndSorted<unsigned int>(ScndSendProcs) );
assert( seq::test::isUniqueAndSorted<unsigned int>(tmpSendProcs) );
for( unsigned int i = 0; i < m_uipSendProcs.size(); i++) {
  std::cout<<m_iRankActive<<" --> "<<m_uipSendProcs[i]<<" (P) "<<m_uipSendCounts[i]<<std::endl;
}
std::cout<<m_iRankActive<<" P-Scatter: "<<m_uipScatterMap.size()<<std::endl;
assert( m_uipSendProcs.size() == m_uipSendCounts.size() );
MPI_Barrier(m_mpiCommActive);
for( unsigned int i = 0; i < ScndSendProcs.size(); i++) {
  std::cout<<m_iRankActive<<" --> "<<ScndSendProcs[i]<<" (S) "<<ScndSendCounts[i]<<std::endl;
}
std::cout<<m_iRankActive<<" S-Scatter: "<<ScndScatterMap.size()<<std::endl;
assert( ScndSendProcs.size() == ScndSendCounts.size() );
MPI_Barrier(m_mpiCommActive);
for( unsigned int i = 0; i < tmpSendProcs.size(); i++) {
  std::cout<<m_iRankActive<<" --> "<<tmpSendProcs[i]<<" (T) "<<tmpSendCounts[i]<<std::endl;
}
assert( tmpSendProcs.size() == tmpSendCounts.size() );
MPI_Barrier(m_mpiCommActive);
if(!m_iRankActive) {
  std::cout<<std::endl;
  std::cout<<"Corrected SendCnts and SendProcs."<<std::endl;
  std::cout<<std::endl;
}
MPI_Barrier(m_mpiCommActive);
#endif

//Assumes primary and secondary Scattermaps to be sorted
//independently for each processors portion. the
//result will also be sorted in chunks.
primaryCnt = 0;
secondaryCnt = 0;
unsigned int primarySz = 0;
unsigned int secondarySz = 0;
while( (primaryCnt < m_uipSendProcs.size()) && (secondaryCnt < ScndSendProcs.size()) ) {
  if( m_uipSendProcs[primaryCnt] < ScndSendProcs[secondaryCnt] ) {
    unsigned int numSent = 0;
    while( numSent < m_uipSendCounts[primaryCnt] ) {
      tmpScatterMap.push_back(oldToNew[m_uipScatterMap[primarySz++]]);
      numSent++;
    }
    primaryCnt++;
  }else if( m_uipSendProcs[primaryCnt] > ScndSendProcs[secondaryCnt] ) {
    unsigned int numSent = 0;
    while( numSent < ScndSendCounts[secondaryCnt] ) {
      tmpScatterMap.push_back(oldToNew[ScndScatterMap[secondarySz++]]);
      numSent++;
    }
    secondaryCnt++;
  }else {
    //Both primary and secondary are sending to the same processor so
    //merge.
    unsigned int nP = 0;
    unsigned int nS = 0;

    //Both are not over.
    while( ( nP < m_uipSendCounts[primaryCnt] ) && ( nS < ScndSendCounts[secondaryCnt] ) ) {
#ifdef __DEBUG_DA_NLIST__
      //The same octant can not be sent to the same processor in both
      //the primary and secondary lists.
      if( m_uipScatterMap[primarySz] == ScndScatterMap[secondarySz] ) {
        std::cout<<m_iRankActive<<" is sending "<<
          in[oldToNew[m_uipScatterMap[primarySz]]]<<" to "
          <<m_uipSendProcs[primaryCnt]<<" in both primary and secondary."<<std::endl;
        assert(false);
      }
#endif
      //Both primary and secondary can only send own elements, the
      //relative ordering of the own elements in the old and new
      //ordering should be the same.
      if( m_uipScatterMap[primarySz] < ScndScatterMap[secondarySz] ) {
#ifdef __DEBUG_DA_NLIST__
        assert( oldToNew[m_uipScatterMap[primarySz]] < oldToNew[ScndScatterMap[secondarySz]] );
#endif
        tmpScatterMap.push_back(oldToNew[m_uipScatterMap[primarySz++]]);
        nP++;
      }else {
#ifdef __DEBUG_DA_NLIST__
        assert( oldToNew[m_uipScatterMap[primarySz]] > oldToNew[ScndScatterMap[secondarySz]] );
#endif
        tmpScatterMap.push_back(oldToNew[ScndScatterMap[secondarySz++]]);
        nS++;
      }
    }

    //Only primary remains
    while( nP < m_uipSendCounts[primaryCnt] ) {
      tmpScatterMap.push_back(oldToNew[m_uipScatterMap[primarySz++]]);
      nP++;
    }

    //only secondary remains
    while( nS < ScndSendCounts[secondaryCnt] ) {
      tmpScatterMap.push_back(oldToNew[ScndScatterMap[secondarySz++]]);
      nS++;
    }

    primaryCnt++;
    secondaryCnt++;
  }//end if-else
}//end while

#ifdef __DEBUG_DA_NLIST__
//ScatterMaps and Procs must finish together.
if( primarySz < m_uipScatterMap.size() ) {
  assert(primaryCnt < m_uipSendProcs.size());
}else {
  assert(primaryCnt == m_uipSendProcs.size());
}

if( secondarySz < ScndScatterMap.size() ) {
  assert(secondaryCnt < ScndSendProcs.size());
}else {
  assert(secondaryCnt == ScndSendProcs.size());
}
#endif

//Only primary remains.
while( primarySz < m_uipScatterMap.size() ) {
  tmpScatterMap.push_back(oldToNew[m_uipScatterMap[primarySz++]]);
}

//Only secondary remains.
while( secondarySz < ScndScatterMap.size() ) {
  tmpScatterMap.push_back(oldToNew[ScndScatterMap[secondarySz++]]);
}

#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
if(!m_iRankActive) {
  std::cout<<std::endl;
  std::cout<<"Merged Primary and Secondary Scatter Maps."<<std::endl;
  std::cout<<std::endl;
}

std::cout<<m_iRankActive<<" pScatterSz: "<<m_uipScatterMap.size()<<" sScatterSz: "
<<ScndScatterMap.size()<<" tScatterSz: "<<tmpScatterMap.size()<<std::endl;

unsigned int debug_pCnt = 0;
unsigned int debug_sCnt = 0;
unsigned int debug_tCnt = 0;
unsigned int debug_pOff = 0;
unsigned int debug_sOff = 0;
unsigned int debug_tOff = 0;
for(unsigned int i = 0; i< m_iNpesActive; i++) {
  std::vector<unsigned int>  pSctList;
  std::vector<unsigned int>  sSctList;
  std::vector<unsigned int>  tSctList;
  if( (debug_pCnt < m_uipSendProcs.size()) && (m_uipSendProcs[debug_pCnt] == i) ) {
    for(unsigned int j = 0; j < m_uipSendCounts[debug_pCnt]; j++) {
      pSctList.push_back(m_uipScatterMap[debug_pOff + j]);
    }
    debug_pOff += m_uipSendCounts[debug_pCnt];
    debug_pCnt++;
  }
  if( (debug_sCnt < ScndSendProcs.size()) && (ScndSendProcs[debug_sCnt] == i) ) {
    for(unsigned int j = 0; j < ScndSendCounts[debug_sCnt]; j++) {
      sSctList.push_back(ScndScatterMap[debug_sOff + j]);
    }
    debug_sOff += ScndSendCounts[debug_sCnt];
    debug_sCnt++;
  }
  if( (debug_tCnt < tmpSendProcs.size()) && (tmpSendProcs[debug_tCnt] == i) ) {
    for(unsigned int j = 0; j < tmpSendCounts[debug_tCnt]; j++) {
      tSctList.push_back(tmpScatterMap[debug_tOff + j]);
    }
    debug_tOff += tmpSendCounts[debug_tCnt];
    debug_tCnt++;
  }

  assert( tSctList.size() == (pSctList.size() + sSctList.size()) );

  assert( seq::test::isUniqueAndSorted<unsigned int>(pSctList) );
  assert( seq::test::isUniqueAndSorted<unsigned int>(sSctList) );

  //Explicitly merge p and s into tt.
  std::vector<unsigned int> ttList;
  for(unsigned int j=0; j < sSctList.size(); j++) {
    ttList.push_back(oldToNew[sSctList[j]]);
  }
  for(unsigned int j=0; j < pSctList.size(); j++) {
    ttList.push_back(oldToNew[pSctList[j]]);
  }

  sort(ttList.begin(),ttList.end());
  assert( seq::test::isUniqueAndSorted<unsigned int>(ttList) );

  // was assert(ttList == tSctList);
  assert(ttList.size() == tSctList.size());
  for (unsigned int j = 0; j < ttList.size(); j++) {
    if (ttList[j] != tSctList[j]) {
      std::cout << m_iRankActive << ": MergeFailed at " << j << std::endl;
      std::cout << ttList[j] << " != "  << tSctList[j] << std::endl;
      for (unsigned int k = 0; k < pSctList.size(); k++) {
        std::cout << ttList[k] << ", " << tSctList[k]
          << ", " << pSctList[k] << std::endl;
      }
      for (unsigned int k = 0; k < sSctList.size(); k++) {
        std::cout << sSctList[k] << std::endl;
      }
      assert(false);
    }
  }//end for j

  ttList.clear();
  pSctList.clear();
  sSctList.clear();
  tSctList.clear();
}//end for i

assert( debug_pCnt == m_uipSendProcs.size() );
assert( debug_sCnt == ScndSendProcs.size() );
assert( debug_tCnt == tmpSendProcs.size() );

assert( tmpScatterMap.size() == (m_uipScatterMap.size() + ScndScatterMap.size()) );
MPI_Barrier(m_mpiCommActive);
#endif

oldToNew.clear();

primaryCnt = 0;
secondaryCnt = 0;
while( (primaryCnt < m_uipRecvProcs.size()) && (secondaryCnt < ScndRecvProcs.size()) ) {
  if( m_uipRecvProcs[primaryCnt] < ScndRecvProcs[secondaryCnt] ) {
    tmpRecvProcs.push_back(m_uipRecvProcs[primaryCnt]);
    tmpRecvCounts.push_back(m_uipRecvCounts[primaryCnt]);
    primaryCnt++;
  }else if( m_uipRecvProcs[primaryCnt] > ScndRecvProcs[secondaryCnt] ) {
    tmpRecvProcs.push_back(ScndRecvProcs[secondaryCnt]);
    tmpRecvCounts.push_back(ScndRecvCounts[secondaryCnt]);
    secondaryCnt++;
  }else {
    //if both are equal select only from primary
    tmpRecvProcs.push_back(m_uipRecvProcs[primaryCnt]);
    tmpRecvCounts.push_back(m_uipRecvCounts[primaryCnt] + ScndRecvCounts[secondaryCnt]);
    primaryCnt++;
    //skip secondary
    secondaryCnt++;
  }
}

while( primaryCnt < m_uipRecvProcs.size() ) {
  tmpRecvProcs.push_back(m_uipRecvProcs[primaryCnt]);
  tmpRecvCounts.push_back(m_uipRecvCounts[primaryCnt]);
  primaryCnt++;
}

while( secondaryCnt < ScndRecvProcs.size() ) {
  tmpRecvProcs.push_back(ScndRecvProcs[secondaryCnt]);
  tmpRecvCounts.push_back(ScndRecvCounts[secondaryCnt]);
  secondaryCnt++;
}

#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
assert( seq::test::isUniqueAndSorted<unsigned int>(m_uipRecvProcs) );
assert( seq::test::isUniqueAndSorted<unsigned int>(ScndRecvProcs) );
assert( seq::test::isUniqueAndSorted<unsigned int>(tmpRecvProcs) );
for( unsigned int i = 0; i < m_uipRecvProcs.size(); i++) {
  std::cout<<m_iRankActive<<" <-- "<<m_uipRecvProcs[i]<<" (P) "<<m_uipRecvCounts[i]<<std::endl;
}
assert( m_uipRecvProcs.size() == m_uipRecvCounts.size() );
MPI_Barrier(m_mpiCommActive);
for( unsigned int i = 0; i < ScndRecvProcs.size(); i++) {
  std::cout<<m_iRankActive<<" <-- "<<ScndRecvProcs[i]<<" (S) "<<ScndRecvCounts[i]<<std::endl;
}
assert( ScndRecvProcs.size() == ScndRecvCounts.size() );
MPI_Barrier(m_mpiCommActive);
for( unsigned int i = 0; i < tmpRecvProcs.size(); i++) {
  std::cout<<m_iRankActive<<" <-- "<<tmpRecvProcs[i]<<" (T) "<<tmpRecvCounts[i]<<std::endl;
}
assert( tmpRecvProcs.size() == tmpRecvCounts.size() );
MPI_Barrier(m_mpiCommActive);
#endif

ScndScatterMap.clear();
ScndSendProcs.clear();
ScndSendCounts.clear();
ScndRecvProcs.clear();
ScndRecvCounts.clear();

m_uipScatterMap  =  tmpScatterMap;
m_uipSendProcs   =  tmpSendProcs;
m_uipSendCounts  =  tmpSendCounts;
m_uipRecvProcs   =  tmpRecvProcs;
m_uipRecvCounts  =  tmpRecvCounts;

tmpScatterMap.clear();
tmpSendProcs.clear();
tmpSendCounts.clear();
tmpRecvProcs.clear();
tmpRecvCounts.clear();

#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
if(!m_iRankActive) {
  std::cout<<std::endl;
  std::cout<<"Corrected RecvCounts and RecvProcs."<<std::endl;
  std::cout<<std::endl;
}
MPI_Barrier(m_mpiCommActive);
#endif

}//end for numFullLoopCtr

// Compute and store the new offsets ...
m_uipSendOffsets.resize(m_uipSendCounts.size());
m_uipRecvOffsets.resize(m_uipRecvCounts.size());

if ( m_uipSendCounts.size() ) {
  m_uipSendOffsets[0] = 0;
  for (unsigned int i=1; i < m_uipSendCounts.size(); i++) {
    m_uipSendOffsets[i] = (m_uipSendCounts[i-1] + m_uipSendOffsets[i-1]);
  }
}

if ( m_uipRecvCounts.size() ) {
  bool adjustedAlready = false;
  if(m_uipRecvProcs[0] < static_cast<unsigned int>(m_iRankActive)) {
    m_uipRecvOffsets[0] = 0;
  }else {
    m_uipRecvOffsets[0] = m_uiPostGhostBegin;
    adjustedAlready = true;
  }
  for (unsigned int i=1; i < m_uipRecvCounts.size(); i++) {
    if( (m_uipRecvProcs[i] < m_iRankActive) || adjustedAlready ) {
      m_uipRecvOffsets[i] = (m_uipRecvCounts[i-1] + m_uipRecvOffsets[i-1]);
    }else {
      m_uipRecvOffsets[i] = m_uiPostGhostBegin;
      adjustedAlready = true;
    }
  }//end for i
}

//Store a copy of the original scattermap first. This will be required for
//communicating ghost elements. Since there are no post-ghost element, we only
//need to send to processors with ranks greater than my rank.
m_uipElemScatterMap.clear();
m_uipElemSendOffsets.clear();
m_uipElemSendProcs.clear();
m_uipElemSendCounts.clear();
for(unsigned int i = 0; i < m_uipSendProcs.size(); i++) {
  if(m_uipSendProcs[i] > static_cast<unsigned int>(m_iRankActive)) {
    unsigned int numElemProcs = static_cast<unsigned int>(m_uipElemSendOffsets.size());
    m_uipElemSendOffsets.push_back(m_uipElemScatterMap.size());
    for(unsigned int j = m_uipSendOffsets[i]; j < (m_uipSendCounts[i] + m_uipSendOffsets[i]); j++) {
      if(m_uipScatterMap[j] < m_uiElementEnd) {
        m_uipElemScatterMap.push_back(m_uipScatterMap[j]);
      }else {
        break;
      }
    }//end for j    
    unsigned int currCount = 
      (static_cast<unsigned int>(m_uipElemScatterMap.size())
       - m_uipElemSendOffsets[numElemProcs]);
    if(currCount) {
      m_uipElemSendProcs.push_back(m_uipSendProcs[i]);
      m_uipElemSendCounts.push_back(currCount);
    }else {
      m_uipElemSendOffsets.resize(numElemProcs);
    }
  }
}//end for i

m_uipElemRecvOffsets.clear();
m_uipElemRecvProcs.clear();
m_uipElemRecvCounts.clear();
for(unsigned int i = 0; i < m_uipRecvProcs.size(); i++) {
  if(m_uipRecvOffsets[i] < m_uiPreGhostElementSize) {
    //Planning to recieve at least 1 pre-ghost element from this processor
    m_uipElemRecvProcs.push_back(m_uipRecvProcs[i]);
    m_uipElemRecvOffsets.push_back(m_uipRecvOffsets[i]);
    if( (m_uipRecvOffsets[i] + m_uipRecvCounts[i]) <= m_uiPreGhostElementSize) {
      m_uipElemRecvCounts.push_back(m_uipRecvCounts[i]);
    } else {
      m_uipElemRecvCounts.push_back(m_uiPreGhostElementSize - m_uipRecvOffsets[i]);
    }
  }
}//end for i

// Correct the scatter-map so that hanging nodes get correct info ...
//  Logic is to find if any entry in the scatter map is hanging, and if
//  so, we simply correct the scatter map so that it points to the
//  correct anchor instead.

for ( unsigned int i = 0; i < m_uipScatterMap.size(); i++) {
  unsigned int idx = m_uipScatterMap[i];
#ifdef __DEBUG_DA_NLIST__
  assert( (idx >= m_uiElementBegin) && (idx < m_uiPostGhostBegin) );
#endif
  // if idx is an elem ... use nlist ...
  if ( idx < nelem) {
    if (m_ucpLutMasks[2*idx + 1] & 1) {
#ifdef __DEBUG_DA_NLIST__
      assert( !(in[idx].getFlag() & ot::TreeNode::NODE) );
#endif
      // anchor is hanging, so let us correct this ...
      m_uipScatterMap[i] = nlist[8*idx];
    }
  } else {
    // if idx points to a boundary node,
    // get the parent of in[idx] ...
    if ( !(in[idx].getFlag() & ot::TreeNode::NODE)) {
      ot::TreeNode tn = in[idx].getParent();
      while ( in[idx] > tn) {
        idx--;
      }
      //At the end of the above while loop, in[idx] will point to 1 element
      //before the first child of tn.
      idx++;
      m_uipScatterMap[i] = idx;
#ifdef __DEBUG_DA_NLIST__
      assert( in[m_uipScatterMap[i]].getAnchor() == tn.getAnchor() );
      assert( (m_uipScatterMap[i] >= m_uiElementBegin) && (m_uipScatterMap[i] < m_uiPostGhostBegin) );
      assert( in[m_uipScatterMap[i]].getFlag() & ot::TreeNode::NODE );
#endif
    }
  }
}//end for i

#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
if(!m_iRankActive) {
  std::cout<<std::endl;
  std::cout<<"Corrected New Scatter Map for hanging anchors."<<std::endl;
  std::cout<<std::endl;
}
MPI_Barrier(m_mpiCommActive);
#endif

//Correct nlist of all FOREIGNs....
for(unsigned int i=0;i<nelem;i++) {
  if(m_ucpLutMasks[2*i+1] == ot::DA_FLAGS::FOREIGN) {
    for(unsigned int j=0;j<8;j++) {
      nlist[8*i + j] = i;
    }
  }
}//end for i

//~~~~~~~~~~~~~~~~~~~~Mark NODES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool *isNode = NULL;
if(m_uiLocalBufferSize) {
  isNode = new bool[m_uiLocalBufferSize];
}
for (unsigned int i = 0; i < m_uiLocalBufferSize; i++) {
  isNode[i] = false;
}
// loop through the LUT, and tag everybody as being nodes or not ...
for ( unsigned int i = 0; i < nelem; i++) {
  if(m_ucpLutMasks[2*i+1] == ot::DA_FLAGS::FOREIGN) {
    continue;
  }
  for (unsigned int j = 0; j < 8; j++) {
#ifdef __DEBUG_DA_NLIST__
    assert(nlist[8*i + j] < m_uiLocalBufferSize);
#endif
    isNode[nlist[8*i + j]] = true;
  }
}

for (unsigned int i = 0; i < m_uiLocalBufferSize; i++) {
  if ( isNode[i] ) {
#ifdef __DEBUG_DA_NLIST__
    if ( (i >= m_uiElementBegin) && (i < m_uiPostGhostBegin) ) {
      assert( (in[i].getFlag() & ot::TreeNode::NODE) );
    }
#endif
    in[i].orFlag( ot::TreeNode::NODE );
  }
}//end for i

#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
if(!m_iRankActive) {
  std::cout<<std::endl;
  std::cout<<"Finished is NODE..."<<std::endl;
  std::cout<<std::endl;
}
MPI_Barrier(m_mpiCommActive);
#endif

// Now compute the elem node size, and boundary node sizes ...
unsigned int elemNodeSz = 0;
unsigned int bndNodeSz = 0;
unsigned int preGhostNodeSz = 0;
unsigned int preBndNodeSz = 0;
unsigned int postGhostNodeSz = 0;
unsigned int postBndNodeSz = 0;

for (unsigned int i = 0; i < m_uiElementBegin; i++) {
  if ( isNode[i]  && (!(in[i].getFlag() & ot::TreeNode::BOUNDARY)) ) {
    preGhostNodeSz++;                    
  } else if (isNode[i]) {
    preBndNodeSz++;
  }
}

//My Elements which are also nodes.
for (unsigned int i = m_uiElementBegin; i < m_uiElementEnd; i++) {
  if (isNode[i]) {
    elemNodeSz++;        
  }
}

for (unsigned int i = m_uiElementEnd; i < m_uiPostGhostBegin; i++) {
  if (isNode[i]) {
    bndNodeSz++;          
  }
}

for (unsigned int i = m_uiPostGhostBegin; i < m_uiLocalBufferSize; i++) {
  if (isNode[i] && (!(in[i].getFlag() & ot::TreeNode::BOUNDARY)) ) {
    postGhostNodeSz++;
  } else if (isNode[i]) {
    postBndNodeSz++;
  }
}

#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
if(!m_iRankActive) {
  std::cout<<std::endl;
  std::cout<<"Finished setting sizes..."<<std::endl;
  std::cout<<std::endl;
}
MPI_Barrier(m_mpiCommActive);
#endif

if(isNode) {
  delete [] isNode;
  isNode = NULL;
}

//Store the sizes.
m_uiPreGhostNodeSize = preGhostNodeSz;
m_uiPreGhostBoundaryNodeSize = preBndNodeSz;
m_uiPostGhostNodeSize = postGhostNodeSz + postBndNodeSz;
m_uiNodeSize = elemNodeSz;
m_uiBoundaryNodeSize = bndNodeSz;

#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
std::cout<<m_iRankActive<<" Just Before Compression...."<<std::endl<<std::endl;
std::cout<<m_iRankActive<<" "<<m_uiPreGhostNodeSize<<" "<<m_uiPreGhostBoundaryNodeSize
<<" "<<m_uiPostGhostNodeSize<<" "<<m_uiNodeSize<<" "
<<m_uiBoundaryNodeSize<<std::endl<<std::endl;
MPI_Barrier(m_mpiCommActive);
#endif

// All done ... Now COMPRESS
// Now compress the node list ...

// allocate a small buffer to unsort ...
if(m_bCompressLut) {
  m_ucpLutRemainders.resize(8*nelem);
  m_ucpSortOrders.resize(nelem);
}

bool foundBeg = false;
bool foundEnd = false;

m_uiIndependentElementSize = 0;
for (unsigned int i = 0; i < nelem; i++) {
  // get basic info ...
  unsigned int ii = 8*i;

  unsigned int x = in[i].getX();
  unsigned int y = in[i].getY();
  unsigned int z = in[i].getZ();

  if(m_bCompressLut) {
    unsigned int d   = in[i].getLevel();
    unsigned int sz = 1u << (m_uiMaxDepth - d);
    //compute and store the sort order for the non-hanging case.
    if(!(m_ucpLutMasks[2*i+1])){
      unsigned int xp = x + sz;
      unsigned int yp = y + sz;
      unsigned int zp = z + sz;

      unsigned int _x = x^xp;
      unsigned int _y = y^yp;
      unsigned int _z = z^zp;

      if (_x > _y) {
        if ( _y > _z) {
          m_ucpSortOrders[i] = ot::DA_FLAGS::ZYX;
        } else if ( _x > _z ) {
          m_ucpSortOrders[i] = ot::DA_FLAGS::YZX;
        } else {
          m_ucpSortOrders[i] = ot::DA_FLAGS::YXZ;
        }
      } else {
        if ( _x > _z) {
          m_ucpSortOrders[i] = ot::DA_FLAGS::ZXY;
        } else if ( _y > _z ) {
          m_ucpSortOrders[i] = ot::DA_FLAGS::XZY;
        } else {
          m_ucpSortOrders[i] = ot::DA_FLAGS::XYZ;
        }
      }
    }else {
      //store the childnumber instead for the hanging cases.
      unsigned int len_par = (unsigned int)(1u << ( m_uiMaxDepth  - d +1 ) );

      unsigned int a = x % len_par;
      unsigned int b = y % len_par;
      unsigned int c = z % len_par;

      a /= sz;
      b /= sz;
      c /= sz;

      m_ucpSortOrders[i]  = (4*c + 2*b + a);
    }//end sortOrderRegular block
  }//end if Lut compressed

  // use this loop to also detect the begining and end of dependent.
  if ( !foundBeg && !(in[i].getFlag() & ot::DA_FLAGS::DEP_ELEM)
      && (m_ucpLutMasks[2*i+1] != ot::DA_FLAGS::FOREIGN) ) {
    //std::cout << GRN"FOUND BEGINING OF INDEPENDENT "NRM << i << std::endl;
    foundBeg = true;
    m_uiIndependentElementBegin = i;
    m_ptIndependentOffset = Point(x,y,z);
  }

  // reverse loop to find end ...
  if ( !foundEnd && !(in[nelem-i-1].getFlag() & ot::DA_FLAGS::DEP_ELEM)
      && (m_ucpLutMasks[2*(nelem-i-1)+1] != ot::DA_FLAGS::FOREIGN) ) {
    //std::cout << GRN"FOUND END OF INDEPENDENT "NRM << nelem - i << std::endl;
    foundEnd = true;
    m_uiIndependentElementEnd = nelem - i;
  }

  //Actual number of Independent elements. In between IndependentElementBegin
  //and IndependentElementEnd, we can also have dependent elements and so
  //simply taking the difference of end and begin will not work
  if( (!(in[i].getFlag() & ot::DA_FLAGS::DEP_ELEM)) &&
      (m_ucpLutMasks[2*i+1] != ot::DA_FLAGS::FOREIGN) ) {
    m_uiIndependentElementSize++;
  }

  if ( i == m_uiElementBegin ) {
    m_uiElementQuotient = static_cast<unsigned int>(m_uspLutQuotients.size());
  }

  if ( i == m_uiIndependentElementBegin ) {
    m_uiIndependentElementQuotient = static_cast<unsigned int>(m_uspLutQuotients.size());
  }

  //First initialize Masks...
  m_ucpLutMasks[2*i] = 0;

  // locally sort ...
  if(m_bCompressLut) {
    // Perform Golomb-Rice encoding
    //Assumes that the highest 8 bits in offset are not significant, i.e. they will be all 0.
    //This means offset can not be more than (2^24-1) = 16,777,215.
    //Since the number of octants on a processor will not be more than 16M,
    //this does not pose any kind of difficulty.

    unsigned short q;
    unsigned int offset;
    std::vector<unsigned int> nl(8);
    for (unsigned int ij = 0; ij < 8; ij++) {
      nl[ij] = nlist[ii+ij];
    }
    std::sort(nl.begin(), nl.end());

    // 0 is special, it computes the offset wrt i.
    // 0, we have a negative offset .. or 0
    offset = i - nl[0];
    q = (offset >> 8);
    m_ucpLutRemainders[8*i] = (offset%256);
    if (q) {
      m_ucpLutMasks[2*i] |= 1;
      m_uspLutQuotients.push_back(q);
    }

    for (unsigned int j = 1; j < 8; j++) {
      offset =  nl[j] - nl[j-1];
      q = (offset >> 8);
      m_ucpLutRemainders[8*i + j] = (offset%256);
      if (q) {
        m_ucpLutMasks[2*i] |= (1 << j);
        m_uspLutQuotients.push_back(q);
      }
    } // for j ...
    nl.clear();
  }//end if compress

} // for i

//Store Nlist if you are not compressing...
if(!m_bCompressLut) {
  m_uiNlist = nlist; 
}

// free up
nlist.clear();

#ifdef __DEBUG_DA_NLIST__
MPI_Barrier(m_mpiCommActive);
std::cout<<m_iRankActive<<" Just After Compression...."<<std::endl<<std::endl;
std::cout<<m_iRankActive<<" "<<m_uiPreGhostNodeSize<<" "
<<m_uiPreGhostBoundaryNodeSize<<" "<<m_uiPostGhostNodeSize
<<" "<<m_uiNodeSize<<" "<<m_uiBoundaryNodeSize<<std::endl<<std::endl;
MPI_Barrier(m_mpiCommActive);
std::cout << m_iRankActive << ": Leaving " << __func__ << std::endl;
MPI_Barrier(m_mpiCommActive);
#endif

PROF_BUILD_NLIST_END

}//end 4-way BuildNode List
}//end namespace ot


