
/**
  @file matProp.C
  @author Rahul Sampath, rahul.sampath@gmail.com
  */

#include "omg.h"
#include "oda.h" 
#include "omgJac.h"
#include "parUtils.h"

#define SQUARE(x) ((x)*(x))

#define ASSIGN_MAT_PROP_FROM_PTS_ELEM_BLOCK {\
  unsigned int idx = da->curr();\
  unsigned int lev = da->getLevel(idx);\
  Point pt = da->getCurrentOffset();\
  ot::TreeNode me(pt.xint(),pt.yint(),pt.zint(),lev,3,maxD);\
  matPropArr[2*idx] = lapBase; /*Laplacian's coefficient*/\
  matPropArr[2*idx+1] = massBase; /*Mass Matrix's coefficient*/\
  if(ptsCtr < (lapJump.size())) {\
    ot::TreeNode lastVisitedPt;\
    /*Assuming pts are sorted in Morton ordering*/\
    /*Assuming pts and local octants are aligned*/\
    bool repeat = false;\
    do {\
      unsigned int ptX = static_cast<unsigned int>(\
          pts[(3*ptsCtr)]*static_cast<double>(1u << (maxD-1)));\
      unsigned int ptY = static_cast<unsigned int>(\
          pts[(3*ptsCtr)+1]*static_cast<double>(1u << (maxD-1)));\
      unsigned int ptZ = static_cast<unsigned int>(\
          pts[(3*ptsCtr)+2]*static_cast<double>(1u << (maxD-1)));\
      ot::TreeNode currPt(ptX,ptY,ptZ,maxD,3,maxD);\
      repeat = false;\
      if( (currPt < me) && (ptsCtr < (lapJump.size()-1)) ) {\
        repeat = true;\
        ptsCtr++;\
      } else {\
        lastVisitedPt = currPt;\
      }\
    } while(repeat);\
    /*There could be more than 1 pt. per octant*/\
    while( (ptsCtr < lapJump.size()) && (me.isAncestor(lastVisitedPt) ||\
          (me == lastVisitedPt))) {\
      matPropArr[2*idx] += lapJump[ptsCtr];\
      ptsCtr++;\
      if(ptsCtr < lapJump.size()) {\
        unsigned int ptX = static_cast<unsigned int>(\
            pts[(3*ptsCtr)]*static_cast<double>(1u << (maxD-1)));\
        unsigned int ptY = static_cast<unsigned int>(\
            pts[(3*ptsCtr)+1]*static_cast<double>(1u << (maxD-1)));\
        unsigned int ptZ = static_cast<unsigned int>(\
            pts[(3*ptsCtr)+2]*static_cast<double>(1u << (maxD-1)));\
        ot::TreeNode currPt(ptX,ptY,ptZ,maxD,3,maxD);\
        lastVisitedPt = currPt;\
      }\
    }\
  }\
  if(matPropArr[2*idx] > maxCoeff) {\
    maxCoeff = matPropArr[2*idx];\
  }\
  if(matPropArr[2*idx] < minCoeff) {\
    minCoeff = matPropArr[2*idx];\
  }\
}

#define ASSIGN_MAT_PROP_1_CUBE_ELEM_BLOCK {\
  Point pt = da->getCurrentOffset();\
  unsigned int idx = da->curr();\
  unsigned int lev = da->getLevel(idx);\
  double h = ((double)(1u << (maxD - lev)))/((double)(1u << (maxD-1)));\
  double x = ((double)(pt.xint()))/((double)(1u << (maxD-1)));\
  double y = ((double)(pt.yint()))/((double)(1u << (maxD-1)));\
  double z = ((double)(pt.zint()))/((double)(1u << (maxD-1)));\
  matPropArr[2*idx] = 1.0; /*Laplacian's coefficient*/\
  matPropArr[2*idx+1] = massBase; /*Mass Matrix's coefficient*/\
  double minX1 = 0.5-h;\
  double minY1 = 0.0;\
  double minZ1 = 0.0;\
  double maxX1 = 0.5+h;\
  double maxY1 = 1.0;\
  double maxZ1 = 1.0;\
  double minX2 = x;\
  double minY2 = y;\
  double minZ2 = z;\
  double maxX2 = (x+h);\
  double maxY2 = (y+h);\
  double maxZ2 = (z+h);\
  bool oneInTwo = ( (minX1 >= minX2) && (minX1 < maxX2) &&  (minY1 >= minY2) && (minY1 < maxY2) &&  (minZ1 >= minZ2) && (minZ1 < maxZ2) );\
  bool twoInOne = ( (minX2 >= minX1) && (minX2 < maxX1) &&  (minY2 >= minY1) && (minY2 < maxY1) &&  (minZ2 >= minZ1) && (minZ2 < maxZ1) );\
  if( oneInTwo || twoInOne ) {\
    double maxOfMinsX = ( (minX2 > minX1) ? minX2 : minX1 );\
    double maxOfMinsY = ( (minY2 > minY1) ? minY2 : minY1 );\
    double maxOfMinsZ = ( (minZ2 > minZ1) ? minZ2 : minZ1 );\
    double minOfMaxsX = ( (maxX2 < maxX1) ? maxX2 : maxX1 );\
    double minOfMaxsY = ( (maxY2 < maxY1) ? maxY2 : maxY1 );\
    double minOfMaxsZ = ( (maxZ2 < maxZ1) ? maxZ2 : maxZ1 );\
    double intersectionVolume = ((minOfMaxsX - maxOfMinsX)*(minOfMaxsY - maxOfMinsY)*(minOfMaxsZ - maxOfMinsZ))/(h*h*h);\
    if(intersectionVolume <= 0.0) {\
      std::cout<<"intesection Volume: "<<intersectionVolume<<" x: "<<x<<" y: "<<y<<" z: "<<z<<" h: "<<h<<" oneInTwo: "<<oneInTwo<<" twoInOne: "<<twoInOne<<std::endl;\
      assert(false);\
    }\
    matPropArr[2*idx] += lapFac*intersectionVolume; /*Laplacian's coefficient*/\
    matPropArr[2*idx+1]+= massFac*intersectionVolume;/*Mass Matrix's coefficient*/\
  }\
  if(matPropArr[2*idx] > maxCoeff) {\
    maxCoeff = matPropArr[2*idx];\
  }\
  if(matPropArr[2*idx] < minCoeff) {\
    minCoeff = matPropArr[2*idx];\
  }\
}

#define ASSIGN_MAT_PROP_MULTIPLE_CUBES_ELEM_BLOCK {\
  unsigned int idx = da->curr();\
  matPropArr[2*idx] = 1.0; /*Laplacian's coefficient*/\
  matPropArr[2*idx+1] = massBase; /*Mass Matrix's coefficient*/\
  unsigned int procElementSize = da->getElementSize();\
  unsigned int numCubesInterval = procElementSize/numCubes;\
  if( (idx % numCubesInterval) == 0 ) {\
    matPropArr[2*idx] += lapFac; /*Laplacian's coefficient*/\
    matPropArr[2*idx+1] +=  massFac; /*Mass Matrix's coefficient*/\
  }\
  if(matPropArr[2*idx] > maxCoeff) {\
    maxCoeff = matPropArr[2*idx];\
  }\
  if(matPropArr[2*idx] < minCoeff) {\
    minCoeff = matPropArr[2*idx];\
  }\
}

#define ASSIGN_MAT_PROP_CHECKER_BOARD_ELEM_BLOCK {\
  Point pt = da->getCurrentOffset();\
  unsigned int idx = da->curr();\
  unsigned int lev = da->getLevel(idx);\
  double h = ((double)(1u << (maxD - lev)))/((double)(1u << (maxD-1)));\
  double x = ((double)(pt.xint()))/((double)(1u << (maxD-1)));\
  double y = ((double)(pt.yint()))/((double)(1u << (maxD-1)));\
  double z = ((double)(pt.zint()))/((double)(1u << (maxD-1)));\
  matPropArr[2*idx+1] = massBase; /*Mass Matrix's coefficient*/\
  if(z >= 0.0 && z < 0.5 ) {\
    if(x >= 0.0 && x < 0.5 && y >= 0.0 && y < 0.5) {\
      matPropArr[2*idx] = lapFac; /*Laplacian's coefficient*/\
    } else if(x >= 0.5 && x < 1.0 && y >= 0.0 && y < 0.5) {\
      matPropArr[2*idx] = 1.0; /*Laplacian's coefficient*/\
    } else if(x >= 0.0 && x < 0.5 && y >= 0.5 && y < 1.0) {\
      matPropArr[2*idx] = 1.0; /*Laplacian's coefficient*/\
    } else {\
      matPropArr[2*idx] = lapFac; /*Laplacian's coefficient*/\
    }\
  } else {\
    if(x >= 0.0 && x < 0.5 && y >= 0.0 && y < 0.5) {\
      matPropArr[2*idx] = 1.0; /*Laplacian's coefficient*/\
    } else if(x >= 0.5 && x < 1.0 && y >= 0.0 && y < 0.5) {\
      matPropArr[2*idx] = lapFac; /*Laplacian's coefficient*/\
    } else if(x >= 0.0 && x < 0.5 && y >= 0.5 && y < 1.0) {\
      matPropArr[2*idx] = lapFac; /*Laplacian's coefficient*/\
    } else {\
      matPropArr[2*idx] = 1.0; /*Laplacian's coefficient*/\
    }\
  }\
  if(matPropArr[2*idx] > maxCoeff) {\
    maxCoeff = matPropArr[2*idx];\
  }\
  if(matPropArr[2*idx] < minCoeff) {\
    minCoeff = matPropArr[2*idx];\
  }\
}

#define ASSIGN_MAT_PROP_ANALYTIC_FN_ELEM_BLOCK {\
  Point pt = da->getCurrentOffset();\
  unsigned int idx = da->curr();\
  unsigned int lev = da->getLevel(idx);\
  double h = ((double)(1u << (maxD - lev)))/((double)(1u << (maxD-1)));\
  double x = ((double)(pt.xint()))/((double)(1u << (maxD-1)));\
  double y = ((double)(pt.yint()))/((double)(1u << (maxD-1)));\
  double z = ((double)(pt.zint()))/((double)(1u << (maxD-1)));\
  matPropArr[2*idx+1] = massBase; /*Mass Matrix's coefficient*/\
  /*Laplacian's coefficient*/\
  matPropArr[2*idx] = 1.0 +\
  (lapFac*(SQUARE(cos(lapFreq*M_PI*(x+(0.5*h))))+\
           SQUARE(cos(lapFreq*M_PI*(y+(0.5*h))))+\
           SQUARE(cos(lapFreq*M_PI*(z+(0.5*h))))));\
  if(matPropArr[2*idx] > maxCoeff) {\
    maxCoeff = matPropArr[2*idx];\
  }\
  if(matPropArr[2*idx] < minCoeff) {\
    minCoeff = matPropArr[2*idx];\
  }\
}

#define ASSIGN_MAT_PROP_FINE_TO_COARSE_ELEM_BLOCK {\
  /*The fine loop is always Writable, but the coarse loop*/\
  /*could be Independent or W_Dependent. Hence the fine counter must*/\
  /*be incremented properly to align with the coarse.*/\
  unsigned int idxC= da->curr();\
  Point Cpt = da->getCurrentOffset();\
  assert(daf->curr() < daf->end<ot::DA_FLAGS::WRITABLE>());\
  while(daf->getCurrentOffset() != Cpt) {\
    daf->next<ot::DA_FLAGS::WRITABLE>();\
    assert(daf->curr() < daf->end<ot::DA_FLAGS::WRITABLE>());\
  }\
  if(daf->getLevel(daf->curr()) == da->getLevel(idxC)) {\
    /*The coarse and fine elements are the same,*/\
    assert(daf->curr() < daf->end<ot::DA_FLAGS::WRITABLE>());\
    unsigned int idxF = daf->curr();\
    matPropArr[2*idxC] = fMatArr[2*idxF];\
    matPropArr[2*idxC+1] = fMatArr[2*idxF+1];\
    daf->next<ot::DA_FLAGS::WRITABLE>();\
  }else {\
    double cLapVal = 0.0;\
    double cMassVal = 0.0;\
    for(unsigned char cNumFine = 0; cNumFine < 8; cNumFine++) {\
      /*The coarse and fine elements are NOT the same. */\
      /*Loop over each of the 8 children of the coarse element.*/\
      /*These are the underlying fine elements.*/\
      assert(daf->curr() < daf->end<ot::DA_FLAGS::WRITABLE>());\
      unsigned int idxF = daf->curr();\
      cLapVal += fMatArr[2*idxF];\
      cMassVal += fMatArr[2*idxF+1];\
      daf->next<ot::DA_FLAGS::WRITABLE>();\
    }\
    matPropArr[2*idxC] = (cLapVal/8.0);\
    matPropArr[2*idxC+1] = (cMassVal/8.0);\
  }\
  if(matPropArr[2*idxC] > maxCoeff) {\
    maxCoeff = matPropArr[2*idxC];\
  }\
  if(matPropArr[2*idxC] < minCoeff) {\
    minCoeff = matPropArr[2*idxC];\
  }\
}

#define CHK_AND_SCATTER_FINE_TO_COARSE_BLOCK {\
  if( damg_i->da_aux == NULL ) {\
    daf = da;\
    fMatPropVec = matPropVecPtr;\
    changedPartition  = false;\
  }else {\
    daf = damg_i->da_aux;\
    /*Need to Scatter Values*/\
    fMatPropVec = new std::vector<double>;\
    changedPartition  = true;\
    /*elemental - non-ghosted*/\
    std::vector<double> tmpVec1;\
    da->createVector<double>(tmpVec1,true,false,2);\
    double *vec1Arr = NULL;\
    /*Elemental,non-Ghosted,Write,2 Dof.*/\
    da->vecGetBuffer<double>(tmpVec1,vec1Arr,true,false,false,2);\
    matPropArr = NULL;\
    /*Elemental,Ghosted,Read-only,2 Dof.*/\
    da->vecGetBuffer<double>((*matPropVecPtr), matPropArr, true, true, true, 2);\
    if(da->iAmActive()) {\
      for(da->init<ot::DA_FLAGS::WRITABLE>();\
          da->curr() < da->end<ot::DA_FLAGS::WRITABLE>();\
          da->next<ot::DA_FLAGS::WRITABLE>()) {\
        unsigned int idx = da->curr();\
        vec1Arr[2*idx] = matPropArr[2*idx];\
        vec1Arr[2*idx+1] = matPropArr[2*idx+1];\
      }\
    }\
    da->vecRestoreBuffer<double>((*matPropVecPtr), matPropArr, true, true, true, 2);\
    da->vecRestoreBuffer<double>(tmpVec1,vec1Arr,true,false,false,2);\
    par::scatterValues<double>(tmpVec1, (*fMatPropVec), (2*(daf->getElementSize())), da->getComm());\
    tmpVec1.clear();\
  }\
}

#define SET_SINGLE_LEVEL_MAT_PROP_BLOCK(ELEM_MAT_PROP_BLOCK, EXTRA_INITS)  {\
  damg_i->user = ctx;\
  da = damg_i->da;\
  comm = da->getCommActive();\
  /*Elem,Ghosted, 2-dof vector.*/\
  /*Note: I am creating a ghosted vector only*/\
  /*because the mat-vec will need it.*/\
  /*So this way, I can avoid mallocs inside the mat-vec.*/\
  da->createVector<double>((*matPropVecPtr), true, true, 2);\
  for(unsigned int i = 0; i < matPropVecPtr->size(); i++) {\
    (*matPropVecPtr)[i] = 0.0;\
  }\
  /*Elemental,Ghosted,Write,2 Dof.*/\
  da->vecGetBuffer<double>((*matPropVecPtr), matPropArr,\
      true, true, false, 2);\
  maxD =  da->getMaxDepth();\
  if(da->iAmActive()) {\
    MPI_Comm_rank(comm,&rank);\
    double maxCoeff = 0.0;\
    double minCoeff = 1.0e+8;\
    double globalMaxCoeff;\
    double globalMinCoeff;\
    /*Any independent element has all its 8*/\
    /*indices pointing to elements the same processor owns.*/\
    /*So, even if any independent element is sent to*/\
    /*another processor as a ghost (This is rare, but might*/\
    /*happen due to a-priori communication and partial*/\
    /*second ring communication), it would have to*/\
    /*be FOREIGN on that processor. So overlapping comm*/\
    /*and comp is not a problem.*/\
    EXTRA_INITS\
    for(da->init<ot::DA_FLAGS::W_DEPENDENT>();\
        da->curr() < da->end<ot::DA_FLAGS::W_DEPENDENT>();\
        da->next<ot::DA_FLAGS::W_DEPENDENT>()) {\
      ELEM_MAT_PROP_BLOCK\
    }\
    da->ReadFromGhostElemsBegin<double>(matPropArr, 2);\
    EXTRA_INITS\
    for(da->init<ot::DA_FLAGS::INDEPENDENT>();\
        da->curr() < da->end<ot::DA_FLAGS::INDEPENDENT>();\
        da->next<ot::DA_FLAGS::INDEPENDENT>()) {\
      ELEM_MAT_PROP_BLOCK\
    }\
    da->ReadFromGhostElemsEnd<double>(matPropArr);\
    par::Mpi_Reduce<double>(&maxCoeff, &globalMaxCoeff, 1, MPI_MAX, 0, comm);\
    par::Mpi_Reduce<double>(&minCoeff, &globalMinCoeff, 1, MPI_MIN, 0, comm);\
    if(!rank) {\
      std::cout<<"Level: "<<(nlevels-(damg_i->nlevels))<<\
      " Max Lap. Coeff: "<<globalMaxCoeff<<\
      " Min Lap. Coeff: "<<globalMinCoeff<<std::endl;\
    }\
    MPI_Barrier(comm);\
  } /*end if active*/\
  da->vecRestoreBuffer<double>((*matPropVecPtr), matPropArr,\
      true, true, false, 2);\
}

#define FINE_TO_COARSE_BLOCK {\
  damg_i = damg[i];\
  damg_i->user = ctx;\
  da = damg_i->da;\
  assert(da->iAmActive() == daf->iAmActive());\
  da->createVector<double>((*matPropVecPtr), true, true, 2);\
  for(unsigned int j = 0; j < matPropVecPtr->size(); j++) {\
    (*matPropVecPtr)[j] = 0.0;\
  }\
  comm = da->getCommActive();\
  /*Elemental,Ghosted,Write,2 Dof.*/\
  matPropArr = NULL;\
  da->vecGetBuffer<double>((*matPropVecPtr), matPropArr,\
      true, true, false, 2);\
  double *fMatArr = NULL;\
  if(changedPartition) {\
    /*Elemental, non-Ghosted, Read-only, 2 Dof.*/\
    daf->vecGetBuffer<double>((*fMatPropVec), fMatArr,\
        true, false, true, 2);\
  }else {\
    /*Elemental, Ghosted, Read-only, 2 Dof.*/\
    daf->vecGetBuffer<double>((*fMatPropVec), fMatArr,\
        true, true, true, 2);\
  }\
  if(da->iAmActive()) {\
    MPI_Comm_rank(comm,&rank);\
    double maxCoeff = 0.0;\
    double minCoeff = 1.0e+8;\
    double globalMaxCoeff;\
    double globalMinCoeff;\
    /*Loop through the coarse and fine simultaneously*/\
    /*Note: If Coarse is Independent, then the*/\
    /*corresponding Fine is also independent.*/\
    /*Hence, overlapping comm with comp is possible.*/\
    /*First, we loop though the dependent elements.*/\
    /*Then we begin the communication and simulatenously*/\
    /*loop over the independent elements.*/\
    for(da->init<ot::DA_FLAGS::W_DEPENDENT>(),\
        daf->init<ot::DA_FLAGS::WRITABLE>();\
        da->curr() < da->end<ot::DA_FLAGS::W_DEPENDENT>();\
        da->next<ot::DA_FLAGS::W_DEPENDENT>()) {\
      ASSIGN_MAT_PROP_FINE_TO_COARSE_ELEM_BLOCK \
    } /*end dependent loop*/\
    da->ReadFromGhostElemsBegin<double>(matPropArr,2);\
    for(da->init<ot::DA_FLAGS::INDEPENDENT>(),\
        daf->init<ot::DA_FLAGS::WRITABLE>();\
        da->curr() < da->end<ot::DA_FLAGS::INDEPENDENT>();\
        da->next<ot::DA_FLAGS::INDEPENDENT>()) {\
      ASSIGN_MAT_PROP_FINE_TO_COARSE_ELEM_BLOCK \
    } /*end Independent loop */\
    da->ReadFromGhostElemsEnd<double>(matPropArr);\
    par::Mpi_Reduce<double>(&maxCoeff, &globalMaxCoeff, 1, MPI_MAX, 0, comm);\
    par::Mpi_Reduce<double>(&minCoeff, &globalMinCoeff, 1, MPI_MIN, 0, comm);\
    if(!rank) {\
      std::cout<<"Level: "<<i<<" Max Lap. Coeff: "\
      <<globalMaxCoeff<<" Min Lap. Coeff: "\
      <<globalMinCoeff<<std::endl;\
    }\
    MPI_Barrier(comm);\
  } /*end check if active*/\
  da->vecRestoreBuffer<double>((*matPropVecPtr),\
      matPropArr, true, true, false, 2);\
  if(changedPartition) {\
    /*Elemental, non-Ghosted, Read-only, 2 Dof.*/\
    daf->vecRestoreBuffer<double>((*fMatPropVec),\
        fMatArr, true, false, true, 2);\
  }else {\
    /*Elemental, Ghosted, Read-only, 2 Dof.*/\
    daf->vecRestoreBuffer<double>((*fMatPropVec),\
        fMatArr, true, true, true, 2);\
  }\
}

void SetUserContexts(ot::DAMG* damg) {

  PetscInt       jacType = 1;
  PetscOptionsGetInt(0, "-jacType", &jacType, 0);

  if(jacType == 1) { return; }

  PetscTruth setMatPropsAtCoarsest;
  PetscOptionsHasName(0,"-setMatPropsAtCoarsest",&setMatPropsAtCoarsest);

  if(setMatPropsAtCoarsest) {
    assert(jacType == 2);
    SetUserContextsCoarsestToFinest(damg);
    return;
  }

  int       nlevels = damg[0]->nlevels; //number of multigrid levels

  //Set Mat Props Finest to Coarsest...
  void * ctx;

  // Set for the finest level first
  if(jacType == 2) {
    ctx = new Jac2MFreeData;
  }else {
    ctx = new Jac3MFreeData;
  }

  ot::DAMG damg_i = damg[nlevels-1];
  int rank;
  MPI_Comm comm;

  std::vector<double> * matPropVecPtr = new std::vector<double>;

  if(jacType == 2) {
    (static_cast<Jac2MFreeData*>(ctx))->matProp = matPropVecPtr;
    (static_cast<Jac2MFreeData*>(ctx))->isFinestLevel = true;
    (static_cast<Jac2MFreeData*>(ctx))->Jmat_private = NULL;
    (static_cast<Jac2MFreeData*>(ctx))->inTmp = NULL;
    (static_cast<Jac2MFreeData*>(ctx))->outTmp = NULL;
  }else {
    (static_cast<Jac3MFreeData*>(ctx))->matProp = matPropVecPtr;
    (static_cast<Jac3MFreeData*>(ctx))->isFinestLevel = true;
    (static_cast<Jac3MFreeData*>(ctx))->isCoarsestLevel = false;
    (static_cast<Jac3MFreeData*>(ctx))->daf = NULL;
    (static_cast<Jac3MFreeData*>(ctx))->matPropFine = NULL;
    (static_cast<Jac3MFreeData*>(ctx))->changedPartition = false;
    (static_cast<Jac3MFreeData*>(ctx))->JmatThisLevel = NULL;
    (static_cast<Jac3MFreeData*>(ctx))->BmatThisLevel = NULL;
    (static_cast<Jac3MFreeData*>(ctx))->Jmat_private = NULL;
    (static_cast<Jac3MFreeData*>(ctx))->inTmp = NULL;
    (static_cast<Jac3MFreeData*>(ctx))->outTmp = NULL;
  }

  //Default values is the same as the const. coeff. case 
  PetscReal lapFac = 0.0;
  PetscReal massFac = 0.0;
  PetscReal massBase = 1.0;
  PetscTruth optFound;
  PetscOptionsGetReal("lap","-MatPropFac",&lapFac,&optFound);
  PetscOptionsGetReal("mass","-MatPropFac",&massFac,&optFound);
  PetscOptionsGetReal("mass","-BaseMatProp",&massBase,&optFound);
  unsigned int maxD ;
  double *matPropArr = NULL;
  ot::DA* da;

  PetscTruth setMatPropFromAnalyticFn;
  PetscOptionsHasName(0,"-setMatPropFromAnalyticFn",&setMatPropFromAnalyticFn);

  if(setMatPropFromAnalyticFn) {
    PetscReal lapFreq = 1.0;
    PetscOptionsGetReal(0,"-lapFreq",&lapFreq,&optFound);
    int dummyInit;
    SET_SINGLE_LEVEL_MAT_PROP_BLOCK(\
        ASSIGN_MAT_PROP_ANALYTIC_FN_ELEM_BLOCK, dummyInit = 0;) 
  } else {
    PetscTruth setCheckerBoardMatProp;
    PetscOptionsHasName(0,"-setCheckerBoardMatProp",&setCheckerBoardMatProp);
    if(setCheckerBoardMatProp) {
      int dummyInit;
      SET_SINGLE_LEVEL_MAT_PROP_BLOCK(\
          ASSIGN_MAT_PROP_CHECKER_BOARD_ELEM_BLOCK, dummyInit = 0;) 
    } else {
      PetscInt numCubes = 1;
      PetscOptionsGetInt(0,"-numCubes",&numCubes,0);
      if(numCubes == 1) {
        int dummyInit;
        SET_SINGLE_LEVEL_MAT_PROP_BLOCK(\
            ASSIGN_MAT_PROP_1_CUBE_ELEM_BLOCK, dummyInit = 0;) 
      } else {
        int dummyInit;
        SET_SINGLE_LEVEL_MAT_PROP_BLOCK(\
            ASSIGN_MAT_PROP_MULTIPLE_CUBES_ELEM_BLOCK, dummyInit = 0;) 
      }
    }
  }

  ot::DA* daf;
  std::vector<double>* fMatPropVec = NULL;
  bool changedPartition;

  if(nlevels > 1) {
    CHK_AND_SCATTER_FINE_TO_COARSE_BLOCK 
  }

  //Coarser levels
  for(int i = (nlevels-2); i >= 0; i--) {
    if(jacType == 2) {
      ctx = new Jac2MFreeData;
    }else {
      ctx = new Jac3MFreeData;
    }

    matPropVecPtr = new std::vector<double>;

    if(jacType == 2) {
      (static_cast<Jac2MFreeData*>(ctx))->matProp = matPropVecPtr;
      (static_cast<Jac2MFreeData*>(ctx))->isFinestLevel = false;
      (static_cast<Jac2MFreeData*>(ctx))->Jmat_private = NULL;
      (static_cast<Jac2MFreeData*>(ctx))->inTmp = NULL;
      (static_cast<Jac2MFreeData*>(ctx))->outTmp = NULL;
    }else {
      (static_cast<Jac3MFreeData*>(ctx))->matProp = matPropVecPtr;
      (static_cast<Jac3MFreeData*>(ctx))->isFinestLevel = false;
      (static_cast<Jac3MFreeData*>(ctx))->isCoarsestLevel = false;
      (static_cast<Jac3MFreeData*>(ctx))->daf = daf;
      //Note, while using matPropFine, first check changedPartition. If the
      //partition was changed, then the vector is a non-ghosted elemental
      //vector. Else it is a ghosted elemental vector.
      (static_cast<Jac3MFreeData*>(ctx))->matPropFine = fMatPropVec;
      //If the partition was changed a new vector would have been created. Else
      // only the pointer is copied.
      (static_cast<Jac3MFreeData*>(ctx))->changedPartition = changedPartition;
      (static_cast<Jac3MFreeData*>(ctx))->JmatThisLevel = NULL;
      (static_cast<Jac3MFreeData*>(ctx))->BmatThisLevel = NULL;
      (static_cast<Jac3MFreeData*>(ctx))->Jmat_private = NULL;
      (static_cast<Jac3MFreeData*>(ctx))->inTmp = NULL;
      (static_cast<Jac3MFreeData*>(ctx))->outTmp = NULL;
    }

    FINE_TO_COARSE_BLOCK

      if(changedPartition && (jacType == 2)) {
        fMatPropVec->clear();
        delete fMatPropVec;
      }

    if(i) {	
      CHK_AND_SCATTER_FINE_TO_COARSE_BLOCK 
    }
  }//end for i

  //To handle the case of a single grid, it is best to set the coarsest level's
  //flag separately at the end.
  if(jacType == 3) {
    (static_cast<Jac3MFreeData*>(damg[0]->user))->isCoarsestLevel = true;
  }
}//end fn.

void SetUserContextsFromPts(ot::DAMG* damg,
    const std::vector<double>& pts,
    const std::vector<double> & lapJump) {

  PetscInt       jacType = 1;
  PetscOptionsGetInt(0, "-jacType", &jacType, 0);

  assert(jacType != 1);

  assert(pts.size() == (3*lapJump.size()));

  int nlevels = damg[0]->nlevels; //number of mg levels
  ot::DAMG damg_i = damg[nlevels-1];

  //Set Mat Props Finest to Coarsest...
  void * ctx;

  // Set for the finest level first
  if(jacType == 2) {
    ctx = new Jac2MFreeData;
  }else {
    ctx = new Jac3MFreeData;
  }

  std::vector<double> * matPropVecPtr = new std::vector<double>;

  if(jacType == 2) {
    (static_cast<Jac2MFreeData*>(ctx))->matProp = matPropVecPtr;
    (static_cast<Jac2MFreeData*>(ctx))->isFinestLevel = true;
    (static_cast<Jac2MFreeData*>(ctx))->Jmat_private = NULL;
    (static_cast<Jac2MFreeData*>(ctx))->inTmp = NULL;
    (static_cast<Jac2MFreeData*>(ctx))->outTmp = NULL;
  }else {
    (static_cast<Jac3MFreeData*>(ctx))->matProp = matPropVecPtr;
    (static_cast<Jac3MFreeData*>(ctx))->isFinestLevel = true;
    (static_cast<Jac3MFreeData*>(ctx))->isCoarsestLevel = false;
    (static_cast<Jac3MFreeData*>(ctx))->daf = NULL;
    (static_cast<Jac3MFreeData*>(ctx))->matPropFine = NULL;
    (static_cast<Jac3MFreeData*>(ctx))->changedPartition = false;
    (static_cast<Jac3MFreeData*>(ctx))->JmatThisLevel = NULL;
    (static_cast<Jac3MFreeData*>(ctx))->BmatThisLevel = NULL;
    (static_cast<Jac3MFreeData*>(ctx))->Jmat_private = NULL;
    (static_cast<Jac3MFreeData*>(ctx))->inTmp = NULL;
    (static_cast<Jac3MFreeData*>(ctx))->outTmp = NULL;
  }

  //Default values is the same as the const. coeff. case 
  PetscReal lapBase = 1.0;
  PetscReal massBase = 1.0;
  PetscTruth optFound;
  PetscOptionsGetReal("lap","-BaseMatProp",&lapBase,&optFound);
  PetscOptionsGetReal("mass","-BaseMatProp",&massBase,&optFound);
  unsigned int maxD ;
  double *matPropArr = NULL;
  ot::DA* da;
  int rank;
  MPI_Comm comm;

  unsigned int ptsCtr;
  SET_SINGLE_LEVEL_MAT_PROP_BLOCK(ASSIGN_MAT_PROP_FROM_PTS_ELEM_BLOCK,
      ptsCtr = 0;) 

    //The coarser levels...
    ot::DA* daf;
  std::vector<double>* fMatPropVec = NULL;
  bool changedPartition;

  if(nlevels > 1) {
    CHK_AND_SCATTER_FINE_TO_COARSE_BLOCK 
  }

  //Coarser levels
  for(int i = (nlevels-2); i >= 0; i--) {
    if(jacType == 2) {
      ctx = new Jac2MFreeData;
    }else {
      ctx = new Jac3MFreeData;
    }

    matPropVecPtr = new std::vector<double>;

    if(jacType == 2) {
      (static_cast<Jac2MFreeData*>(ctx))->matProp = matPropVecPtr;
      (static_cast<Jac2MFreeData*>(ctx))->isFinestLevel = false;
      (static_cast<Jac2MFreeData*>(ctx))->Jmat_private = NULL;
      (static_cast<Jac2MFreeData*>(ctx))->inTmp = NULL;
      (static_cast<Jac2MFreeData*>(ctx))->outTmp = NULL;
    }else {
      (static_cast<Jac3MFreeData*>(ctx))->matProp = matPropVecPtr;
      (static_cast<Jac3MFreeData*>(ctx))->isFinestLevel = false;
      (static_cast<Jac3MFreeData*>(ctx))->isCoarsestLevel = false;
      (static_cast<Jac3MFreeData*>(ctx))->daf = daf;
      //Note, while using matPropFine, first check changedPartition. If the
      //partition was changed, then the vector is a non-ghosted elemental
      //vector. Else it is a ghosted elemental vector.
      (static_cast<Jac3MFreeData*>(ctx))->matPropFine = fMatPropVec;
      //If the partition was changed a new vector would have been created. Else
      // only the pointer is copied.
      (static_cast<Jac3MFreeData*>(ctx))->changedPartition = changedPartition;
      (static_cast<Jac3MFreeData*>(ctx))->JmatThisLevel = NULL;
      (static_cast<Jac3MFreeData*>(ctx))->BmatThisLevel = NULL;
      (static_cast<Jac3MFreeData*>(ctx))->Jmat_private = NULL;
      (static_cast<Jac3MFreeData*>(ctx))->inTmp = NULL;
      (static_cast<Jac3MFreeData*>(ctx))->outTmp = NULL;
    }

    FINE_TO_COARSE_BLOCK

      if(changedPartition && (jacType == 2)) {
        fMatPropVec->clear();
        delete fMatPropVec;
      }

    if(i) {	
      CHK_AND_SCATTER_FINE_TO_COARSE_BLOCK 
    }
  }//end for i

  //To handle the case of a single grid, it is best to set the coarsest level's
  //flag separately at the end.
  if(jacType == 3) {
    (static_cast<Jac3MFreeData*>(damg[0]->user))->isCoarsestLevel = true;
  }
}//end fn.

#undef ASSIGN_MAT_PROP_FINE_TO_COARSE_ELEM_BLOCK 
#undef CHK_AND_SCATTER_FINE_TO_COARSE_BLOCK 
#undef FINE_TO_COARSE_BLOCK

#define ASSIGN_MAT_PROP_COARSE_TO_FINE_ELEM_BLOCK {\
  /*The fine loop is always Writable, but the coarse loop*/\
  /*could be Independent or W_Dependent. Hence the fine counter must*/\
  /*be incremented properly to align with the coarse.*/\
  unsigned int idxC= dac->curr();\
  Point Cpt = dac->getCurrentOffset();\
  assert(daf->curr() < daf->end<ot::DA_FLAGS::WRITABLE>());\
  while(daf->getCurrentOffset() != Cpt) {\
    daf->next<ot::DA_FLAGS::WRITABLE>();\
    assert(daf->curr() < daf->end<ot::DA_FLAGS::WRITABLE>());\
  }\
  if(daf->getLevel(daf->curr()) == dac->getLevel(idxC)) {\
    /*The coarse and fine elements are the same,*/\
    assert(daf->curr() < daf->end<ot::DA_FLAGS::WRITABLE>());\
    unsigned int idxF = daf->curr();\
    fMatArr[2*idxF] = cMatArr[2*idxC];\
    fMatArr[(2*idxF)+1] = cMatArr[(2*idxC)+1];\
    if(fMatArr[2*idxF] > maxCoeff) {\
      maxCoeff = fMatArr[2*idxF];\
    }\
    if(fMatArr[2*idxF] < minCoeff) {\
      minCoeff = fMatArr[2*idxF];\
    }\
    daf->next<ot::DA_FLAGS::WRITABLE>();\
  }else {\
    for(unsigned char cNumFine = 0; cNumFine < 8; cNumFine++) {\
      /*The coarse and fine elements are NOT the same. */\
      /*Loop over each of the 8 children of the coarse element.*/\
      /*These are the underlying fine elements.*/\
      assert(daf->curr() < daf->end<ot::DA_FLAGS::WRITABLE>());\
      unsigned int idxF = daf->curr();\
      fMatArr[2*idxF] = cMatArr[2*idxC];\
      fMatArr[(2*idxF)+1] = cMatArr[(2*idxC)+1];\
      if(fMatArr[2*idxF] > maxCoeff) {\
        maxCoeff = fMatArr[2*idxF];\
      }\
      if(fMatArr[2*idxF] < minCoeff) {\
        minCoeff = fMatArr[2*idxF];\
      }\
      daf->next<ot::DA_FLAGS::WRITABLE>();\
    }\
  }\
}

#define COARSE_TO_FINE_BLOCK {\
  /*The finer levels... */\
  for(int i = 1; i < nlevels; i++) {\
    ot::DA* dac = damg[i-1]->da;\
    ot::DA* daf;\
    bool changedPartition;\
    if(damg[i]->da_aux) {\
      daf = damg[i]->da_aux;\
      changedPartition = true;\
    }else {\
      daf = damg[i]->da;\
      changedPartition = false;\
    }\
    assert(dac->iAmActive() == daf->iAmActive());\
    comm = daf->getCommActive();\
    if(daf->iAmActive()) {\
      MPI_Comm_rank(comm,&rank);\
    }\
    std::vector<double>* cMatPropVec = \
    (static_cast<Jac2MFreeData*>(damg[i-1]->user))->matProp;\
    std::vector<double>* fMatPropVec = new std::vector<double>;\
    if(changedPartition) {\
      /*Elemental and Non-Ghosted*/\
      daf->createVector<double>((*fMatPropVec), true, false, 2);\
    } else {\
      /*Elemental and Ghosted*/\
      daf->createVector<double>((*fMatPropVec), true, true, 2);\
    }\
    for(unsigned int j = 0; j < fMatPropVec->size(); j++) {\
      (*fMatPropVec)[j] = 0.0;\
    }\
    double *cMatArr = NULL;\
    double *fMatArr = NULL;\
    double maxCoeff = 0.0;\
    double minCoeff = 1.0e+8;\
    double globalMaxCoeff;\
    double globalMinCoeff;\
    /*Read-only buffer*/\
    dac->vecGetBuffer<double>((*cMatPropVec), cMatArr,true, true, true, 2);\
    if(changedPartition) {\
      if(daf->iAmActive()) {\
        /*Writable buffer*/\
        daf->vecGetBuffer<double>((*fMatPropVec), fMatArr,\
            true, false, false, 2);\
        /*Can not overlap comm and comp here. So direct WRITABLE loop*/\
        for(dac->init<ot::DA_FLAGS::WRITABLE>(),\
            daf->init<ot::DA_FLAGS::WRITABLE>();\
            dac->curr() < dac->end<ot::DA_FLAGS::WRITABLE>();\
            dac->next<ot::DA_FLAGS::WRITABLE>()) {\
          ASSIGN_MAT_PROP_COARSE_TO_FINE_ELEM_BLOCK \
        } /*end loop*/\
        daf->vecRestoreBuffer<double>((*fMatPropVec), fMatArr,\
            true, false, false, 2);\
      }/*end if active*/\
      std::vector<double> tmpVecForScatter;\
      /*Scatter from fMatPropVec to*/\
      /*tmpVecForScatter (created in the function)*/\
      par::scatterValues<double>((*fMatPropVec), tmpVecForScatter,\
          (2*(damg[i]->da->getElementSize())),\
          daf->getComm());\
      fMatPropVec->clear();\
      delete fMatPropVec;\
      std::vector<double> *tmpFineVec = new std::vector<double>;\
      /*Elemental and Ghosted*/\
      damg[i]->da->createVector<double>((*tmpFineVec), true, true, 2);\
      if(damg[i]->da->iAmActive()) {\
        /*Writable buffer*/\
        damg[i]->da->vecGetBuffer<double>((*tmpFineVec), fMatArr,\
            true, true, false, 2);\
        /*Read-only buffer*/\
        double *tmpFmatArr = NULL;\
        damg[i]->da->vecGetBuffer<double>(tmpVecForScatter, tmpFmatArr,\
            true, false, true, 2);\
        /*Copy from tmpVecForScatter to tmpFineVec*/\
        /*W_DEPENDENT loop*/\
        for(damg[i]->da->init<ot::DA_FLAGS::W_DEPENDENT>();\
            damg[i]->da->curr() < damg[i]->da->end<ot::DA_FLAGS::W_DEPENDENT>();\
            damg[i]->da->next<ot::DA_FLAGS::W_DEPENDENT>()) {\
          unsigned int idxCurr = damg[i]->da->curr();\
          fMatArr[2*idxCurr] = tmpFmatArr[2*idxCurr]; \
          fMatArr[(2*idxCurr) + 1] = tmpFmatArr[(2*idxCurr) + 1]; \
        } /*end dependent loop*/\
        /*Begin read from ghost elements on the da grid*/      \
        damg[i]->da->ReadFromGhostElemsBegin<double>(fMatArr,2);\
        /*Overlap communication with Independent loop*/\
        for(damg[i]->da->init<ot::DA_FLAGS::INDEPENDENT>();\
            damg[i]->da->curr() < damg[i]->da->end<ot::DA_FLAGS::INDEPENDENT>();\
            damg[i]->da->next<ot::DA_FLAGS::INDEPENDENT>()) {\
          unsigned int idxCurr = damg[i]->da->curr();\
          fMatArr[2*idxCurr] = tmpFmatArr[2*idxCurr]; \
          fMatArr[(2*idxCurr) + 1] = tmpFmatArr[(2*idxCurr) + 1]; \
        } /*end Independent loop */\
        /*End read from ghost elements on the da grid*/      \
        damg[i]->da->ReadFromGhostElemsEnd<double>(fMatArr);\
        damg[i]->da->vecRestoreBuffer<double>(tmpVecForScatter, tmpFmatArr,\
            true, false, true, 2);\
        tmpVecForScatter.clear();\
        damg[i]->da->vecRestoreBuffer<double>((*tmpFineVec), fMatArr,\
            true, true, false, 2);\
      }/*end if active*/\
      fMatPropVec = tmpFineVec;\
      tmpFineVec = NULL;      \
    } else {\
      if(daf->iAmActive()) {\
        /*Writable buffer*/\
        daf->vecGetBuffer<double>((*fMatPropVec), fMatArr,\
            true, true, false, 2);\
        /*W_DEPENDENT loop*/\
        for(dac->init<ot::DA_FLAGS::W_DEPENDENT>(),\
            daf->init<ot::DA_FLAGS::WRITABLE>();\
            dac->curr() < dac->end<ot::DA_FLAGS::W_DEPENDENT>();\
            dac->next<ot::DA_FLAGS::W_DEPENDENT>()) {\
          ASSIGN_MAT_PROP_COARSE_TO_FINE_ELEM_BLOCK \
        } /*end dependent loop*/\
        daf->ReadFromGhostElemsBegin<double>(fMatArr,2);\
        /*Overlap communication with INDEPEDENT*/\
        for(dac->init<ot::DA_FLAGS::INDEPENDENT>(),\
            daf->init<ot::DA_FLAGS::WRITABLE>();\
            dac->curr() < dac->end<ot::DA_FLAGS::INDEPENDENT>();\
            dac->next<ot::DA_FLAGS::INDEPENDENT>()) {\
          ASSIGN_MAT_PROP_COARSE_TO_FINE_ELEM_BLOCK \
        } /*end Independent loop */\
        daf->ReadFromGhostElemsEnd<double>(fMatArr);\
        daf->vecRestoreBuffer<double>((*fMatPropVec), fMatArr,\
            true, true, false, 2);\
      } /*end check if active*/\
    }/*end if changedPartition*/\
    dac->vecRestoreBuffer<double>((*cMatPropVec), cMatArr,\
        true, true, true, 2);\
    ctx = new Jac2MFreeData;\
    (static_cast<Jac2MFreeData*>(ctx))->matProp = fMatPropVec;\
    (static_cast<Jac2MFreeData*>(ctx))->isFinestLevel = (i == (nlevels-1));\
    (static_cast<Jac2MFreeData*>(ctx))->Jmat_private = NULL;\
    (static_cast<Jac2MFreeData*>(ctx))->inTmp = NULL;\
    (static_cast<Jac2MFreeData*>(ctx))->outTmp = NULL;\
    damg[i]->user = ctx;\
    if(daf->iAmActive()) {\
      par::Mpi_Reduce<double>(&maxCoeff, &globalMaxCoeff, 1, MPI_MAX, 0, comm);\
      par::Mpi_Reduce<double>(&minCoeff, &globalMinCoeff, 1, MPI_MIN, 0, comm);\
      if(!rank) {\
        std::cout<<"Level: "<<i<<" Max Lap. Coeff: "\
        <<globalMaxCoeff<<" Min Lap. Coeff: "\
        <<globalMinCoeff<<std::endl;\
      }\
      MPI_Barrier(comm);\
    }/*end if active*/\
  }/*end for i*/\
}

void SetUserContextsCoarsestToFinest(ot::DAMG* damg) {
  //jactype =2 only
  int       nlevels = damg[0]->nlevels; //number of multigrid levels

  //coarsest first...
  void * ctx = new Jac2MFreeData;   
  ot::DAMG damg_i = damg[0];

  std::vector<double> * matPropVecPtr = new std::vector<double>;
  (static_cast<Jac2MFreeData*>(ctx))->matProp = matPropVecPtr;
  (static_cast<Jac2MFreeData*>(ctx))->isFinestLevel = (nlevels == 1);
  (static_cast<Jac2MFreeData*>(ctx))->Jmat_private = NULL;
  (static_cast<Jac2MFreeData*>(ctx))->inTmp = NULL;
  (static_cast<Jac2MFreeData*>(ctx))->outTmp = NULL;

  int rank;
  MPI_Comm comm;

  //Default values is the same as the const. coeff. case 
  PetscReal lapFac = 0.0;
  PetscReal massFac = 0.0;
  PetscReal massBase = 1.0;
  PetscTruth optFound;
  PetscOptionsGetReal("lap","-MatPropFac",&lapFac,&optFound);
  PetscOptionsGetReal("mass","-MatPropFac",&massFac,&optFound);
  PetscOptionsGetReal("mass","-BaseMatProp",&massBase,&optFound);
  unsigned int maxD;
  ot::DA* da;
  double *matPropArr = NULL;

  PetscTruth setMatPropFromAnalyticFn;
  PetscOptionsHasName(0,"-setMatPropFromAnalyticFn",&setMatPropFromAnalyticFn);

  if(setMatPropFromAnalyticFn) {
    PetscReal lapFreq = 1.0;
    PetscOptionsGetReal(0,"-lapFreq",&lapFreq,&optFound);
    int dummyInit;
    SET_SINGLE_LEVEL_MAT_PROP_BLOCK(\
        ASSIGN_MAT_PROP_ANALYTIC_FN_ELEM_BLOCK, dummyInit = 0;) 
  } else {
    PetscTruth setCheckerBoardMatProp;
    PetscOptionsHasName(0,"-setCheckerBoardMatProp",&setCheckerBoardMatProp);
    if(setCheckerBoardMatProp) {
      int dummyInit;
      SET_SINGLE_LEVEL_MAT_PROP_BLOCK(ASSIGN_MAT_PROP_CHECKER_BOARD_ELEM_BLOCK, 
          dummyInit = 0;) 
    } else {
      PetscInt numCubes = 1;
      PetscOptionsGetInt(0,"-numCubes",&numCubes,0);
      if(numCubes == 1) {
        int dummyInit;
        SET_SINGLE_LEVEL_MAT_PROP_BLOCK(ASSIGN_MAT_PROP_1_CUBE_ELEM_BLOCK, 
            dummyInit = 0;) 
      } else {
        int dummyInit;
        SET_SINGLE_LEVEL_MAT_PROP_BLOCK(ASSIGN_MAT_PROP_MULTIPLE_CUBES_ELEM_BLOCK,
            dummyInit = 0;) 
      }
    }
  }

  COARSE_TO_FINE_BLOCK

}//end fn.

void SetCoarseToFineFromPts(ot::DAMG* damg,
    const std::vector<double>& pts,
    const std::vector<double> & lapJump) {

  assert(pts.size() == (3*lapJump.size()));

  //jactype =2 only
  int       nlevels = damg[0]->nlevels; //number of multigrid levels

  //coarsest first...
  void * ctx = new Jac2MFreeData;   
  ot::DAMG damg_i = damg[0];

  std::vector<double> * matPropVecPtr = new std::vector<double>;
  (static_cast<Jac2MFreeData*>(ctx))->matProp = matPropVecPtr;
  (static_cast<Jac2MFreeData*>(ctx))->isFinestLevel = (nlevels == 1);
  (static_cast<Jac2MFreeData*>(ctx))->Jmat_private = NULL;
  (static_cast<Jac2MFreeData*>(ctx))->inTmp = NULL;
  (static_cast<Jac2MFreeData*>(ctx))->outTmp = NULL;

  int rank;
  MPI_Comm comm;

  //Default values is the same as the const. coeff. case 
  PetscReal lapBase = 1.0;
  PetscReal massBase = 1.0;
  PetscTruth optFound;
  PetscOptionsGetReal("lap","-BaseMatProp",&lapBase,&optFound);
  PetscOptionsGetReal("mass","-BaseMatProp",&massBase,&optFound);
  unsigned int maxD;
  ot::DA* da;
  double *matPropArr = NULL;

  unsigned int ptsCtr;
  SET_SINGLE_LEVEL_MAT_PROP_BLOCK(ASSIGN_MAT_PROP_FROM_PTS_ELEM_BLOCK,
      ptsCtr = 0;) 

    COARSE_TO_FINE_BLOCK

}//end fn.

#undef ASSIGN_MAT_PROP_FROM_PTS_ELEM_BLOCK 
#undef ASSIGN_MAT_PROP_1_CUBE_ELEM_BLOCK 
#undef ASSIGN_MAT_PROP_MULTIPLE_CUBES_ELEM_BLOCK 
#undef ASSIGN_MAT_PROP_CHECKER_BOARD_ELEM_BLOCK 

#undef SET_SINGLE_LEVEL_MAT_PROP_BLOCK

#undef ASSIGN_MAT_PROP_COARSE_TO_FINE_ELEM_BLOCK 
#undef COARSE_TO_FINE_BLOCK 

#undef ASSIGN_MAT_PROP_ANALYTIC_FN_ELEM_BLOCK 
#undef SQUARE

void DestroyUserContexts(ot::DAMG* damg) {

  PetscInt       jacType = 1;
  PetscOptionsGetInt(0,"-jacType",&jacType,0);

  if(jacType == 1) { return; }

  int       nlevels = damg[0]->nlevels; //number of multigrid levels

  if(jacType == 2) {
    for(int i = 0; i < nlevels; i++) {
      Jac2MFreeData* ctx = (static_cast<Jac2MFreeData*>(damg[i]->user));
      ctx->matProp->clear();
      delete ctx->matProp;
      ctx->matProp = NULL;
      if(ctx->Jmat_private) {
        MatDestroy(ctx->Jmat_private);
        ctx->Jmat_private = NULL;
      }
      if(ctx->inTmp) {
        VecDestroy(ctx->inTmp);
        ctx->inTmp = NULL;
      }
      if(ctx->outTmp) {
        VecDestroy(ctx->outTmp);
        ctx->outTmp = NULL;
      }
      delete ctx;
      ctx = NULL;
    }
  }//end type-2

  if(jacType == 3) {
    for(int i = 0; i < nlevels; i++) {
      Jac3MFreeData* ctx = (static_cast<Jac3MFreeData*>(damg[i]->user));
      //The coarsest level will not create matProp. It will only use the finer
      //level's material properties.
      if(ctx->matProp) {
        ctx->matProp->clear();
        delete ctx->matProp;
        ctx->matProp = NULL;
      }
      //Sometimes, memory is allocated for the fine material property vector.
      //Sometimes it is just a copy of pointers.
      if(ctx->changedPartition) {
        ctx->matPropFine->clear();
        delete ctx->matPropFine;
        ctx->matPropFine = NULL;
      }
      if(ctx->Jmat_private) {
        MatDestroy(ctx->Jmat_private);
        ctx->Jmat_private = NULL;
      }
      if(ctx->inTmp) {
        VecDestroy(ctx->inTmp);
        ctx->inTmp = NULL;
      }
      if(ctx->outTmp) {
        VecDestroy(ctx->outTmp);
        ctx->outTmp = NULL;
      }
      delete ctx;
      ctx = NULL;
    }
  }//end type-3

}//end fn.



