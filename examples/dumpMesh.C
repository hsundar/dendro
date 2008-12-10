
#include "mpi.h"
#include "petsc.h"
#include "sys.h"
#include "parUtils.h"
#include "TreeNode.h"
#include "oda.h"
#include "colors.h"
#include "externVars.h"
#include "dendro.h"

//Don't want time to be synchronized. Need to check load imbalance.
#ifdef MPI_WTIME_IS_GLOBAL
#undef MPI_WTIME_IS_GLOBAL
#endif

int main(int argc, char ** argv ) {	
  int size, rank;
  char bFile[50];
  DendroIntL locSz, totalSz;

  PetscInitialize(&argc,&argv,0,NULL);
  ot::RegisterEvents();

  if(argc < 2) {
    std::cerr << "Usage: " << argv[0] << "inpfile " << std::endl;
    return -1;
  }

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  sprintf(bFile,"%s%d_%d.ot",argv[1],rank,size);

  if(!rank){
    std::cout << " reading  "<<bFile<<std::endl; // Point size
  }
  std::vector<ot::TreeNode> balOct;
  ot::readNodesFromFile (bFile,balOct);

  MPI_Barrier(MPI_COMM_WORLD);

  if(!rank){
    std::cout << " finished reading  "<<bFile<<std::endl; // Point size
  }

  // compute total inp size and output size
  locSz = balOct.size();
  par::Mpi_Reduce<DendroIntL>(&locSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  if(!rank) {
    std::cout << "# of Balanced Octants: "<< totalSz << std::endl;       
  }

  //ODA ...
  MPI_Barrier(MPI_COMM_WORLD);	

  assert(!(balOct.empty()));
  ot::DA da(balOct,MPI_COMM_WORLD, MPI_COMM_WORLD, false,NULL,NULL);

  balOct.clear();

  MPI_Barrier(MPI_COMM_WORLD);

  if(!rank) {
    std::cout << "Finished Meshing "<< std::endl;       
  }

  MPI_Barrier(MPI_COMM_WORLD);

  unsigned int maxDepth = (da.getMaxDepth()-1);

  //ComputeLocalToGlobalElemMappings
  if( !(da.computedLocalToGlobalElems()) ) {
    da.computeLocalToGlobalElemMappings();
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if(!rank) {
    std::cout << "Computed Local To Global (Nodes) Maps "<< std::endl;       
  }

  MPI_Barrier(MPI_COMM_WORLD);

  //ComputeLocalToGlobalNodeMappings
  if( !(da.computedLocalToGlobal()) ) {
    da.computeLocalToGlobalMappings();
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if(!rank) {
    std::cout << "Computed Local To Global (Elements) Maps "<< std::endl;       
  }

  if(!rank) {
    std::cout << "Dumping Mesh... "<< std::endl;       
  }

  DendroIntL* localToGlobalMap = da.getLocalToGlobalMap();
  DendroIntL* localToGlobalElemsMap = da.getLocalToGlobalElemsMap();

  unsigned int localElemSize = da.getElementSize();
  unsigned int localNodeSize = da.getNodeSize();

  std::vector<char> nodeVisited;
  da.createVector<char>(nodeVisited, false, false, 1);

  //Initialize
  for(unsigned int i = 0; i < nodeVisited.size(); i++) {
    nodeVisited[i] = 0;
  }

  char* nodeVisitedBuffer = NULL;
  da.vecGetBuffer<char>(nodeVisited, nodeVisitedBuffer, false, false, true, 1);

  unsigned int idxElemBegin = da.getIdxElementBegin();
  unsigned int idxPostGhostBegin = da.getIdxPostGhostBegin();

  char nodeListFileName[100];
  sprintf(nodeListFileName,"%s_NL_%d_%d.out",argv[1],rank,size);
  FILE* outNodeListFile = fopen(nodeListFileName,"wb");
  fwrite((&localNodeSize), sizeof(unsigned int), 1, outNodeListFile);

  unsigned int nodeCnt = 0;
  if(da.iAmActive()) {
    for( da.init<ot::DA_FLAGS::ALL>(); 
        da.curr() < da.end<ot::DA_FLAGS::ALL>();
        da.next<ot::DA_FLAGS::ALL>() ) {
      Point currPt = da.getCurrentOffset();
      unsigned int Lev = da.getLevel(da.curr());
      unsigned int xCurr = currPt.xint();
      unsigned int yCurr = currPt.yint();
      unsigned int zCurr = currPt.zint();

      ot::TreeNode current(xCurr, yCurr, zCurr, Lev, 3, (maxDepth+1));
      ot::TreeNode parent = current.getParent();

      unsigned char hnType = da.getHangingNodeIndex(da.curr());
      unsigned int xyz[8][3];

      //node 0
      if(hnType & (1 << 0)) {
        xyz[0][0] = parent.minX();
        xyz[0][1] = parent.minY();
        xyz[0][2] = parent.minZ();
      } else {
        xyz[0][0] = current.minX();
        xyz[0][1] = current.minY();
        xyz[0][2] = current.minZ();
      }

      //node 1
      if(hnType & (1 << 1)) {
        xyz[1][0] = parent.maxX();
        xyz[1][1] = parent.minY();
        xyz[1][2] = parent.minZ();
      } else {
        xyz[1][0] = current.maxX();
        xyz[1][1] = current.minY();
        xyz[1][2] = current.minZ();
      }

      //node 2
      if(hnType & (1 << 2)) {
        xyz[2][0] = parent.minX();
        xyz[2][1] = parent.maxY();
        xyz[2][2] = parent.minZ();
      } else {
        xyz[2][0] = current.minX();
        xyz[2][1] = current.maxY();
        xyz[2][2] = current.minZ();
      }

      //node 3
      if(hnType & (1 << 3)) {
        xyz[3][0] = parent.maxX();
        xyz[3][1] = parent.maxY();
        xyz[3][2] = parent.minZ();
      } else {
        xyz[3][0] = current.maxX();
        xyz[3][1] = current.maxY();
        xyz[3][2] = current.minZ();
      }

      //node 4
      if(hnType & (1 << 4)) {
        xyz[4][0] = parent.minX();
        xyz[4][1] = parent.minY();
        xyz[4][2] = parent.maxZ();
      } else {
        xyz[4][0] = current.minX();
        xyz[4][1] = current.minY();
        xyz[4][2] = current.maxZ();
      }

      //node 5
      if(hnType & (1 << 5)) {
        xyz[5][0] = parent.maxX();
        xyz[5][1] = parent.minY();
        xyz[5][2] = parent.maxZ();
      } else {
        xyz[5][0] = current.maxX();
        xyz[5][1] = current.minY();
        xyz[5][2] = current.maxZ();
      }

      //node 6
      if(hnType & (1 << 6)) {
        xyz[6][0] = parent.minX();
        xyz[6][1] = parent.maxY();
        xyz[6][2] = parent.maxZ();
      } else {
        xyz[6][0] = current.minX();
        xyz[6][1] = current.maxY();
        xyz[6][2] = current.maxZ();
      }

      //node 7
      if(hnType & (1 << 7)) {
        xyz[7][0] = parent.maxX();
        xyz[7][1] = parent.maxY();
        xyz[7][2] = parent.maxZ();
      } else {
        xyz[7][0] = current.maxX();
        xyz[7][1] = current.maxY();
        xyz[7][2] = current.maxZ();
      }

      unsigned int indices[8];
      da.getNodeIndices(indices);

      for(unsigned int j = 0; j < 8; j++) {
        if( (indices[j] >= idxElemBegin) && (indices[j] < idxPostGhostBegin) &&
            (nodeVisitedBuffer[indices[j]] == 0) ) {
          fwrite(xyz[j], sizeof(unsigned int), 3,  outNodeListFile);
          fwrite((&(localToGlobalMap[indices[j]])), sizeof(DendroIntL), 1, outNodeListFile);
          nodeVisitedBuffer[indices[j]] = 1;
          nodeCnt++;
        }
      }//end for j
    }//end loop
  }//end if active

  da.vecRestoreBuffer<char>(nodeVisited, nodeVisitedBuffer, false, false, true, 1);
  nodeVisited.clear();

  assert(nodeCnt == localNodeSize);
  fclose(outNodeListFile);

  MPI_Barrier(MPI_COMM_WORLD);
  if(!rank) {
    std::cout<<"Finished Writing Node List."<<std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);


  char elemListFileName[100];
  sprintf(elemListFileName,"%s_EL_%d_%d.out",argv[1],rank,size);
  FILE* outElemListFile = fopen(elemListFileName,"wb");
  fwrite(&maxDepth, sizeof(unsigned int), 1, outElemListFile);
  fwrite(&localElemSize, sizeof(unsigned int), 1, outElemListFile);

  unsigned int writableCnt = 0;
  if(da.iAmActive()) {
    for( da.init<ot::DA_FLAGS::WRITABLE>(); da.curr() < da.end<ot::DA_FLAGS::WRITABLE>(); da.next<ot::DA_FLAGS::WRITABLE>() ) {
      Point currPt = da.getCurrentOffset();
      DendroIntL xyzdi[5];
      xyzdi[0] = currPt.xint();
      xyzdi[1] = currPt.yint();
      xyzdi[2] = currPt.zint();
      xyzdi[3] = ( (da.getLevel(da.curr())) - 1 );
      xyzdi[4] = localToGlobalElemsMap[da.curr()];

      fwrite(xyzdi, sizeof(DendroIntL), 5, outElemListFile);

      writableCnt++; 
    }
  }

  assert(writableCnt == localElemSize); 
  fclose(outElemListFile);

  MPI_Barrier(MPI_COMM_WORLD);
  if(!rank) {
    std::cout<<"Finished Writing Element List."<<std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  char meshFileName[100];
  sprintf(meshFileName,"%s_M_%d_%d.out",argv[1],rank,size);
  FILE* outMeshFile = fopen(meshFileName,"wb");
  fwrite(&localElemSize, sizeof(unsigned int), 1, outMeshFile);

  if(da.iAmActive()) {
    for( da.init<ot::DA_FLAGS::WRITABLE>();
        da.curr() < da.end<ot::DA_FLAGS::WRITABLE>();
        da.next<ot::DA_FLAGS::WRITABLE>() ) {
      unsigned int indices[8];
      da.getNodeIndices(indices);

      DendroIntL record[9];
      record[0] = localToGlobalElemsMap[da.curr()];

      for(unsigned int j = 0; j < 8; j++) {
        record[j+1] = localToGlobalMap[indices[j]];
      }

      fwrite(record, sizeof(DendroIntL), 9, outMeshFile);

      writableCnt++; 
    }
  }

  fclose(outMeshFile);

  MPI_Barrier(MPI_COMM_WORLD);
  if(!rank) {
    std::cout<<"Finished Dumping Mesh."<<std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (!rank) {
    std::cout << GRN << "Finalizing PETSC" << NRM << std::endl;
  }

  PetscFinalize();
}//end function

