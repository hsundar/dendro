
#include "mpi.h"
#include "nodeAndValues.h"
#include "TreeNode.h"
#include "parUtils.h"
#include "externVars.h"

int main(int argc, char** argv) {

  MPI_Init(&argc, &argv);

  ot::NodeAndValues<double, 3> obj;
  ot::NodeAndValues<double, 3> trueObj;

  //Default
  obj.node = ot::TreeNode(3, 30);
  obj.values[0] = 0.0;
  obj.values[1] = 0.0;
  obj.values[2] = 0.0;

  int rank, npes;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  ot::TreeNode myNode(3, 30);

  for(int i = 0; i < 5; i++) {
    std::vector<ot::TreeNode> octs;
    myNode.addChildren(octs);
    myNode = octs[1];
  }

  trueObj.node = myNode;
  trueObj.values[0] = 10.5;
  trueObj.values[1] = 220.2;
  trueObj.values[2] = 1324.111;

  assert(obj != trueObj);
  assert(!obj.Equals(trueObj));
  assert(obj < trueObj);
  assert(trueObj > obj);

  if(!rank) {
    obj = trueObj;
  }
  
  par::Mpi_Bcast<ot::NodeAndValues<double, 3> >(&obj, 1, 0, MPI_COMM_WORLD);

  assert(obj == trueObj);
  assert(obj.Equals(trueObj));

  MPI_Barrier(MPI_COMM_WORLD);

  std::cout<<"Success!"<<std::endl;

  MPI_Finalize();
}


