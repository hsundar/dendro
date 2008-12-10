#include <mpi.h>
#include <iostream>
#include <vector>
#include <cassert>
#include <stdlib.h>

int main(int argc, char**argv) {

  int len = atoi(argv[1]);
  
  double stTime = MPI_Wtime();
  std::vector<int> arr1;  
  for(int i = 0; i < len; i++) {
    arr1.push_back(i);
  } 
  double endTime = MPI_Wtime();
  std::cout<<"Initial push_back: "<<(endTime - stTime)<<std::endl;
  arr1.clear();
  
  stTime = MPI_Wtime();
  std::vector<int> arr2(len);
  for(int i = 0; i < len; i++) {
    arr2[i] = i;
  } 
  endTime = MPI_Wtime();
  std::cout<<"Alloc and set: "<<(endTime - stTime)<<std::endl;

  stTime = MPI_Wtime();
  int* arr2ptr = (&(*(arr2.begin())));
  for(int i = 0; i < len; i++) {
    arr2ptr[i] = i;
  } 
  endTime = MPI_Wtime();
  std::cout<<"Set using ptr: "<<(endTime - stTime)<<std::endl;

  stTime = MPI_Wtime();
  std::vector<int> arr3;
  for(int i = 0; i < len; i++) {
    arr3.push_back(arr2[i]);
  }  
  endTime = MPI_Wtime();
  std::cout<<"copy by push_back: "<<(endTime - stTime)<<std::endl;  
  arr2.clear();

  std::vector<int> arr4;
  assert(arr4.begin() == arr4.end());
  stTime = MPI_Wtime();
  arr4.insert(arr4.begin(), arr3.begin(), (arr3.begin() + len));
  endTime = MPI_Wtime();
  std::cout<<"copy by insert: "<<(endTime - stTime) <<std::endl;
  arr3.clear();
  arr4.clear();  
}

