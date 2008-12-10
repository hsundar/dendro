#include <vector>
#include <cstdlib>
#include <time.h>
#include <iostream>

/*
Compare the speeds for using push_back, insertion into a vector and insertion into a array.
Input Argument: size of the list.
Author: Rahul S. Sampath
*/

int main(int argc, char ** argv) {
  time_t secs1, secs2;
  int n = atoi(argv[1]);
  std::vector<double> vec;

  vec.resize(n);

  secs1 = time(NULL);
  for(int i=0;i<n;i++) {
    vec[i] = (2.0);
  }
  secs2 = time(NULL);

  std::cout<<"Inserting "<<n<<" elements in a vector takes "<<(secs2-secs1)<<" seconds."<<std::endl;

  vec.clear();

  secs1 = time(NULL);
  for(int i=0;i<n;i++) {
    vec.push_back(2.0);
  }
  secs2 = time(NULL);

  vec.clear();

  std::cout<<"Using push_back for "<<n<<" elements takes "<<(secs2-secs1)<<" seconds."<<std::endl;

  secs1 = time(NULL);
  double * arr = new double[n];
  secs2 = time(NULL);
  std::cout<<"Malloc for "<<n<<" elements takes "<<(secs2-secs1)<<" seconds."<<std::endl;

  secs1 = time(NULL);
  for(int i=0;i<n;i++) {
    arr[i] = (2.0);
  }
  secs2 = time(NULL);

  std::cout<<"Inserting "<<n<<" elements in an array takes "<<(secs2-secs1)<<" seconds."<<std::endl;

  secs1 = time(NULL);
  delete [] arr;
  secs2 = time(NULL);
  std::cout<<"Freeing memory for "<<n<<" elements takes "<<(secs2-secs1)<<" seconds."<<std::endl;

  //Iterative variable length insertions...
  secs1 = time(NULL);
   for(int i=0;i<n;i++) {
     std::vector<double> tmp;
     for(int j=(i*i);j<((i+1)*(i+1));j++) {
       tmp.push_back(2.0);
     }
     vec.push_back(tmp.begin(),tmp.end());
     tmp.clear();
   }
  secs2 = time(NULL);

  std::cout<<"Using push_back (for varaible lenghts) "<<vec.size()
    <<" elements takes "<<(secs2-secs1)<<" seconds."<<std::endl;
  vec.clear();
}
