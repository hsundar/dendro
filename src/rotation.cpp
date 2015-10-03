///*
// *@author: Milinda Fernando
// *@date: 09/04/2015 // This is refactored code from HilbertBenchmark code.
// *School of Computing, University of Utah
// *
// * Contains the code to gernerate the Hilbert Table
// * Later we can hardcorde the table to the dendro header file.
// *
// */
//
//#include "../include/rotation.h"
//
//// conntains the rotations permutations for the Hilbert curve.
//
//std::vector<Rotation3D> rotations;
//
//
//
//
//char *HILBERT_2D_TABLE;
//char *HILBERT_3D_TABLE;
//char * HILBERT_TABLE;
//
//
//void initializeHilbetTable(int dim)
//{
//  unsigned int num_children=1<<dim;
//
//  if(dim==2)
//  {
//
//    // NOTE: Assuming the rotations vector is sorted.
//    generateRotationPermutations<Rotation3D>(dim,rotations);
//    Rotation3D temp;
//
//    HILBERT_TABLE=new char[rotations.size()*num_children];
//
//
////     for(int i=0;i<rotations.size();i++){
////       std::cout<<"Rotation:"<<rotations[i].rot_perm_str()<<"\t Rot_Index:"<<rotations[i].rot_index_str()<<std::endl;
////     }
//
//    for(int i=0;i<rotations.size();i++)
//    {
//
//      for(int j=0;j<num_children;j++)
//      {
//	temp=rotations[i];
//	rotate(j,temp.rot_perm,temp.rot_index,dim,false);
//
//	bool found=false;
//	int index=0;
//
//	for(int w=0;w<rotations.size();w++)
//	{
//	  if(rotations[w]==temp){
//	    found=true;
//	    index=w;
//	    break;
//	  }
//	}
//	if(found==false)
//	{
//	  std::cout<<"Rotation Permutations Error. Found a rotations permutation which is not in the list"<<std::endl;
//	}else
//	{
//
//	  HILBERT_TABLE[i*num_children+j]=index;
//      //std::cout<<"HILBERT_TABLE["<<(i*num_children+j)<<"]="<<(int)HILBERT_TABLE[i*num_children+j]<<";"<<std::endl;
//	}
//      }
//    }
//
//  }else if(dim==3)
//  {
//
//    // NOTE: Assuming the rotations vector is sorted.
//    generateRotationPermutations<Rotation3D>(dim,rotations);
////     for(int i=0;i<rotations.size();i++){
////       std::cout<<"Rotation:"<<rotations[i].rot_perm_str()<<"\t Rotation_index:"<<rotations[i].rot_index_str()<<std::endl;
////     }
//
//
//    HILBERT_TABLE=new char[rotations.size()*num_children];
//    Rotation3D temp;
//
//    for(int i=0;i<rotations.size();i++)
//    {
//
//      for(int j=0;j<num_children;j++)
//      {
//	temp=rotations[i];
//	rotate(j,temp.rot_perm,temp.rot_index,dim,false);
//	bool found=false;
//	int index=0;
//	for(int w=0;w<rotations.size();w++)
//	{
//	  if(rotations[w]==temp){
//	    found=true;
//	    index=w;
//	    break;
//	  }
//	}
//	if(found==false)
//	{
//	  std::cout<<"Rotation Permutations Error. Found a rotations permutation which is not in the list"<<std::endl;
//	}else
//	{
//	  HILBERT_TABLE[i*num_children+j]=index;
//	  //std::cout<<"HILBERT_TABLE["<<(i*num_children+j)<<"]="<<(int)HILBERT_TABLE[i*num_children+j]<<";"<<std::endl;
//	}
//
//      }
//    }
//  }
//
//
//
//
//}
//
//
//void rotate_table_based(int index,int& current_rot,int dim)
//{
//  int index_temp=0;
//  if(dim==2)
//  {
//
//     Rotation3D temp=rotations[current_rot];
//      index_temp=temp.rot_index[index];
//     current_rot=HILBERT_TABLE[current_rot*4+index_temp];
//  }else if(dim==3)
//  {
//
//      Rotation3D temp=rotations[current_rot];
//      index_temp=temp.rot_index[index];
//      current_rot=HILBERT_TABLE[current_rot*8+index_temp];
//
//  }
//
//
//}
