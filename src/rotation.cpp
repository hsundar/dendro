/*
 *@author: Milinda Fernando
 *@date: 09/04/2015 // This is refactored code from HilbertBenchmark code.
 *School of Computing, University of Utah
 * 
 * Contains the code to gernerate the Hilbert Table
 * Later we can hardcorde the table to the dendro header file. 
 * 
 */

#include "../include/rotation.h"

// conntains the rotations permutations for the Hilbert curve. 

std::vector<Rotation3D> rotations_3d;
std::vector<Rotation2D> rotations_2d;



char *HILBERT_2D_TABLE;
char *HILBERT_3D_TABLE;
char * HILBERT_TABLE;


void initializeHilbetTable(int dim)
{
  unsigned int num_children=1<<dim;
  
  if(dim==2)
  {
 
    // NOTE: Assuming the rotations vector is sorted. 
    generateRotationPermutations<Rotation2D>(dim,rotations_2d);
    Rotation2D temp;
    
    HILBERT_TABLE=new char[rotations_2d.size()*num_children];
    
    
//     for(int i=0;i<rotations_2d.size();i++){
//       std::cout<<"Rotation:"<<rotations_2d[i].rot_perm_str()<<"\t Rot_Index:"<<rotations_2d[i].rot_index_str()<<std::endl;
//     }

    for(int i=0;i<rotations_2d.size();i++)
    {
      
      for(int j=0;j<num_children;j++)
      {
	temp=rotations_2d[i];
	rotate(j,temp.rot_perm,temp.rot_index,dim,false);
	
	bool found=false;
	int index=0;
	
	for(int w=0;w<rotations_2d.size();w++)
	{
	  if(rotations_2d[w]==temp){
	    found=true;
	    index=w;
	    break;
	  }
	}
	if(found==false)
	{
	  std::cout<<"Rotation Permutations Error. Found a rotations permutation which is not in the list"<<std::endl;
	}else
	{
	  
	  HILBERT_TABLE[i*num_children+j]=index;
	  //std::cout<<"HILBERT_TABLE["<<(i*num_children+j)<<"]="<<(int)HILBERT_TABLE[i*num_children+j]<<";"<<std::endl;
	}
      }
    }
    
  }else if(dim==3)
  {
    
    // NOTE: Assuming the rotations vector is sorted. 
    generateRotationPermutations<Rotation3D>(dim,rotations_3d);
     for(int i=0;i<rotations_3d.size();i++){
       std::cout<<"Rotation:"<<rotations_3d[i].rot_perm_str()<<"\t Rotation_index:"<<rotations_3d[i].rot_index_str()<<std::endl;
     }
  

    HILBERT_TABLE=new char[rotations_3d.size()*num_children];
    Rotation3D temp;
    
    for(int i=0;i<rotations_3d.size();i++)
    {
     
      for(int j=0;j<num_children;j++)
      { 
	temp=rotations_3d[i];
	rotate(j,temp.rot_perm,temp.rot_index,dim,false);
	bool found=false;
	int index=0;
	for(int w=0;w<rotations_3d.size();w++)
	{
	  if(rotations_3d[w]==temp){
	    found=true;
	    index=w;
	    break;
	  }
	}
	if(found==false)
	{
	  std::cout<<"Rotation Permutations Error. Found a rotations permutation which is not in the list"<<std::endl;
	}else
	{
	  HILBERT_TABLE[i*num_children+j]=index;
	  std::cout<<"HILBERT_TABLE["<<(i*num_children+j)<<"]="<<(int)HILBERT_TABLE[i*num_children+j]<<";"<<std::endl;
	}
	
      }
    }
  }
  
  
  
  
}


void rotate_table_based(int index,int& current_rot,int dim)
{
  int index_temp=0;
  if(dim==2)
  {
    
     Rotation2D temp=rotations_2d[current_rot];
      index_temp=temp.rot_index[index];
     current_rot=HILBERT_TABLE[current_rot*4+index_temp];
  }else if(dim==3)
  {

      Rotation3D temp=rotations_3d[current_rot];
      index_temp=temp.rot_index[index];
      current_rot=HILBERT_TABLE[current_rot*8+index_temp];
  
  }
  
  
}
