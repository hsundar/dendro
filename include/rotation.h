#ifndef ROTATION_H
#define ROTATION_H

#include "sfc.h"
#include <cstring>
#include <vector>
#include <algorithm>
#include <cstdio>



extern char **HILBERT_2D_TABLE;
extern char **HILBERT_3D_TABLE;



/*
 *
 * Rotations for 2D and 3D case. 
 *
 */

struct Rotation2D
{
  char rot_perm[4];
  char rot_index[4];
  
  Rotation2D()
  {
    for(int i=0;i<4;i++)
    {
      rot_perm[i]=i;
      rot_index[i]=i;
    }
  }
  
  Rotation2D(char* p_rot_perm,char* p_rot_index)
  {
    memcpy(rot_perm,p_rot_perm,sizeof(char)*4);
    memcpy(rot_index,p_rot_index,sizeof(char)*4);
  }
  
  
  Rotation2D& operator =(const Rotation2D& other)
  {
    memcpy(this->rot_perm,other.rot_perm,sizeof(char)*4);
    memcpy(this->rot_index,other.rot_index,sizeof(char)*4);
    return *this;
  }
  
  bool operator==(Rotation2D& other)
  {
    //return (std::string(this->rot_perm)==std::string(other.rot_perm));
    for(int i=0;i<4;i++)
    {
      if(this->rot_perm[i]!=other.rot_perm[i])
	return false;
    }
    return true;
  }
  
  bool operator <(const Rotation2D & other) const
  {
    //return (std::string(this->rot_perm)<std::string(other.rot_perm));
    return (rot_perm_str()<other.rot_perm_str());
    
    
  }
  
  std::string rot_perm_str() const
  {
    char str[4];
    std::sprintf(str,"%d%d%d%d",this->rot_perm[0],this->rot_perm[1],this->rot_perm[2],this->rot_perm[3]);
    return std::string(str);
  }
  
  std::string rot_index_str() const
  {
    char str[4];
    std::sprintf(str,"%d%d%d%d",this->rot_index[0],this->rot_index[1],this->rot_index[2],this->rot_index[3]);
    return std::string(str);
  }
  
  
};


struct Rotation3D
{
  char rot_perm[8];
  char rot_index[8];
  
  Rotation3D()
  {
    for(int i=0;i<8;i++)
    {
      rot_perm[i]=i;
      rot_index[i]=i;
    }
  }
  
  Rotation3D(char* p_rot_perm,char* p_rot_index)
  {
    memcpy(rot_perm,p_rot_perm,sizeof(char)*8);
    memcpy(rot_index,p_rot_index,sizeof(char)*8);
  }
  
  Rotation3D& operator =(const Rotation3D& other)
  {
    memcpy(this->rot_perm,other.rot_perm,sizeof(char)*8);
    memcpy(this->rot_index,other.rot_index,sizeof(char)*8);
    return *this;
    
  }
  
  bool operator==(Rotation3D& other)
  {
    //return (std::string(this->rot_perm)==std::string(other.rot_perm));
    
    for(int i=0;i<8;i++)
    {
      if(this->rot_perm[i]!=other.rot_perm[i])
	return false;
    }
    return true;
  }
  
  bool operator <(const Rotation3D& other) const
  {
   return (rot_perm_str()<other.rot_perm_str());
  }
  
  std::string rot_perm_str() const
  {
    char str[8];
    std::sprintf(str,"%d%d%d%d%d%d%d%d",this->rot_perm[0],this->rot_perm[1],this->rot_perm[2],this->rot_perm[3],this->rot_perm[4],this->rot_perm[5],this->rot_perm[6],this->rot_perm[7]);
    return std::string(str);
  }
  
  std::string rot_index_str() const
  {
    char str[8];
    std::sprintf(str,"%d%d%d%d%d%d%d%d",this->rot_index[0],this->rot_index[1],this->rot_index[2],this->rot_index[3],this->rot_index[4],this->rot_index[5],this->rot_index[6],this->rot_index[7]);
    return std::string(str);
  }
  
  
  
};

extern std::vector<Rotation2D> rotations_2d;
extern std::vector<Rotation3D> rotations_3d;

template<typename T>
void insert_unique(T& item,std::vector<T>& unique_rot_patterns,std::vector<T>& rot_patterns)
{
  for (int i=0;i<unique_rot_patterns.size();i++)
  {
    if(unique_rot_patterns[i]==item){
      return;
    }
  }
  
  rot_patterns.push_back(item);
  unique_rot_patterns.push_back(item);
 
}

template <typename T>
void generateRotationPermutations(int dim,std::vector<T>& rotation_table)
{
  rotation_table.clear();
  std::vector<T> temp_rot;
  if(dim==2)
  {
     int rot_count=0;
     T rot[4];
     char rot_current[4]={0,1,2,3};
     char rot_index[4]={0,1,2,3};
     
     T default_rotation(rot_current,rot_index);
     //Rotation2D current_rot;
     insert_unique<T>(default_rotation,rotation_table,temp_rot);
       
      while(temp_rot.size()!=0)
      {
	
	T rotation=temp_rot[0];
	temp_rot.erase(temp_rot.begin());
			
	for(int i=0;i<4;i++){
	  
	  rot[i]=rotation;//T(rotation.rot_perm,rotation.rot_index);
	  
	  rotate(i,rot[i].rot_perm,rot[i].rot_index,2);
	  insert_unique(rot[i],rotation_table,temp_rot);
	}
      
	
      }
      std::sort(rotation_table.begin(),rotation_table.end());
     
  }else if(dim==3)
  {
      //default rotation pattern inserted.
      
      
      int rot_count=0;
      T rot[8];
      
      char rot_current[8]={0,1,2,3,4,5,6,7};
      char rot_index[8]={0,1,2,3,4,5,6,7};
      
      T default_rotation(rot_current,rot_index);
            
      insert_unique<T>(default_rotation,rotation_table,temp_rot);
      
      while(temp_rot.size()!=0)
      {
	
	T rotation=temp_rot[0];
	temp_rot.erase(temp_rot.begin());
	
	for(int i=0;i<8;i++){
	
	  rot[i]=rotation;//(rotation.rot_perm,rotation.rot_index);
	  
	  rotate(i,rot[i].rot_perm,rot[i].rot_index,3);
	  insert_unique(rot[i],rotation_table,temp_rot);
	}
      	
      }
     
      std::sort(rotation_table.begin(),rotation_table.end());
  }
}

void initializeHilbetTable(int dim);

void rotate_table_based(int index,int& current_rot,int dim);

#endif