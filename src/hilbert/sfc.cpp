#include "../../include/hilbert/sfc.h"

int G_MAX_DEPTH=0;
int G_dim=0;

void rotate(int index,char* current,char * rot_index,int dim,bool cal_true_index)
{
  
  if(dim==2)
  {
    if(cal_true_index){
      index=rot_index[index];
    }
     if(index==0)
    {
      rot_index[current[1]]=3;
      rot_index[current[3]]=1;
      SWAP(current[1],current[3]); // RIGHT Rotate and flip orientation
      
      
    }else if (index==3)
    {
      rot_index[current[0]]=2;
      rot_index[current[2]]=0;
      SWAP(current[0],current[2]); //LEFT Rotate and flip orientation: 
      
    }
    
  }else if(dim==3)
  {	
    if(cal_true_index){
      index=rot_index[index];
    }
    if(index==0)
    {
      rot_index[current[1]]=7;
      rot_index[current[7]]=1;
      SWAP(current[1],current[7]);
      rot_index[current[4]]=2;
      rot_index[current[2]]=4;
      SWAP(current[2],current[4]);
      
    }else if(index==1)
    {
      rot_index[current[3]]=7;
      rot_index[current[7]]=3;
      SWAP(current[3],current[7]);
      rot_index[current[2]]=6;
      rot_index[current[6]]=2;
      SWAP(current[2],current[6]);
      
    }else if(index==3)
    {
      rot_index[current[3]]=5;
      rot_index[current[5]]=3;
      SWAP(current[3],current[5]);
      
      rot_index[current[3]]=7;
      rot_index[current[7]]=3;
      SWAP(current[3],current[7]);
      
      rot_index[current[2]]=6;
      rot_index[current[6]]=2;
      SWAP(current[2],current[6]);
      
      rot_index[current[0]]=2;
      rot_index[current[2]]=0;
      SWAP(current[0],current[2]);
    }else if(index==4)
    {
      rot_index[current[1]]=7;
      rot_index[current[7]]=1;
      SWAP(current[1],current[7]);
      
      rot_index[current[1]]=5;
      rot_index[current[5]]=1;
      SWAP(current[1],current[5]);
      
      rot_index[current[0]]=4;
      rot_index[current[4]]=0;
      SWAP(current[0],current[4]);
      
      rot_index[current[0]]=2;
      rot_index[current[2]]=0;
      SWAP(current[0],current[2]);
      
    }else if(index==6)
    {
      rot_index[current[1]]=5;
      rot_index[current[5]]=1;
      SWAP(current[1],current[5]);
      
      rot_index[current[0]]=4;
      rot_index[current[4]]=0;
      SWAP(current[0],current[4]);
      
    }else if(index==7)
    {
      
      rot_index[current[0]]=6;
      rot_index[current[6]]=0;
      SWAP(current[0],current[6]);
      
      rot_index[current[3]]=5;
      rot_index[current[5]]=3;
      SWAP(current[3],current[5]);
    }
    
  }
  
}

