/*
 *@author: Milinda Fernando
 *@date: 09/04/2015 // This is refactored code from HilbertBenchmark code.
 *School of Computing, University of Utah
 * 
 * Contains the Ordering Implementations of Hilbert Recursive method and NCA based Method. 
 * 
 */
/*
 * Recursive Hilbert order implementation. 2D and 3D
 */
#include "hilbert.h"

int findIndex(Point * pt, int x, int y, int z,int len)
{
    for(int i=0;i<len;i++)
    {
      if(pt[i].xint()==x && pt[i].yint()==y && pt[i].zint()==z)
      {
	return i;
      }
    } 
  
}

bool hilbert_order(const Point & p1,const Point& p2)
{
  
 
  int x1=p1.xint();
  int y1=p1.yint();
  int z1=p1.zint();
  
  int x2=p2.xint();
  int y2=p2.yint();
  int z2=p2.zint();
  
  if(x1==x2 && y1==y2 && z1==z2)
  {
    return false;
  }
  
  int index1=0;
  int index2=0;
  int min_x,min_y,min_z,max_x,max_y,max_z;
  int MAX_LIMIT=1<<G_MAX_DEPTH;
  int len=MAX_LIMIT;
  int deapth=0;
  min_x=0;
  min_y=0;
  min_z=0;
   
  max_x=len;
  max_y=len;
  max_z=len;
  
  unsigned int dim;

  if(G_dim==2)
  {
    Point pt_hilbert[4];
    Point pt_hilbert_new[4];
  
    pt_hilbert[0]=Point ((int)min_x,(int)min_y,(int)0);
    pt_hilbert[1]=Point ((int)min_x,(int)max_y,(int)0);
    pt_hilbert[2]=Point ((int)max_x,(int)max_y,(int)0);
    pt_hilbert[3]=Point ((int)max_x,(int)min_y,(int)0); 
    
    
  
    while (len>1 && deapth<G_MAX_DEPTH)
    {
      
      int xl=pt_hilbert[0].xint();
      int yl=pt_hilbert[0].yint();
      int nca_index=0;
      index1=0;
      index2=0;
      for (int i=1;i<4;i++)
      {
	
	if(xl>pt_hilbert[i].xint())
	{
	  xl=pt_hilbert[i].xint();
	  yl=pt_hilbert[i].yint();
	  nca_index=i;
	}else if(xl==pt_hilbert[i].xint())
	{
	  if(yl>pt_hilbert[i].yint())
	  {
	    yl=pt_hilbert[i].yint();
	    nca_index=i;
	  }
	}
	 
	  
      }
      int len_nca=len/2;
      // Checking the membership cell.
      if((x1-xl)<len_nca && (y1-yl)<len_nca)
      {	
	index1=findIndex(pt_hilbert,xl,yl,0,4);
      }else if((x1-xl)<len_nca && (y1-yl)>=len_nca)
      {
	index1=findIndex(pt_hilbert,xl,yl+len,0,4);;
      }else if ((x1-xl)>=len_nca && (y1-yl)>=len_nca)
      {
	index1=findIndex(pt_hilbert,xl+len,yl+len,0,4);;
      }else if ((x1-xl)>=len_nca && (y1-yl)<len_nca)
      {
	index1=findIndex(pt_hilbert,xl+len,yl,0,4);;
      }
      
      if((x2-xl)<len_nca && (y2-yl)<len_nca)
      {
	index2=findIndex(pt_hilbert,xl,yl,0,4);
      }else if((x2-xl)<len_nca && (y2-yl)>=len_nca)
      {
	index2=findIndex(pt_hilbert,xl,yl+len,0,4);
      }else if ((x2-xl)>=len_nca && (y2-yl)>=len_nca)
      {
	index2=findIndex(pt_hilbert,xl+len,yl+len,0,4);
      }else if ((x2-xl)>=len_nca && (y2-yl)<len_nca)
      {
	index2=findIndex(pt_hilbert,xl+len,yl,0,4);
      }
               
    if (index1<index2){
	  return true;
    }
    else if (index1>index2){
	return false;
    }
    
   switch (index1)
   {
    case 0:
          pt_hilbert_new[0]=pt_hilbert[0];
	  pt_hilbert_new[1]=(pt_hilbert[0]+pt_hilbert[3])/2;
	  pt_hilbert_new[2]=(pt_hilbert[0]+pt_hilbert[2])/2;
	  pt_hilbert_new[3]=(pt_hilbert[0]+pt_hilbert[1])/2;
	   break;
   case 1:
	  pt_hilbert_new[0]=(pt_hilbert[0]+pt_hilbert[1])/2;
	  pt_hilbert_new[1]=pt_hilbert[1];
	  pt_hilbert_new[2]=(pt_hilbert[1]+pt_hilbert[2])/2;
	  pt_hilbert_new[3]=(pt_hilbert[1]+pt_hilbert[3])/2;
          break;
   case 2:
	  pt_hilbert_new[0]=(pt_hilbert[0]+pt_hilbert[2])/2;
	  pt_hilbert_new[1]=(pt_hilbert[1]+pt_hilbert[2])/2;
	  pt_hilbert_new[2]=pt_hilbert[2];
	  pt_hilbert_new[3]=(pt_hilbert[2]+pt_hilbert[3])/2;
	  break;
   case 3:
          pt_hilbert_new[0]=(pt_hilbert[3]+pt_hilbert[2])/2;
	  pt_hilbert_new[1]=(pt_hilbert[1]+pt_hilbert[3])/2;
	  pt_hilbert_new[2]=(pt_hilbert[3]+pt_hilbert[0])/2;
	  pt_hilbert_new[3]=pt_hilbert[3];
	  break;
   default:
	  std::cout<<"Hilbert ordering error:Invalid nearest cell"<<std::endl;
	  break;
	    
    }
    
    
    pt_hilbert[0]=pt_hilbert_new[0];
    pt_hilbert[1]=pt_hilbert_new[1];
    pt_hilbert[2]=pt_hilbert_new[2];
    pt_hilbert[3]=pt_hilbert_new[3];
    len=len/2;
    deapth+=1;
    
   }
   return false;
  }
  else if(G_dim==3)
  {
    Point pt_hilbert[8];
    Point pt_hilbert_new[8];
    
    pt_hilbert[0]=Point ((int)min_x,(int)min_y,(int)min_z);
    pt_hilbert[1]=Point ((int)min_x,(int)max_y,(int)min_z);
    pt_hilbert[2]=Point ((int)max_x,(int)max_y,(int)min_z);
    pt_hilbert[3]=Point ((int)max_x,(int)min_y,(int)min_z);
      
    pt_hilbert[4]=Point ((int)max_x,(int)min_y,(int)max_z);
    pt_hilbert[5]=Point ((int)max_x,(int)max_y,(int)max_z);
    pt_hilbert[6]=Point ((int)min_x,(int)max_y,(int)max_z);
    pt_hilbert[7]=Point ((int)min_x,(int)min_y,(int)max_z);
  
      
      while(len>1 && deapth<G_MAX_DEPTH)
      {
	
	
      int xl=pt_hilbert[0].xint();
      int yl=pt_hilbert[0].yint();
      int zl=pt_hilbert[0].zint();
      int nca_index=0;
      index1=0;
      index2=0;
      for (int i=1;i<8;i++)
      {
	
	if(xl>pt_hilbert[i].xint())
	{
	  xl=pt_hilbert[i].xint();
	  yl=pt_hilbert[i].yint();
	  nca_index=i;
	}else if(xl==pt_hilbert[i].xint())
	{
	  if(yl>pt_hilbert[i].yint())
	  {
	    yl=pt_hilbert[i].yint();
	    nca_index=i;
	  }else if(yl==pt_hilbert[i].yint())
	  {
	    if(zl>pt_hilbert[i].zint())
	    {
	      zl=pt_hilbert[i].zint();
	      nca_index=i;
	    }
	  }
	  
	}
	 
	  
      }
      int len_nca=len/2;
           
      if((x1-xl)<len_nca && (y1-yl)<len_nca && (z1-zl)<len_nca)
      {
	index1=findIndex(pt_hilbert,xl,yl,zl,8);
	
      }else if((x1-xl)<len_nca && (y1-yl)>=len_nca && (z1-zl)<len_nca)
      {
	index1=findIndex(pt_hilbert,xl,yl+len,zl,8);
      }else if((x1-xl)>=len_nca && (y1-yl)>=len_nca && (z1-zl)<len_nca)
      {
	index1=findIndex(pt_hilbert,xl+len,yl+len,zl,8);
      }else if((x1-xl)>=len_nca && (y1-yl)<len_nca && (z1-zl)<len_nca)
      {
	index1=findIndex(pt_hilbert,xl+len,yl,zl,8);
      }else if ((x1-xl)>=len_nca && (y1-yl)<len_nca && (z1-zl)>=len_nca)
      {
	index1=findIndex(pt_hilbert,xl+len,yl,zl+len,8);
      }else if ((x1-xl)>=len_nca && (y1-yl)>=len_nca && (z1-zl)>=len_nca)
      {
	index1=findIndex(pt_hilbert,xl+len,yl+len,zl+len,8);
      }else if ((x1-xl)<len_nca && (y1-yl)>=len_nca && (z1-zl)>=len_nca)
      {
	index1=findIndex(pt_hilbert,xl,yl+len,zl+len,8);
      }else if ((x1-xl)<len_nca && (y1-yl)<len_nca && (z1-zl)>=len_nca)
      {
	index1=findIndex(pt_hilbert,xl,yl,zl+len,8);
      }
      
       if((x2-xl)<len_nca && (y2-yl)<len_nca && (z2-zl)<len_nca)
      {
	index2=findIndex(pt_hilbert,xl,yl,zl,8);
	
      }else if((x2-xl)<len_nca && (y2-yl)>=len_nca && (z2-zl)<len_nca)
      {
	index2=findIndex(pt_hilbert,xl,yl+len,zl,8);
      }else if((x2-xl)>=len_nca && (y2-yl)>=len_nca && (z2-zl)<len_nca)
      {
	index2=findIndex(pt_hilbert,xl+len,yl+len,zl,8);
      }else if((x2-xl)>=len_nca && (y2-yl)<len_nca && (z2-zl)<len_nca)
      {
	index2=findIndex(pt_hilbert,xl+len,yl,zl,8);
      }else if ((x2-xl)>=len_nca && (y2-yl)<len_nca && (z2-zl)>=len_nca)
      {
	index2=findIndex(pt_hilbert,xl+len,yl,zl+len,8);
      }else if ((x2-xl)>=len_nca && (y2-yl)>=len_nca && (z2-zl)>=len_nca)
      {
	index2=findIndex(pt_hilbert,xl+len,yl+len,zl+len,8);
      }else if ((x2-xl)<len_nca && (y2-yl)>=len_nca && (z2-zl)>=len_nca)
      {
	index2=findIndex(pt_hilbert,xl,yl+len,zl+len,8);
      }else if ((x2-xl)<len_nca && (y2-yl)<len_nca && (z2-zl)>=len_nca)
      {
	index2=findIndex(pt_hilbert,xl,yl,zl+len,8);
      }
	
	
	if (index1<index2){
	  return true;
	}
	else if (index1>index2){
	  return false;
	}
	// means that index1 ==index2
	switch(index1)
	{
	      case 0:
		      pt_hilbert_new[0]=pt_hilbert[0];
		      pt_hilbert_new[1]=(pt_hilbert[0]+pt_hilbert[7])/2;
		      pt_hilbert_new[2]=(pt_hilbert[0]+pt_hilbert[4])/2;
		      pt_hilbert_new[3]=(pt_hilbert[0]+pt_hilbert[3])/2;
		      
		      pt_hilbert_new[4]=(pt_hilbert[0]+pt_hilbert[2])/2;
		      pt_hilbert_new[5]=(pt_hilbert[0]+pt_hilbert[5])/2;
		      pt_hilbert_new[6]=(pt_hilbert[0]+pt_hilbert[6])/2;
		      pt_hilbert_new[7]=(pt_hilbert[0]+pt_hilbert[1])/2;
		      break;
	      case 1:          
		      pt_hilbert_new[0]=(pt_hilbert[0]+pt_hilbert[1])/2;
		      pt_hilbert_new[1]=pt_hilbert[1];
		      pt_hilbert_new[2]=(pt_hilbert[1]+pt_hilbert[6])/2;
		      pt_hilbert_new[3]=(pt_hilbert[1]+pt_hilbert[7])/2;
		      
		      pt_hilbert_new[4]=(pt_hilbert[1]+pt_hilbert[4])/2;
		      pt_hilbert_new[5]=(pt_hilbert[1]+pt_hilbert[5])/2;
		      pt_hilbert_new[6]=(pt_hilbert[1]+pt_hilbert[2])/2;
		      pt_hilbert_new[7]=(pt_hilbert[1]+pt_hilbert[3])/2;
		      break;
	      case 2:
		      pt_hilbert_new[0]=(pt_hilbert[2]+pt_hilbert[0])/2;
		      pt_hilbert_new[1]=(pt_hilbert[2]+pt_hilbert[1])/2;
		      pt_hilbert_new[2]=pt_hilbert[2];
		      pt_hilbert_new[3]=(pt_hilbert[2]+pt_hilbert[3])/2;
		      
		      pt_hilbert_new[4]=(pt_hilbert[2]+pt_hilbert[4])/2;
		      pt_hilbert_new[5]=(pt_hilbert[2]+pt_hilbert[5])/2;
		      pt_hilbert_new[6]=(pt_hilbert[2]+pt_hilbert[6])/2;
		      pt_hilbert_new[7]=(pt_hilbert[2]+pt_hilbert[7])/2;
		      
		      break;
	      case 3:
		      pt_hilbert_new[0]=(pt_hilbert[3]+pt_hilbert[6])/2;
		      pt_hilbert_new[1]=(pt_hilbert[3]+pt_hilbert[1])/2;
		      pt_hilbert_new[2]=(pt_hilbert[3]+pt_hilbert[0])/2;
		      pt_hilbert_new[3]=(pt_hilbert[3]+pt_hilbert[7])/2;
		      
		      pt_hilbert_new[4]=(pt_hilbert[3]+pt_hilbert[4])/2;
		      pt_hilbert_new[5]=pt_hilbert[3];
		      pt_hilbert_new[6]=(pt_hilbert[3]+pt_hilbert[2])/2;
		      pt_hilbert_new[7]=(pt_hilbert[3]+pt_hilbert[5])/2;
		      break;
		      
	      case 4: pt_hilbert_new[0]=(pt_hilbert[4]+pt_hilbert[2])/2;
		      pt_hilbert_new[1]=(pt_hilbert[4]+pt_hilbert[5])/2;
		      pt_hilbert_new[2]=pt_hilbert[4];
		      pt_hilbert_new[3]=(pt_hilbert[4]+pt_hilbert[3])/2;
		      
		      pt_hilbert_new[4]=(pt_hilbert[4]+pt_hilbert[0])/2;
		      pt_hilbert_new[5]=(pt_hilbert[4]+pt_hilbert[7])/2;
		      pt_hilbert_new[6]=(pt_hilbert[4]+pt_hilbert[6])/2;
		      pt_hilbert_new[7]=(pt_hilbert[4]+pt_hilbert[1])/2;
		      break;
	      case 5: 
		      pt_hilbert_new[0]=(pt_hilbert[5]+pt_hilbert[0])/2;
		      pt_hilbert_new[1]=(pt_hilbert[5]+pt_hilbert[1])/2;
		      pt_hilbert_new[2]=(pt_hilbert[5]+pt_hilbert[2])/2;
		      pt_hilbert_new[3]=(pt_hilbert[5]+pt_hilbert[3])/2;
		      
		      pt_hilbert_new[4]=(pt_hilbert[5]+pt_hilbert[4])/2;
		      pt_hilbert_new[5]=pt_hilbert[5];
		      pt_hilbert_new[6]=(pt_hilbert[5]+pt_hilbert[6])/2;
		      pt_hilbert_new[7]=(pt_hilbert[5]+pt_hilbert[7])/2;
		      break;
	      case 6: pt_hilbert_new[0]=(pt_hilbert[6]+pt_hilbert[4])/2;
		      pt_hilbert_new[1]=(pt_hilbert[6]+pt_hilbert[5])/2;
		      pt_hilbert_new[2]=(pt_hilbert[6]+pt_hilbert[2])/2;
		      pt_hilbert_new[3]=(pt_hilbert[6]+pt_hilbert[3])/2;
		      
		      pt_hilbert_new[4]=(pt_hilbert[6]+pt_hilbert[0])/2;
		      pt_hilbert_new[5]=(pt_hilbert[6]+pt_hilbert[1])/2;
		      pt_hilbert_new[6]=pt_hilbert[6];
		      pt_hilbert_new[7]=(pt_hilbert[6]+pt_hilbert[7])/2;
		      break;
	      case 7: pt_hilbert_new[0]=(pt_hilbert[7]+pt_hilbert[6])/2;
		      pt_hilbert_new[1]=(pt_hilbert[7]+pt_hilbert[1])/2;
		      pt_hilbert_new[2]=(pt_hilbert[7]+pt_hilbert[2])/2;
		      pt_hilbert_new[3]=(pt_hilbert[7]+pt_hilbert[5])/2;
		      
		      pt_hilbert_new[4]=(pt_hilbert[7]+pt_hilbert[4])/2;
		      pt_hilbert_new[5]=(pt_hilbert[7]+pt_hilbert[3])/2;
		      pt_hilbert_new[6]=(pt_hilbert[7]+pt_hilbert[0])/2;
		      pt_hilbert_new[7]=pt_hilbert[7];
		      break;
	      default:
		      std::cout<<"Hilbert ordering error:Invalid nearest cell"<<std::endl;
		      break;
	  
	  
	 }
	    
	    
	pt_hilbert[0]=pt_hilbert_new[0];
	pt_hilbert[1]=pt_hilbert_new[1];
	pt_hilbert[2]=pt_hilbert_new[2];
	pt_hilbert[3]=pt_hilbert_new[3];
	
	pt_hilbert[4]=pt_hilbert_new[4];
	pt_hilbert[5]=pt_hilbert_new[5];
	pt_hilbert[6]=pt_hilbert_new[6];
	pt_hilbert[7]=pt_hilbert_new[7];
	
	len=len/2;
	deapth+=1;
    
      }
      return false;
  
     
    
  }
     
}



bool hilbert_order_NCA(const Point& p1,const Point& p2)
{
  
  unsigned int x1 = p1.xint();
  unsigned int x2 = p2.xint();
  
  unsigned int y1 = p1.yint();
  unsigned int y2 = p2.yint();
  
  unsigned int z1 = p1.zint();
  unsigned int z2 = p2.zint();
    
  if(x1==x2 && y1==y2 && z1==z2)
  {
	return false;
  }
 
  unsigned int maxDepth = G_MAX_DEPTH;
  unsigned int maxDiff = (unsigned int)(std::max((std::max((x1^x2),(y1^y2))),(z1^z2)));
  int dim=G_dim;
  
  unsigned int maxDiffBinLen = binOp::binLength(maxDiff);
  //Eliminate the last maxDiffBinLen bits.
  unsigned int ncaX = ((x1>>maxDiffBinLen)<<maxDiffBinLen);
  unsigned int ncaY = ((y1>>maxDiffBinLen)<<maxDiffBinLen);
  unsigned int ncaZ = ((z1>>maxDiffBinLen)<<maxDiffBinLen);
  unsigned int ncaLev = (maxDepth - maxDiffBinLen);
  
  //std:cout<<"nca_x:"<<ncaX<<" nca_y:"<<ncaY<<" nca_z:"<<ncaZ<<"NCA Level:"<<ncaLev<<std::endl;
  
  unsigned int xl=0;
  unsigned int yl=0;
  unsigned int zl=0;
  
  unsigned int len=1<<G_MAX_DEPTH;
  int count=0;
  unsigned int index1=0;
  unsigned int index2=0;
  
  //Rotation2D temp_2d;
  //Rotation3D temp_3d;
  int index_temp=0;
  
  if(G_dim==2)
  {
     char rotation[4]={0,1,2,3};
     char rot_index[4]={0,1,2,3};
     int current_rot=0;
      while ((xl!=ncaX || yl!=ncaY || zl!=ncaZ || count!=ncaLev ))
      {
	    len=len>>1;
	  
	    index1 = 0;
	    if ( ncaX >= (len + xl) ) {
	      index1 += 1;
	      xl += len;
	      if (ncaY < (len + yl)) index1 += 2;
	    }
	    if ( ncaY >= (len + yl) ) { index1 += 1; yl += len; }
	  
	  //rotate(index1,rotation,rot_index,dim);
	  
	  //rotate_table_based(index1,current_rot,dim);
	  //temp_2d=rotations_2d[current_rot];
	  index_temp=rotations_2d[current_rot].rot_index[index1];
	  current_rot=HILBERT_2D_TABLE[current_rot][index_temp];
	  count++;
	  
      }
      
      len=len>>1;
      //Rotation2D temp=Rotation2D(rotation,rot_index);
      Rotation2D temp=rotations_2d[current_rot];
      
      if((x1-ncaX)<len && (y1-ncaY)<len)
      { // index 0
        index1= temp.rot_index[0];

      }else if((x1-ncaX)<len && (y1-ncaY)>=len)
      { // index 1
	index1=temp.rot_index[1];
      }else if((x1-ncaX)>=len && (y1-ncaY)>=len)
      { // index 2
	index1=temp.rot_index[2];
      }else if ((x1-ncaX)>=len && (y1-ncaY)<len)
      { // index 3
	index1=temp.rot_index[3];
	
      }
      
      
      if((x2-ncaX)<len && (y2-ncaY)<len)
      { // index 0
	index2=temp.rot_index[0];
	
      }else if((x2-ncaX)<len && (y2-ncaY)>=len)
      { // index 1
	index2=temp.rot_index[1];
	
      }else if((x2-ncaX)>=len && (y2-ncaY)>=len)
      { // index 2
	index2=temp.rot_index[2];
      }else if ((x2-ncaX)>=len && (y2-ncaY)<len)
      { // index 3
	index2=temp.rot_index[3];
	
      }
    
    
  }
  else if(G_dim==3)
  {
      char rotation[8]={0,1,2,3,4,5,6,7}; // Initial rotation
      char rot_index[8]={0,1,2,3,4,5,6,7}; // Initial rotation indices 
      int current_rot=0;
      while ((xl!=ncaX || yl!=ncaY || zl!=ncaZ || count!=ncaLev ) /*&& len >0*/)
      {

	len >>=1;
	
	index1 = 0;
	if ( ncaZ < (len + zl) ) {
	    if ( ncaX >= (len + xl) ) {
	      index1 += 1;
	      xl += len;
	      if (ncaY < (len + yl)) 
		index1 += 2;
	    }
	    if ( ncaY >= (len + yl) ) { 
	      index1 += 1; 
	      yl += len; 
	    }
	} else {
	    index1 = 4;
	    zl += len;
	    if ( ncaX < (len + xl) ){ 
	      index1 += 1;
	      if (ncaY < (len + yl)) 
		index1 += 2;
	    }
	    else {
	      xl += len;
	    }
	    if ( ncaY >= (len + yl) ) { 
	      index1 += 1; 
	      yl += len; 
	    }
      }
      
      //rotate(index1,rotation,rot_index,dim);
      //rotate_table_based(index1,current_rot,dim);
      
      //temp_3d=rotations_3d[current_rot];
      index_temp=rotations_3d[current_rot].rot_index[index1];
      current_rot=HILBERT_3D_TABLE[current_rot][index_temp];
      
      count++;
	
      }
      Rotation3D temp=rotations_3d[current_rot];      
      //Rotation3D temp=Rotation3D(rotation,rot_index);
      
      len >>=1;
        
       if((x1-ncaX)<len && (y1-ncaY)<len && (z1-ncaZ)<len)
	{ 
	  index1=temp.rot_index[0];
	}else if ((x1-ncaX)<len && (y1-ncaY)>=len && (z1-ncaZ)<len)
	{ 
	  index1=temp.rot_index[1]; 
	}else if((x1-ncaX)>=len && (y1-ncaY)>=len && (z1-ncaZ)<len)
	{
	  index1=temp.rot_index[2];
	}else if((x1-ncaX)>=len && (y1-ncaY)<len && (z1-ncaZ)<len)
	{ 
	  index1=temp.rot_index[3];
	}else if ((x1-ncaX)>=len && (y1-ncaY)<len && (z1-ncaZ)>=len)
	{
	  index1=temp.rot_index[4];
	}else if((x1-ncaX)>=len && (y1-ncaY)>=len && (z1-ncaZ)>=len)
	{ 
	  index1=temp.rot_index[5];
	}else if((x1-ncaX)<len && (y1-ncaY)>=len && (z1-ncaZ)>=len)
	{ 
	  index1=temp.rot_index[6];
	}else if((x1-ncaX)<len && (y1-ncaY)<len && (z1-ncaZ)>=len)
	{
	  index1=temp.rot_index[7];
	}
	
	
	if((x2-ncaX)<len && (y2-ncaY)<len && (z2-ncaZ)<len)
	{ 
	  index2=temp.rot_index[0];
	}else if ((x2-ncaX)<len && (y2-ncaY)>=len && (z2-ncaZ)<len)
	{ 
	  index2=temp.rot_index[1];
	}else if((x2-ncaX)>=len && (y2-ncaY)>=len && (z2-ncaZ)<len)
	{
	  index2=temp.rot_index[2];
	}else if((x2-ncaX)>=len && (y2-ncaY)<len && (z2-ncaZ)<len)
	{
	  index2=temp.rot_index[3];
	}else if ((x2-ncaX)>=len && (y2-ncaY)<len && (z2-ncaZ)>=len)
	{
	  index2=temp.rot_index[4];
	}else if((x2-ncaX)>=len && (y2-ncaY)>=len && (z2-ncaZ)>=len)
	{ 
	  index2=temp.rot_index[5];
	}else if((x2-ncaX)<len && (y2-ncaY)>=len && (z2-ncaZ)>=len)
	{ 
	  index2=temp.rot_index[6];
	}else if((x2-ncaX)<len && (y2-ncaY)<len && (z2-ncaZ)>=len)
	{
	  index2=temp.rot_index[7];
	}
	
	
    
  }
  
  return index1<index2; 
  
}

// #endif
