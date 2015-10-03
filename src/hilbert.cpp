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
#include "../include/hilbert.h"

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


  //NOTE: This code is directly implemented in the treenode < operator rather than calling as a function. This is here for the reference only.
  // WARNING: if you use this function you need to initialize the extern variables G_MAX_DEPTH and G_dim.

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



  unsigned int index1=0;
  unsigned int index2=0;
  unsigned int num_children=1u<<dim; // This is basically the hilbert table offset
  unsigned int rot_offset=num_children<<1;
  char index_temp=0;
  int current_rot=0;

  //unsigned int b_x,b_y,b_z;
  //unsigned int a,b,c;
  unsigned int mid_bit=G_MAX_DEPTH;

  for(int i=0; i<ncaLev;i++)
  {
    mid_bit=G_MAX_DEPTH-i-1;
//	   b_x=((ncaX&(1<<mid_bit))>>mid_bit);
//	   b_y=((ncaY&(1<<mid_bit))>>mid_bit);
//	   b_z=((ncaZ&(1<<mid_bit))>>mid_bit);
    //index1=(b_z<<2) + ((b_x^b_z)<<1) + (b_x^b_y^b_z);

    index1= (((ncaZ&(1<<mid_bit))>>mid_bit)<<2)|( (((ncaX&(1<<mid_bit))>>mid_bit)^((ncaZ&(1<<mid_bit))>>mid_bit)) <<1)|(((ncaX&(1<<mid_bit))>>mid_bit)^((ncaY&(1<<mid_bit))>>mid_bit)^((ncaZ&(1<<mid_bit))>>mid_bit));


    index_temp=rotations[rot_offset*current_rot+num_children+index1]-'0';
    current_rot=HILBERT_TABLE[current_rot*num_children+index_temp];

  }

  mid_bit--;
  index1= (((z1&(1<<mid_bit))>>mid_bit)<<2)|( (((x1&(1<<mid_bit))>>mid_bit)^((z1&(1<<mid_bit))>>mid_bit)) <<1)|(((x1&(1<<mid_bit))>>mid_bit)^((y1&(1<<mid_bit))>>mid_bit)^((z1&(1<<mid_bit))>>mid_bit));
  index2= (((z2&(1<<mid_bit))>>mid_bit)<<2)|( (((x2&(1<<mid_bit))>>mid_bit)^((z2&(1<<mid_bit))>>mid_bit)) <<1)|(((x2&(1<<mid_bit))>>mid_bit)^((y2&(1<<mid_bit))>>mid_bit)^((z2&(1<<mid_bit))>>mid_bit));



  return rotations[rot_offset*current_rot+num_children+index1] < rotations[rot_offset*current_rot+num_children+index2];



}

// #endif