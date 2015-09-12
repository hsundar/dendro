/*
 *@author: Milinda Fernando
 *@date: 09/04/2015 //// This is refactored code from HilbertBenchmark code.
 *School of Computing, University of Utah
 * 
 * Contains the Ordering Implementations of Morton Naive implementation and NCA based Method. 
 * 
 */


/*
 * Baseline Dendro Morton implementation
 * This is the current implementation in the dendro. 
 */

#include "../../include/hilbert/morton.h"

bool morton_order(Point p1, Point p2)
{
   
    if(p1==p2)
    {
	return false;
    }
   
    unsigned int x = (p1.xint()^p2.xint());
    unsigned int y = (p1.yint()^p2.yint());
    unsigned int z = (p1.zint()^p2.zint());

    //Default pref: z > y > x.
    unsigned int maxC = z;
    unsigned int yOrx = y;
    if(yOrx < x) { if( (x^yOrx) >= yOrx ) {yOrx = x;} }
    if(maxC < yOrx) { if( (maxC^yOrx) >= maxC ) {maxC = yOrx;} }

    if(maxC == z) { 
      if (p1.z() < p2.z())
	return true;
      else
	return false;
      
    }
    else if(maxC == y) { 
      if(p1.y() < p2.y())
	return true;
      else 
	return false;
      
    }
    else {  
      if(p1.x() < p2.x())
	return true;
      else
	return false;
      
    }  
   
   
 
}

bool morton_order_NCA(const Point& p1,const Point& p2)
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
  
  unsigned int maxDiffBinLen = binOp::binLength(maxDiff);
  //Eliminate the last maxDiffBinLen bits.
  unsigned int ncaX = ((x1>>maxDiffBinLen)<<maxDiffBinLen);
  unsigned int ncaY = ((y1>>maxDiffBinLen)<<maxDiffBinLen);
  unsigned int ncaZ = ((z1>>maxDiffBinLen)<<maxDiffBinLen);
  unsigned int ncaLev = (maxDepth - maxDiffBinLen);
  
  unsigned int xl=0;
  unsigned int yl=0;
  unsigned int zl=0;
  
  unsigned int len=1<<G_MAX_DEPTH;
   
  len=len/(1<<(ncaLev+1));
  unsigned int index1=0;
  unsigned int index2=0;
  
  if(G_dim==2){
   
  
    index1 = 0;
    if ( x1>=(len+ncaX) ) {
	index1 += 1;
      if (y1>=(len+ncaY)){ 
	  index1 += 2;
      }
    }else if ( y1>=(len+ncaY) ) { 
      index1 += 2;
      
    }

    index2 = 0;
    if ( x2>=(len+ncaX) ) {
	index2 += 1;
      if (y2>=(len+ncaY)){
	index2 += 2;
      }
    }else if ( y2>=(len+ncaY) ) {
      index2 += 2;
      
    }
    
    
  }else if(G_dim==3)
  {
      
      index1=0;
      if ( z1 < (len + ncaZ) ) {
	 if ( x1 >= (len + ncaX) ) {
	    index1 += 1;
	    if (y1 >= (len + ncaY)) 
		index1 += 2;
	 }else if ( y1>= (len + ncaY) ) { 
	    index1 += 2; 
	 }
      } else {
	    index1 = 4;
	 if ( x1 >= (len + ncaX) ) {
	    index1 += 1;
	    if (y1 >= (len + ncaY)) 
	      index1 += 2;
	 }else if ( y1>= (len + ncaY) ) { 
	    index1 += 2; 
	 }
      }
      
      index2=0;
      if ( z2 < (len + ncaZ) ) {
	 if ( x2 >= (len + ncaX) ) {
	    index2 += 1;
	    if (y2 >= (len + ncaY)) 
		index2 += 2;
	 }else if ( y2>= (len + ncaY) ) { 
	    index2 += 2; 
	 }
      } else {
	    index2 = 4;
	 if ( x2 >= (len + ncaX) ) {
	    index2 += 1;
	    if (y2 >= (len + ncaY)) 
	      index2 += 2;
	 }else if ( y2>= (len + ncaY) ) { 
	    index2 += 2; 
	 }
      }
    

  }
  return index1<index2;
  
  
}
