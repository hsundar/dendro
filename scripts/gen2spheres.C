// this program builds the sequence of points on the surfaces of 2 concentric spheres in the unit cube. Center of bothspheres is (0.5,0.5,0.5)
// Command line parameters: <fileName> <maxDepth> <regularDepth> [<innerRaduis>] [<outerRadius>]
// Default for innerRadius is 0.2
// Default for outerRadius is 0.4
//

#include "TreeNode.h"
#include <iostream>
#include <algorithm>
#include "externVars.h"

inline double MaxDistSq(
    const ot::TreeNode & oct, 
    const double c_x,
    const double c_y,
    const double c_z 
    )
{
  using namespace std;
  
  double s; 
  double max_dist = 0; /* will store square of distance from center to farthest point on cube  */ 

  double octLen = ldexp(1.0,-oct.getLevel());
  unsigned maxD = oct.getMaxDepth();
  double octX = ldexp(static_cast<double>(oct.getX()),-maxD);
  double octY = ldexp(static_cast<double>(oct.getY()),-maxD);
  double octZ = ldexp(static_cast<double>(oct.getZ()),-maxD);
    
  s = max(fabs(c_x - octX), fabs(octX+octLen - c_x));
  max_dist += s*s;
    
  s = max(fabs(c_y - octY), fabs(octY+octLen - c_y));
  max_dist += s*s;
  
  s = max(fabs(c_z - octZ), fabs(octZ+octLen - c_z));
  max_dist += s*s;
  
  return max_dist;
}

inline double MinDistSq(
    const ot::TreeNode & oct, 
    const double c_x,
    const double c_y,
    const double c_z 
    )
{
  using namespace std;
  
  double s; 
  double min_dist = 0; /* will store square distance from center to nearest point on cube   */

  double octLen = ldexp(1.0,-oct.getLevel());
  unsigned maxD = oct.getMaxDepth();
  double octX = ldexp(static_cast<double>(oct.getX()),-maxD);
  double octY = ldexp(static_cast<double>(oct.getY()),-maxD);
  double octZ = ldexp(static_cast<double>(oct.getZ()),-maxD);
    
  if (c_x < octX)
  {
    s = octX - c_x;
    min_dist += s*s;
  } 
  else if (c_x > octX+octLen)
  {
    s = c_x - octX - octLen;
    min_dist += s*s;
  }
    
  if (c_y < octY)
  {
    s = octY - c_y;
    min_dist += s*s;
  } 
  else if (c_y > octY+octLen)
  {
    s = c_y - octY - octLen;
    min_dist += s*s;
  }
  
  if (c_z < octZ)
  {
    s = octZ - c_z;
    min_dist += s*s;
  } 
  else if (c_z > octZ+octLen)
  {
    s = c_z - octZ - octLen;
    min_dist += s*s;
  }
  
  return min_dist ;

}


int main(int argc, char ** argv )
{
  using namespace std;
  
  unsigned maxD, regD;
  double inneR=0.2, outeR=0.4;

  // center of sphere
  const double c_x = 0.5;
  const double c_y = 0.5;
  const double c_z = 0.5;
  
  char * fName;
  
  switch (argc)
  {
    case 6:
      outeR=atof(argv[5]);
    case 5:
      inneR=atof(argv[4]);
    case 4:
      regD=atoi(argv[3]);
      maxD=atoi(argv[2]);
      fName = argv[1];
      break;
    
    default:
      cout<<"Usage: "<<argv[0]<<" <fileName> <maxDepth> <regularDepth> [<innerRaduis>] [<outerRadius>]"<<endl;
      return(1);
  }

  cout<<"Using these settings: maxDepth="<<maxD<<" regularDepth="<<regD<<" innerRaduis="<<inneR<<" outerRadius="<<outeR<<endl;
  
  if (outeR>1 || inneR<=0 || outeR<=inneR || regD>maxD)
  {
    cout<<"Erroneous arguments"<<endl;
    return(1);
  }

  // we don't do any pre-allocation;  hope STL will do a good job guessing when to preallocate
  vector<ot::TreeNode> buffer;
  vector<ot::TreeNode> output;
  
  // push root octant into buffer
  buffer.push_back(ot::TreeNode(3,maxD));
  
  double outeS=outeR*outeR;
  double inneS=inneR*inneR;
  
  while (!buffer.empty())
  {
    ot::TreeNode curr = buffer.back();
    buffer.pop_back();
    
    double currMaxDist =  MaxDistSq( curr, c_x, c_y, c_z);  // SQUARE of distance
    if(currMaxDist<=inneS)
      continue; // octant inside (maybe touches) inner sphere - drop it
    
    double currMinDist =  MinDistSq( curr, c_x, c_y, c_z);  // SQUARE!
    if (currMinDist>=outeS)
      continue; // octant outside (maybe touches) outer sphere -- drop it

    if (currMaxDist>outeS || currMinDist<inneS) // if the octant intersects (not just touches) either sphere
    {
      if (curr.getLevel()<maxD)
	curr.addChildren(buffer);
      else
	output.push_back(curr);
      
      continue;
    }
    
    // the contant lies between two spheres (maybe touches one or both)
    if (curr.getLevel()<regD)
      curr.addChildren(buffer);
    else
      output.push_back(curr);
  }

  cout<<"Produced "<<output.size()<<" octants."<<endl; 
  ot::writeNodesToFile(fName, output);
}

