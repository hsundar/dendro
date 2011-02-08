
/**
  @file TreeNode.txx
  @author Rahul S. Sampath, rahul.sampath@gmail.com
 */

#ifdef __DEBUG__
#ifndef __DEBUG_TN__
#define __DEBUG_TN__
#endif
#endif

#include "colors.h"
#include <cassert>

namespace ot {
  // Inline Functions ...
  inline unsigned char TreeNode::getChildNumber() const {
    unsigned int len = (1u << (m_uiMaxDepth - getLevel()));
    unsigned int len_par = (1u << (m_uiMaxDepth - getLevel() + 1u));

    unsigned int i = (m_uiX % len_par);
    unsigned int j = (m_uiY % len_par);
    unsigned int k = (m_uiZ % len_par);
    i /= len;
    j /= len;
    k /= len;

    unsigned char childNum = static_cast<unsigned char>(4*k + 2*j + i);
    return childNum;

  }//end function

  inline int TreeNode::getAnchor(unsigned int &x, unsigned int&y, unsigned int&z) const {
    x = m_uiX;
    y = m_uiY;
    z = m_uiZ;
    return 1;
  }

  inline bool TreeNode  :: operator == ( TreeNode   const & other)  const{
#ifdef __DEBUG_TN__
    if( ((this->m_uiDim) != (other.m_uiDim)) || ((this->m_uiMaxDepth) != (other.m_uiMaxDepth)) ) {
      std::cout<<"Me: "<<(*this)<<" Other: "<<other<<std::endl;
      std::cout<<"My Dim: "<<m_uiDim<<" OthDim: "<<other.m_uiDim<<" My MaxD: "<<m_uiMaxDepth<<" othMD: "<<other.m_uiMaxDepth<<std::endl;
      assert( false );
    }
#endif
    if( (m_uiX == other.m_uiX) && (m_uiY == other.m_uiY) && (m_uiZ == other.m_uiZ) &&
        ( (m_uiLevel & ot::TreeNode::MAX_LEVEL) == (other.m_uiLevel & ot::TreeNode::MAX_LEVEL) ) ) {
      return true;
    }else {
      return false;
    }
  }//end fn.

  inline bool TreeNode  :: operator != (TreeNode  const  & other)  const{
#ifdef __DEBUG_TN__
    if( ((this->m_uiDim) != (other.m_uiDim)) || ((this->m_uiMaxDepth) != (other.m_uiMaxDepth)) ) {
      std::cout<<"Me: "<<(*this)<<" Other: "<<other<<std::endl;
      std::cout<<"My Dim: "<<m_uiDim<<" OthDim: "<<other.m_uiDim<<" My MaxD: "<<m_uiMaxDepth<<" othMD: "<<other.m_uiMaxDepth<<std::endl;
      assert(false);
    }
#endif
    return (!((*this) == other));
  }//end fn.

  inline bool TreeNode  :: operator  < (TreeNode   const & other)  const{
#ifdef __DEBUG_TN__
    if( ((this->m_uiDim) != (other.m_uiDim)) || ((this->m_uiMaxDepth) != (other.m_uiMaxDepth)) ) {
      std::cout<<"Me: "<<(*this)<<" Other: "<<other<<std::endl;
      std::cout<<"My Dim: "<<m_uiDim<<" OthDim: "<<other.m_uiDim<<" My MaxD: "<<m_uiMaxDepth<<" othMD: "<<other.m_uiMaxDepth<<std::endl;
      assert(false);
    }
#endif
    // first compare the x, y, and z to determine which one dominates ...
    //Ancestor is smaller.
    if( (this->m_uiX == other.m_uiX) && (this->m_uiY == other.m_uiY) && (this->m_uiZ == other.m_uiZ) ) {
      return ( (this->m_uiLevel & ot::TreeNode::MAX_LEVEL) < (other.m_uiLevel & ot::TreeNode::MAX_LEVEL) );	
    }//end if	

    unsigned int x = (m_uiX^other.m_uiX);
    unsigned int y = (m_uiY^other.m_uiY);
    unsigned int z = (m_uiZ^other.m_uiZ);

    //Default pref: z > y > x.
    unsigned int maxC = z;
    unsigned int yOrx = y;
    if(yOrx < x) { if( (x^yOrx) >= yOrx ) {yOrx = x;} }
    if(maxC < yOrx) { if( (maxC^yOrx) >= maxC ) {maxC = yOrx;} }

    if(maxC == z) { return (m_uiZ < other.m_uiZ); }
    else if(maxC == y) {  return (m_uiY < other.m_uiY); }
    else {  return (m_uiX < other.m_uiX); }
  }//end function

  inline bool TreeNode  :: operator  <= ( TreeNode  const  & other)  const{
#ifdef __DEBUG_TN__
    if( ((this->m_uiDim) != (other.m_uiDim)) || ((this->m_uiMaxDepth) != (other.m_uiMaxDepth)) ) {
      std::cout<<"Me: "<<(*this)<<" Other: "<<other<<std::endl;
      std::cout<<"My Dim: "<<m_uiDim<<" OthDim: "<<other.m_uiDim<<" My MaxD: "<<m_uiMaxDepth<<" othMD: "<<other.m_uiMaxDepth<<std::endl;
      assert(false);
    }
#endif
    return ( ((*this) < other) || ((*this) == other) );
  }//end fn.

  inline bool TreeNode  :: operator  > (TreeNode  const  & other)  const{
#ifdef __DEBUG_TN__
    if( ((this->m_uiDim) != (other.m_uiDim)) || ((this->m_uiMaxDepth) != (other.m_uiMaxDepth)) ) {
      std::cout<<"Me: "<<(*this)<<" Other: "<<other<<std::endl;
      std::cout<<"My Dim: "<<m_uiDim<<" OthDim: "<<other.m_uiDim<<" My MaxD: "<<m_uiMaxDepth<<" othMD: "<<other.m_uiMaxDepth<<std::endl;
      assert(false);
    }
#endif
    return ( (!((*this) < other)) && (!((*this) == other)) );
  }//end fn.

  inline bool TreeNode  :: operator  >= ( TreeNode  const & other)  const{
#ifdef __DEBUG_TN__
    if( ((this->m_uiDim) != (other.m_uiDim)) || ((this->m_uiMaxDepth) != (other.m_uiMaxDepth)) ) {
      std::cout<<"Me: "<<(*this)<<" Other: "<<other<<std::endl;
      std::cout<<"My Dim: "<<m_uiDim<<" OthDim: "<<other.m_uiDim<<" My MaxD: "<<m_uiMaxDepth<<" othMD: "<<other.m_uiMaxDepth<<std::endl;
      assert(false);
    }
#endif
    return ( !((*this) < other) ) ;
  }//end fn.

  inline TreeNode TreeNode  ::getParent() const {
    //For any node at level l, the last (maxD-l) bits are 0. 
    //By convention, root's parent is also root.
    unsigned int parX,parY,parZ;
    unsigned int parLev = (( (m_uiLevel & ot::TreeNode::MAX_LEVEL) > 0)
        ? ((m_uiLevel & ot::TreeNode::MAX_LEVEL)-1) : 0);
    parX = ( ( m_uiX >> ( m_uiMaxDepth - parLev ) ) << ( m_uiMaxDepth - parLev ) );
    parY = ( ( m_uiY >> ( m_uiMaxDepth - parLev ) ) << ( m_uiMaxDepth - parLev ) );
    parZ = ( ( m_uiZ >> ( m_uiMaxDepth - parLev ) ) << ( m_uiMaxDepth - parLev ) );
    return TreeNode (1,parX,parY,parZ,parLev,m_uiDim,m_uiMaxDepth);
  }//end function

  inline TreeNode TreeNode  ::getAncestor(unsigned int ancLev) const {
    //For any node at level l, the last (maxD-l) bits are 0. 
    unsigned int ancX,ancY,ancZ;
    ancX = ( ( m_uiX >> ( m_uiMaxDepth - ancLev ) ) << ( m_uiMaxDepth - ancLev ) );
    ancY = ( ( m_uiY >> ( m_uiMaxDepth - ancLev ) ) << ( m_uiMaxDepth - ancLev ) );
    ancZ = ( ( m_uiZ >> ( m_uiMaxDepth - ancLev ) ) << ( m_uiMaxDepth - ancLev ) );
    return TreeNode (1,ancX,ancY,ancZ,ancLev,m_uiDim,m_uiMaxDepth);
  }//end function

  inline int TreeNode::setWeight(unsigned int w) {
    m_uiWeight = w;
    return 1;
  }    

  inline int TreeNode::addWeight(unsigned int w) {
    m_uiWeight += w;
    return 1;
  }    

  inline unsigned int TreeNode::getDim() const {
    return m_uiDim;
  }

  inline unsigned int TreeNode::getMaxDepth() const {
    return m_uiMaxDepth;
  }

  inline unsigned int TreeNode::getLevel() const { 
    return (m_uiLevel & ot::TreeNode::MAX_LEVEL);   
  }

  inline unsigned int TreeNode::getFlag() const { 
    return m_uiLevel;   
  }

  inline int TreeNode::orFlag(unsigned int w) { 
    m_uiLevel = (m_uiLevel | w);
    return 1;   
  }

  inline int TreeNode::setFlag(unsigned int w) { 
    m_uiLevel = w;
    return 1;   
  }

  inline unsigned int TreeNode::getWeight() const {
    return m_uiWeight;
  }   

  inline unsigned int TreeNode::getX() const {
    return m_uiX;
  }

  inline unsigned int TreeNode::getY() const {
    return m_uiY;
  }

  inline unsigned int TreeNode::getZ() const {
    return m_uiZ;
  }

  inline unsigned int TreeNode::getParentX() const {
    return getParent().getX();
  }

  inline unsigned int TreeNode::getParentY() const {
    return getParent().getY();
  }

  inline unsigned int TreeNode::getParentZ() const {
    return getParent().getZ();
  }

  inline TreeNode TreeNode  ::getDFD() const {
    TreeNode dfd(1,m_uiX,m_uiY,m_uiZ,m_uiMaxDepth,m_uiDim,m_uiMaxDepth);
    return dfd;
  }//end function

  inline TreeNode TreeNode  ::getDLD() const {
    TreeNode dld(1,maxX()-1,maxY()-1,maxZ()-1,m_uiMaxDepth,m_uiDim,m_uiMaxDepth);
    return dld;
  }//end function

  inline TreeNode TreeNode::getNext() const {

    TreeNode m=*this;
    unsigned int mask = (1u << (m_uiMaxDepth - getLevel()));

    int i;
    for(i=m.m_uiLevel;i>=0;i--){
      m.m_uiX=(m.m_uiX^mask);
      if((m.m_uiX & mask))
	break;
      m.m_uiY=(m.m_uiY^mask);
      if((m.m_uiY & mask))
	break;
      m.m_uiZ=(m.m_uiZ^mask);
      if((m.m_uiZ & mask))
	break;
      mask=(mask<<1);
    }
    m.m_uiLevel=i;
    return m;

  }

  inline TreeNode TreeNode::getFirstChild() const{
    TreeNode m=*this;
    m.m_uiLevel++;
    return m;
  }

  inline bool TreeNode :: isAncestor(TreeNode const & other) const {
#ifdef __DEBUG_TN__
    if( ((this->m_uiDim) != (other.m_uiDim)) || ((this->m_uiMaxDepth) != (other.m_uiMaxDepth)) )  {
      std::cout<<"Me: "<<(*this)<<" Other: "<<other<<std::endl;
      std::cout<<"My Dim: "<<m_uiDim<<" OthDim: "<<other.m_uiDim<<" My MaxD: "<<m_uiMaxDepth<<" othMD: "<<other.m_uiMaxDepth<<std::endl;
      assert(false);
    }
#endif
    return ( (other > (*this)) && (other <= (this->getDLD())) );
  }//end function

  inline unsigned int TreeNode   :: minX() const{	
    return getX();
  }//end fn.

  inline  unsigned int TreeNode   ::  minY()  const{	
    if(m_uiDim < 2) {return 0;}
    return getY();
  }//end fn.

  inline  unsigned int TreeNode   ::  minZ() const {	
    if(m_uiDim < 3) {return 0;}
    return getZ();
  }//end fn.

  inline  unsigned int TreeNode   ::  maxX() const { 	
    unsigned int len = (1u << (m_uiMaxDepth - getLevel()));
    return (minX() + len);
  }//end fn.

  inline  unsigned int TreeNode   ::  maxY() const {	
    if(m_uiDim < 2) {return 1;}
    unsigned int len = (1u << (m_uiMaxDepth - getLevel()));
    return (minY() + len);
  }//end fn.

  inline  unsigned int TreeNode   ::  maxZ() const {	
    if(m_uiDim < 3) {return 1;}
    unsigned int len = (1u << (m_uiMaxDepth - getLevel()));
    return (minZ() + len);
  }//end fn.

  inline  TreeNode   TreeNode   ::  getLeft() const {
    //-ve X
    if( minX() > 0 ) {
      unsigned int len =  (1u << (m_uiMaxDepth - getLevel()));
      unsigned int xres = (minX() - len);
      unsigned int yres = minY();
      unsigned int zres = minZ();
      TreeNode res(1, xres, yres, zres, getLevel(), m_uiDim, m_uiMaxDepth);
      return res;
    }else {
      TreeNode res(m_uiDim, m_uiMaxDepth);
      return res;
    }
  }//end fn.

  inline  TreeNode   TreeNode   ::  getRight() const  {
    //+ve X
    if( maxX() < (1u << m_uiMaxDepth) ) {
    unsigned int xres = maxX();
     unsigned int yres = minY();
     unsigned int zres = minZ();
      TreeNode  res(1, xres, yres, zres, getLevel(), m_uiDim, m_uiMaxDepth);
      return res;
    }else {
      TreeNode res(m_uiDim, m_uiMaxDepth);
      return res;
    }
  }//end fn.
  
  inline  TreeNode   TreeNode   ::  getBack() const {
    //+ve Y
    if( (m_uiDim > 1) && (maxY() < (1u << m_uiMaxDepth)) ) {
     unsigned int xres = minX();
     unsigned int yres = maxY();
     unsigned int zres = minZ();
      TreeNode res(1, xres, yres, zres, getLevel(), m_uiDim, m_uiMaxDepth);
      return res;
    }else {
      TreeNode res (m_uiDim, m_uiMaxDepth);
      return res;
    }
  }//end fn.

  inline TreeNode   TreeNode   ::  getFront() const {
    //-ve Y
    if( minY() > 0 ) {
     unsigned int len =  (1u << (m_uiMaxDepth - getLevel()));
    unsigned int xres = minX();
     unsigned int yres = minY() - len;
     unsigned int zres = minZ();
      TreeNode res(1, xres, yres, zres, getLevel(), m_uiDim, m_uiMaxDepth);
      return res;
    }else {
      TreeNode res(m_uiDim, m_uiMaxDepth);
      return res;
    }
  }//end fn.

  inline TreeNode   TreeNode   ::  getBottom() const {
    //-ve Z
    if( minZ() > 0 ) {
    unsigned int len = (1u << (m_uiMaxDepth - getLevel()));
     unsigned int xres = minX();
     unsigned int yres = minY();
     unsigned int zres = minZ() - len;
      TreeNode res(1, xres, yres, zres, getLevel(), m_uiDim, m_uiMaxDepth);
      return res;
    }else {
      TreeNode res(m_uiDim, m_uiMaxDepth);
      return res;
    }
  }//end fn.

  inline  TreeNode   TreeNode   ::  getTop() const {
    //+ve Z
    if( (m_uiDim == 3) && (maxZ() < (1u << m_uiMaxDepth)) ) {
    unsigned int xres = minX();
     unsigned int yres = minY() ;
     unsigned int zres = maxZ();
      TreeNode res(1, xres, yres, zres, getLevel(), m_uiDim, m_uiMaxDepth);
      return res;
    }else {
      TreeNode res(m_uiDim, m_uiMaxDepth);
      return res;
    }
  }//end fn.

  inline  TreeNode   TreeNode   ::  getLeftBack() const {
    return (getLeft().getBack());
  }//end fn.

  inline  TreeNode    TreeNode   :: getRightBack()  const{
    return (getRight().getBack());
  }//end fn.

  inline  TreeNode   TreeNode   ::  getLeftFront() const {
     return (getLeft().getFront());
  }//end fn.

  inline TreeNode   TreeNode   ::  getRightFront() const {
     return (getRight().getFront());
  }//end fn.

  inline TreeNode   TreeNode   ::  getBottomLeft() const {
    return (getBottom().getLeft());
  }//end fn.

  inline  TreeNode   TreeNode   ::  getBottomRight() const {
    return (getBottom().getRight());
  }//end fn.

  inline  TreeNode   TreeNode   ::  getBottomBack() const {
    return (getBottom().getBack());
  }//end fn.
  
  inline TreeNode   TreeNode   ::  getBottomFront() const {
    return (getBottom().getFront());
  }//end fn.

  inline TreeNode   TreeNode   ::  getBottomLeftBack() const {
    return (getBottomLeft().getBack());
  }//end fn.

  inline  TreeNode   TreeNode   ::  getBottomRightBack() const {
    return (getBottomRight().getBack());
  }//end fn.

  inline  TreeNode   TreeNode   ::  getBottomLeftFront() const{
    return (getBottomLeft().getFront());
  }//end fn.

  inline TreeNode   TreeNode   ::  getBottomRightFront()  const{
    return (getBottomRight().getFront());
  }//end fn.

  inline  TreeNode   TreeNode   ::  getTopLeft() const {
    return (getTop().getLeft());
  }//end fn.

  inline   TreeNode    TreeNode   :: getTopRight() const {
    return (getTop().getRight());
  }//end fn.

  inline TreeNode    TreeNode   :: getTopBack() const {
    return (getTop().getBack());
  }//end fn.

  inline  TreeNode   TreeNode   ::  getTopFront() const {
    return (getTop().getFront());
  }//end fn.

  inline  TreeNode   TreeNode   ::  getTopLeftBack() const {
    return (getTopLeft().getBack());
  }//end fn.

  inline TreeNode    TreeNode   :: getTopRightBack() const {
    return (getTopRight().getBack());
  }//end fn.

  inline  TreeNode   TreeNode   ::  getTopLeftFront() const {
    return (getTopLeft().getFront());
  }//end fn.

  inline TreeNode   TreeNode   ::  getTopRightFront() const {
    return (getTopRight().getFront());
  }//end fn.
  
  ///////////////////////////////////////////////////////////////////////////////////////
  //Returns a Vector of Candidates for Balanced Neigbors in a particular direction.
  //The results are returned in the following order: 1) Same level, 2) Parent's
  //level and 3) Children's level.
  //If there is no candidate at the parent's or child's level, then the candidate at the
  //same level is returned instead for that position
  //Note: For the neighbours at the parent's level, the dimension is not
  //considered explicitly. It is implicit in the childNum.

  inline std::vector<TreeNode >  TreeNode   :: getB_Left() const {
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"L"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(6);
    }else if(dim == 2) {
      it.resize(4);
    }else {
      it.resize(3);
    }

    it[0] = getLeft();
    unsigned int childNum = getChildNumber();

    //0,2,4,6
    if((childNum == 0) || (childNum == 2) || (childNum == 4) || (childNum == 6) ){
      it[1] = getParent().getLeft();	
    }else {
      it[1] = it[0]; 
    }	

    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[0].getLeft();
    if(dim > 1) {
      it[3] = children[2].getLeft();
    }
    if(dim > 2) {
      it[4] = children[4].getLeft();	
      it[5] = children[6].getLeft();
    }
    children.clear();
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_Right() const {
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"R"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(6);
    }else if(dim == 2) {
      it.resize(4);
    }else {
      it.resize(3);
    }
    //1,3,5,7
    it[0] = getRight();
    unsigned int childNum = getChildNumber();
    if((childNum == 1) || (childNum == 3) || (childNum == 5) || (childNum == 7)){
      it[1] = getParent().getRight();
    }else {
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[1].getRight();
    if(dim > 1) {
      it[3] = children[3].getRight();
    }
    if(dim >2) {
      it[4] = children[5].getRight();	
      it[5] = children[7].getRight();
    }
    children.clear();
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_Top() const {
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"T"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(6);
    }else if(dim == 2) {
      it.resize(0);
      return it;
    }else {
      it.resize(0);
      return it;
    }
    //4,5,6,7
    it[0] = getTop();
    unsigned int childNum = getChildNumber();
    if( (childNum  == 4) || (childNum == 5) || (childNum  == 6) || (childNum == 7) ) {
      it[1] = getParent().getTop();	
    }else {
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[4].getTop();	
    it[3] = children[5].getTop();
    it[4] = children[6].getTop();	
    it[5] = children[7].getTop();
    children.clear();
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_Bottom() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"Bo"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(6);
    }else if(dim == 2) {
      it.resize(0);
      return it;
    }else {
      it.resize(0);
      return it;
    }
    //0,1,2,3
    it[0] = getBottom();
    unsigned int childNum = getChildNumber();
    if( (childNum == 0) || (childNum == 1) || ( childNum == 2) || ( childNum == 3) ) {
      it[1] = getParent().getBottom();	
    }else {
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[0].getBottom();	 
    it[3] = children[1].getBottom();
    it[4] = children[2].getBottom();	
    it[5] = children[3].getBottom();
    children.clear();
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_Front() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"F"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(6);
    }else if(dim == 2) {
      it.resize(4);
    }else {
      it.resize(0);
      return it;
    }
    //0,1,4,5
    it[0] = getFront();
    unsigned int childNum = getChildNumber();
    if( (childNum == 0) || (childNum == 1) ||(childNum == 4) || (childNum == 5) ) {
      it[1] = getParent().getFront();	
    }else {
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[0].getFront();	
    it[3] = children[1].getFront();	
    if(dim >2) {
      it[4] = children[4].getFront();	
      it[5] = children[5].getFront();
    }
    children.clear();
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_Back() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"Bk"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(6);
    }else if(dim == 2) {
      it.resize(4);
    }else {
      it.resize(0);
      return it;
    }
    //2,3,6,7
    it[0] = getBack();
    unsigned int childNum = getChildNumber();
    if( (childNum == 2) || (childNum == 3) || (childNum == 6)||(childNum == 7) ) {
      it[1] = getParent().getBack();
    }else {
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[2].getBack();	
    it[3] = children[3].getBack();	
    if(dim == 3) {
      it[4] = children[6].getBack();	
      it[5] = children[7].getBack();
    }
    children.clear(); 
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_TopLeft() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"TL"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(4);
    }else if(dim == 2) {
      it.resize(0);
      return it;
    }else {
      it.resize(0);
      return it;
    }
    //4,6 -> TL
    //5,7->T but not L
    //0,2->L but not T
    it[0] = getTopLeft();
    unsigned int childNum = getChildNumber();
    if( (childNum == 4) || (childNum == 6) ) {
      it[1] = getParent().getTopLeft();	
    }else if( (childNum == 5) || (childNum == 7) ){
      it[1] = getParent().getTop();	
    }else if( (childNum == 0) || (childNum == 2) ){
      it[1] = getParent().getLeft();	
    }else {
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[4].getTopLeft();
    it[3] = children[6].getTopLeft();
    children.clear();
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_TopRight() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"TR"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(4);
    }else if(dim == 2) {
      it.resize(0);
      return it;
    }else {
      it.resize(0);
      return it;
    }
    //5,7 TR
    //4,6 T but not R
    //1,3 R but not T
    it[0] = getTopRight();
    unsigned int childNum = getChildNumber();
    if( (childNum == 5) || (childNum == 7) ) {
      it[1] = getParent().getTopRight();	
    }else if( (childNum == 4) || (childNum == 6) ){
      it[1] = getParent().getTop();	
    }else if( (childNum == 1) || (childNum == 3) ){
      it[1] = getParent().getRight();	
    }else {	
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[5].getTopRight();
    it[3] = children[7].getTopRight();
    children.clear(); 
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_BottomLeft() const{ 
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"BoL"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(4);
    }else if(dim == 2) {
      it.resize(0);
      return it;
    }else {
      it.resize(0);
      return it;
    }
    //0,2 BL
    //4,6 L but not B
    //1,3 B but not L
    it[0] = getBottomLeft();
    unsigned int childNum = getChildNumber();
    if( (childNum == 0) || (childNum == 2) ) {
      it[1] = getParent().getBottomLeft();	
    }else if( (childNum == 1) || (childNum == 3) ){
      it[1] = getParent().getBottom();	
    }else if( (childNum == 4) || (childNum == 6) ){
      it[1] = getParent().getLeft();	
    }else {
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[0].getBottomLeft();
    it[3] = children[2].getBottomLeft();
    children.clear();
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_BottomRight() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"BoR"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(4);
    }else if(dim == 2) {
      it.resize(0);
      return it;
    }else {
      it.resize(0);
      return it;
    }
    //1,3 BR
    //0,2 B but not R
    //5,7 R but not B
    it[0] = getBottomRight();
    unsigned int childNum = getChildNumber();
    if( (childNum == 1) || (childNum == 3) ) {
      it[1] = getParent().getBottomRight();	
    }else if( (childNum == 0) || (childNum == 2) ){
      it[1] = getParent().getBottom();	
    }else if( (childNum == 5) || (childNum == 7) ){
      it[1] = getParent().getRight();	
    }else {
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[1].getBottomRight();
    it[3] = children[3].getBottomRight();
    children.clear();
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_LeftFront() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"LF"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(4);
    }else if(dim == 2) {
      it.resize(3);
    }else {
      it.resize(0);
      return it;
    }
    //0,4 LF
    //1,5 F but not L
    //2,6 L but not F
    it[0] = getLeftFront();
    unsigned int childNum = getChildNumber();

    if((childNum == 0) || (childNum == 4)) {
      it[1] = getParent().getLeftFront();
    }else if((childNum == 1) || (childNum == 5)) {
      it[1] = getParent().getFront();
    }else if((childNum == 2) || (childNum == 6)) {
      it[1] = getParent().getLeft();
    }
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[0].getLeftFront();
    if(dim >2) {	
      it[3] = children[4].getLeftFront();
    }
    children.clear();
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_RightFront() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"RF"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(4);
    }else if(dim == 2) {
      it.resize(3);
    }else {
      it.resize(0);
      return it;
    }
    //1,5 RF
    //0,4 F but not R
    //3,7 R but not F
    it[0] = getRightFront();
    unsigned int childNum = getChildNumber();
    if( (childNum == 1) || (childNum == 5) ) {
      it[1] = getParent().getRightFront();	
    }else if( (childNum == 0) || (childNum == 4) ) {
      it[1] = getParent().getFront();	
    }else if((childNum == 3) || (childNum == 7)) {
      it[1] = getParent().getRight();	
    }	

    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[1].getRightFront();
    if(dim > 2) {
      it[3] = children[5].getRightFront();
    }
    children.clear();
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_TopFront() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"TF"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(4);
    }else if(dim == 2) {
      it.resize(0);
      return it;
    }else {
      it.resize(0);
      return it;
    }
    //4,5 TF
    //6,7 T but not F
    //0,1 F but not T
    it[0] = getTopFront();
    unsigned int childNum = getChildNumber();
    if( (childNum == 4) || (childNum == 5) ) {
      it[1] = getParent().getTopFront();	
    }else if( (childNum == 6) || (childNum == 7) ){
      it[1] = getParent().getTop();	
    }else if( (childNum == 0) || (childNum == 1) ){
      it[1] = getParent().getFront();	
    }else {
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[4].getTopFront();	
    it[3] = children[5].getTopFront();
    children.clear();
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_BottomFront() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"BoF"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(4);
    }else if(dim == 2) {
      it.resize(0);
      return it;
    }else {
      it.resize(0);
      return it;
    }
    //0,1 BoF
    //2,3 Bo but not F
    //4,5 F but not Bo
    it[0] = getBottomFront();
    unsigned int childNum = getChildNumber();
    if( (childNum == 0) || (childNum == 1) ) {
      it[1] = getParent().getBottomFront();	
    }else if( (childNum == 2) || (childNum == 3) ){
      it[1] = getParent().getBottom();	
    }else if( (childNum == 4) || (childNum == 5) ){
      it[1] = getParent().getFront();	
    }else {
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[0].getBottomFront();	 
    it[3] = children[1].getBottomFront();
    children.clear();
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_TopLeftFront() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"TLF"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(3);
    }else if(dim == 2) {
      it.resize(0);
      return it;
    }else {
      it.resize(0);
      return it;
    }
    //4 TLF
    //6 TL but not F
    //0 LF but not T
    //5 TF but not L
    //7 T but not LF
    //2 L but not TF
    //1 F but not TL
    it[0] = getTopLeftFront();
    unsigned int childNum = getChildNumber();
    if( (childNum == 4) ) {
      it[1] = getParent().getTopLeftFront();	
    }else if( (childNum == 6) ) {
      it[1] = getParent().getTopLeft();	
    }else if( (childNum == 0) ) {
      it[1] = getParent().getLeftFront();	
    }else if( (childNum == 5) ) {
      it[1] = getParent().getTopFront();	
    }else if( (childNum == 7) ){
      it[1] = getParent().getTop();	
    }else if( (childNum == 2) ){
      it[1] = getParent().getLeft();	
    }else if( (childNum == 1) ){
      it[1] = getParent().getFront();	
    }else {
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[4].getTopLeftFront();
    children.clear();
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_TopRightFront() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"TRF"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(3);
    }else if(dim == 2) {
      it.resize(0);
      return it;
    }else {
      it.resize(0);
      return it;
    }
    //5 TRF
    //7 TR but not F
    //1 RF but not T
    //4 TF but not R
    //6 T but not RF
    //3 R but not TF
    //0 F but not TR
    it[0] = getTopRightFront();
    unsigned int childNum = getChildNumber();
    if( (childNum == 5) ) {
      it[1] = getParent().getTopRightFront();	
    }else if( (childNum == 7) ) {
      it[1] = getParent().getTopRight();	
    }else if( (childNum == 1) ) {
      it[1] = getParent().getRightFront();	
    }else if( (childNum == 4) ) {
      it[1] = getParent().getTopFront();	
    }else if( (childNum == 6) ){
      it[1] = getParent().getTop();	
    }else if( (childNum == 3) ){
      it[1] = getParent().getRight();	
    }else if( (childNum == 0) ){
      it[1] = getParent().getFront();	
    }else {
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[5].getTopRightFront();
    children.clear();
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_BottomLeftFront() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"BoLF"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(3);
    }else if(dim == 2) {
      it.resize(0);
      return it;
    }else {
      it.resize(0);
      return it;
    }
    //0 BoLF
    //2 BoL but not F
    //4 LF but not Bo
    //1 BoF but not L
    //3 Bo but not LF
    //6 L but not BoF
    //5 F but not BoL
    it[0] = getBottomLeftFront();
    unsigned int childNum = getChildNumber();
    if( (childNum == 0) ) {
      it[1] = getParent().getBottomLeftFront();	
    }else if( (childNum == 2) ) {
      it[1] = getParent().getBottomLeft();	
    }else if( (childNum == 4) ){
      it[1] = getParent().getLeftFront();	
    }else if( (childNum == 1) ) {
      it[1] = getParent().getBottomFront();	
    }else if( (childNum == 3) ) {
      it[1] = getParent().getBottom();	
    }else if( (childNum == 6) ){
      it[1] = getParent().getLeft();	
    }else if( (childNum == 5) ) {
      it[1] = getParent().getFront();	
    }else {
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[0].getBottomLeftFront();
    children.clear();
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_BottomRightFront() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"BoRF"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(3);
    }else if(dim == 2) {
      it.resize(0);
      return it;
    }else {
      it.resize(0);
      return it;
    }
    //1 BoRF
    //3 BoR
    //5 RF
    //0 BoF
    //2 Bo
    //7 R
    //4 F
    it[0] = getBottomRightFront();
    unsigned int childNum = getChildNumber();
    if( (childNum == 1) ) {
      it[1] = getParent().getBottomRightFront();	
    }else if( (childNum == 3) ) {
      it[1] = getParent().getBottomRight();	
    }else if( (childNum == 5) ) {
      it[1] = getParent().getRightFront();	
    }else if( (childNum == 0) ) {
      it[1] = getParent().getBottomFront();	
    }else if( (childNum == 2) ){
      it[1] = getParent().getBottom();	
    }else if( (childNum == 7) ){
      it[1] = getParent().getRight();	
    }else if( (childNum == 4) ){
      it[1] = getParent().getFront();	
    }else {
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[1].getBottomRightFront();
    children.clear(); 
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_LeftBack() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"LBk"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(4);
    }else if(dim == 2) {
      it.resize(3);
    }else {
      it.resize(0);
      return it;
    }
    //2,6 LBk
    //0,4 L not Bk
    //3,7 Bk not L
    it[0] = getLeftBack();
    unsigned int childNum = getChildNumber();
    if((childNum == 2) || (childNum == 6)) {
      it[1] = getParent().getLeftBack();	
    }else if( (childNum == 0) || (childNum == 4)) {
      it[1] = getParent().getLeft();	
    }else if( (childNum == 3) || (childNum == 7 )) {
      it[1] = getParent().getBack();	
    }else {
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[2].getLeftBack();	
    if(dim == 3) {	
      it[3] = children[6].getLeftBack();
    }
    children.clear(); 
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_RightBack() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"RBk"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(4);
    }else if(dim == 2) {
      it.resize(3);
    }else {
      it.resize(0);
      return it;
    }
    //3,7 RBk
    //1,5 R
    //2,6 Bk
    it[0] = getRightBack();
    unsigned int childNum = getChildNumber();
    if((childNum == 3) || (childNum == 7)) {
      it[1] = getParent().getRightBack();	
    }else if( (childNum == 1) || (childNum == 5) ) {
      it[1] = getParent().getRight();	
    } else if( (childNum == 2) || (childNum == 6)) {
      it[1] = getParent().getBack();	
    }else {
      it[1] = it[0]; 
    }
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[3].getRightBack();	
    if(dim == 3) {
      it[3] = children[7].getRightBack();
    }
    children.clear();
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_TopBack() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"TBk"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(4);
    }else if(dim == 2) {
      it.resize(0);
      return it;
    }else {
      it.resize(0);
      return it;
    }
    //6,7 TBk
    //4,5 T
    //2,3 Bk
    it[0] = getTopBack();
    unsigned int childNum = getChildNumber();
    if( (childNum == 6) || (childNum == 7) ) {
      it[1] = getParent().getTopBack();	
    }else if( (childNum == 4) || (childNum == 5) ){
      it[1] = getParent().getTop();	
    }else if( (childNum == 2) || (childNum == 3) ){
      it[1] = getParent().getBack();	
    }else {
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[6].getTopBack();	
    it[3] = children[7].getTopBack();
    children.clear();
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_BottomBack() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"BoBk"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(4);
    }else if(dim == 2) {
      it.resize(0);
      return it;
    }else {
      it.resize(0);
      return it;
    }
    //2,3 BoBk
    //0,1 Bo
    //6,7 Bk
    it[0] = getBottomBack();
    unsigned int childNum = getChildNumber();
    if( (childNum == 2) || (childNum == 3) ) {
      it[1] = getParent().getBottomBack();	
    }else if( (childNum == 0) || (childNum == 1) ){
      it[1] = getParent().getBottom();	
    }else if( (childNum == 6) || (childNum == 7) ){
      it[1] = getParent().getBack();	
    }else {
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[2].getBottomBack();	 
    it[3] = children[3].getBottomBack();
    children.clear();
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_TopLeftBack() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"TLBk"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(3);
    }else if(dim == 2) {
      it.resize(0);
      return it;
    }else {
      it.resize(0);
      return it;
    }
    //6 TLBk
    //4 TL
    //2 LBk
    //7 TBk
    //5 T
    //0 L
    //3 Bk
    it[0] = getTopLeftBack();
    unsigned int childNum = getChildNumber();
    if( (childNum == 6) ) {
      it[1] = getParent().getTopLeftBack();	
    }else if( (childNum == 4) ) {
      it[1] = getParent().getTopLeft();	
    }else if( (childNum == 2) ) {
      it[1] = getParent().getLeftBack();	
    }else if( (childNum == 7) ) {
      it[1] = getParent().getTopBack();	
    }else if( (childNum == 5) ){
      it[1] = getParent().getTop();	
    }else if( (childNum == 0) ){
      it[1] = getParent().getLeft();	
    }else if( (childNum == 3) ){
      it[1] = getParent().getBack();	
    }else {
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[6].getTopLeftBack();
    children.clear();
    return it;
  }//end fn. 

  inline std::vector<TreeNode >  TreeNode   :: getB_TopRightBack() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"TRBk"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(3);
    }else if(dim == 2) {
      it.resize(0);
      return it;
    }else {
      it.resize(0);
      return it;
    }
    //7 TRBk
    //5 TR
    //3 RBk
    //6 TBk
    //4 T
    //1 R
    //2 Bk
    it[0] = getTopRightBack();
    unsigned int childNum = getChildNumber();
    if( (childNum == 7) ) {
      it[1] = getParent().getTopRightBack();	
    }else if( (childNum == 5) ) {
      it[1] = getParent().getTopRight();	
    }else if( (childNum == 3) ) {
      it[1] = getParent().getRightBack();	
    }else if( (childNum == 6) ) {
      it[1] = getParent().getTopBack();	
    }else if( (childNum == 4) ){
      it[1] = getParent().getTop();	
    }else if( (childNum == 1) ){
      it[1] = getParent().getRight();	
    }else if( (childNum == 2) ){
      it[1] = getParent().getBack();	
    }else {
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[7].getTopRightBack();
    children.clear();
    return it;
  }//end fn. 


  inline std::vector<TreeNode >  TreeNode   :: getB_BottomLeftBack() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"BoLBk"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(3);
    }else if(dim == 2) {
      it.resize(0);
      return it;
    }else {
      it.resize(0);
      return it;
    }
    //2 BoLBk
    //0 BoL
    //6 LBk
    //3 BoBk
    //1 Bo
    //4 L
    //7 Bk
    it[0] = getBottomLeftBack();
    unsigned int childNum = getChildNumber();
    if( (childNum == 2) ) {
      it[1] = getParent().getBottomLeftBack();	
    }else if( (childNum == 0) ) {
      it[1] = getParent().getBottomLeft();	
    }else if( (childNum == 6) ) {
      it[1] = getParent().getLeftBack();	
    }else if( (childNum == 3) ) {
      it[1] = getParent().getBottomBack();	
    }else if( (childNum == 1) ){
      it[1] = getParent().getBottom();	
    }else if( (childNum == 4) ){
      it[1] = getParent().getLeft();	
    }else if( (childNum == 7) ){
      it[1] = getParent().getBack();	
    }else {
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[2].getBottomLeftBack();
    children.clear();
    return it;
  }//end fn. 

  inline std::vector<TreeNode >  TreeNode   :: getB_BottomRightBack() const{
    unsigned int dim = m_uiDim;

#ifdef __DEBUG_TN__
    std::cout<<RED<<"BoRBk"<<NRM<<std::endl;
#endif
    std::vector<TreeNode > it;
    if(dim == 3) {
      it.resize(3);
    }else if(dim == 2) {
      it.resize(0);
      return it;
    }else {
      it.resize(0);
      return it;
    }
    //3 BoRBk
    //1 BoR
    //7 RBk
    //2 BoBk
    //0 Bo
    //5 R
    //6 Bk
    it[0] = getBottomRightBack();
    unsigned int childNum = getChildNumber();
    if( (childNum == 3) ) {
      it[1] = getParent().getBottomRightBack();	
    }else if( (childNum == 1) ) {
      it[1] = getParent().getBottomRight();	
    }else if( (childNum == 7) ) {
      it[1] = getParent().getRightBack();	
    }else if( (childNum == 2) ) {
      it[1] = getParent().getBottomBack();	
    }else if( (childNum == 0) ){
      it[1] = getParent().getBottom();	
    }else if( (childNum == 5) ){
      it[1] = getParent().getRight();	
    }else if( (childNum == 6) ){
      it[1] = getParent().getBack();	
    }else {
      it[1] = it[0]; 
    }	
    std::vector<TreeNode> children;
    addChildren(children);
    it[2] = children[3].getBottomRightBack();
    children.clear();
    return it;
  }//end fn. 

}//end namespace ot

