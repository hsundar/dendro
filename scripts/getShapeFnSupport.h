
#ifndef _GET_SH_FN_SUPPORT_H_
#define _GET_SH_FN_SUPPORT_H_

#include<vector>
#include<iostream>
#include <cstdio>
#include <cmath>

namespace shFn {
  class Point {
    protected:
      double x ,y ,z;
    public:
      Point() : x(0.0), y(0.0), z(0.0) {
      }

      Point(double X, double Y, double Z) : x(X), y(Y), z(Z) {
      }

      //copy constructor	
      Point (const Point & other) {
        this->x = other.x; 
        this->y = other.y; 
        this->z = other.z; 
      }

      //assignment operator
      Point & operator = (Point const  & other) {
        if(this == (&other)) {return *this;}	
        this->x = other.x;
        this->y = other.y; 
        this->z = other.z; 
        return *this;
      }

      //compound assignment
      Point & operator += (const Point & other) {
        if(this == (&other)) {
          Point copy(other);
          this->x += copy.x;        
          this->y += copy.y;        
          this->z += copy.z;        
          return *this;
        }	
        this->x += other.x;        
        this->y += other.y;        
        this->z += other.z;        
        return *this;
      }

      Point & operator -= (const Point & other) {
        if(this == (&other)) {
          Point copy(other);
          this->x -= copy.x;        
          this->y -= copy.y;        
          this->z -= copy.z;        
          return *this;
        }	
        this->x -= other.x;        
        this->y -= other.y;        
        this->z -= other.z;        
        return *this;
      }

      Point & operator *= (double fac) {
        this->x *= fac;        
        this->y *= fac;        
        this->z *= fac;        
        return *this;
      }

      //binary addition
      const Point operator + (const Point & other) const {
        Point result = (*this);
        result += other;
        return result;      
      }

      //binary subtraction
      const Point operator - (const Point & other) const {
        Point result = (*this);
        result -= other;
        return result;      
      }

      const Point operator * (double fac) const {
        Point result = (*this);
        result *= fac;
        return result;
      }

      bool operator == ( Point const  &other) const {
        return ( (this->x == other.x) && (this->y == other.y) && (this->z == other.z) );
      }

      //Precendence: x > y > z
      bool operator < ( Point const & other ) const {
        if( this->x != other.x ) { 
          return ( this->x < other.x );
        } else if ( this->y != other.y ) {
          return ( this->y < other.y );
        } else {
          return ( this->z < other.z );
        }		
      }

      bool operator >= ( Point const & other ) const {
        return (!((*this) < other));
      }

      bool operator <= ( Point const & other ) const {
        return ( ((*this) < other) || ((*this) == other) );
      }

      bool operator > ( Point const & other) const {
        return (!((*this) <= other));
      }

      bool operator != ( Point const  &other) const {
        return (!((*this) == other));
      }

      double X() const {
        return this->x;
      }

      double Y() const {
        return this->y;
      }

      double Z() const {
        return this->z;
      }

      friend std::ostream & operator << (std::ostream & os, Point const & p);

  };

  class Box {
    protected:
      Point anchor;
      double size; 
    public:
      Box() : anchor(), size(0.0) {
      }

      Box(Point pt, double sz) : anchor(pt), size(sz) {
      }

      //copy constructor	
      Box (const Box & other) {
        this->anchor  = other.anchor;
        this->size  = other.size;
      }//end function

      //assignment operator
      Box & operator = (Box const  & other) {
        if(this == (&other)) {return *this;}	
        this->anchor = other.anchor;
        this->size = other.size;
        return *this;
      }//end fn.


      //compound assignment
      Box & operator += (const Point & disp) {
        this->anchor += disp;
        return *this;
      }

      Box & operator -= (const Point & disp) {
        this->anchor -= disp;        
        return *this;
      }

      Box & operator *= (double fac) {
        this->anchor *= fac;
        this->size *= fac;
        return *this;
      }

      //binary addition
      const Box operator + (const Point & disp) const {
        Box result = (*this);
        result += disp;
        return result;      
      }

      //binary subtraction
      const Box operator - (const Point & disp) const {
        Box result = (*this);
        result -= disp;
        return result;      
      }

      const Box operator * (double fac) const{
        Box result = (*this);
        result *= fac;
        return result;
      }

      bool  operator == ( Box const  &other) const {
        return ( (this->size == other.size) && (this->anchor == other.anchor) );
      }

      bool  operator != (Box const  &other) const {
        return (!((*this) == other));
      }

      const Point Anchor() const {
        return this->anchor;
      }

      double Size() const {
        return this->size;
      }

      const Point Vertex(unsigned int i) const {
        Point vtx;
        switch(i) {
          case 0 : {
                     Point disp(0,0,0);
                     vtx = (this->anchor) + disp;
                     break;
                   }
          case 1 : {
                     Point disp(this->size,0,0);
                     vtx = (this->anchor) + disp;
                     break;
                   }
          case 2 : {
                     Point disp(0,this->size,0);
                     vtx = (this->anchor) + disp;
                     break;
                   }
          case 3 : {
                     Point disp(this->size,this->size,0);
                     vtx = (this->anchor) + disp;
                     break;
                   }
          case 4 : {
                     Point disp(0,0,this->size);
                     vtx = (this->anchor) + disp;
                     break;
                   }
          case 5 : {
                     Point disp(this->size,0,this->size);
                     vtx = (this->anchor) + disp;
                     break;
                   }
          case 6 : {
                     Point disp(0,this->size,this->size);
                     vtx = (this->anchor) + disp;
                     break;
                   }
          case 7 : {
                     Point disp(this->size,this->size,this->size);
                     vtx = (this->anchor) + disp;
                     break;
                   }
          default : {
                      assert(false);
                    }
        }
        return vtx;
      }

      bool Intersects(const Box & other) const {
        bool result = false;
        double minX1 = (this->anchor).X();
        double minY1 = (this->anchor).Y();
        double minZ1 = (this->anchor).Z();

        Point myLastVtx = this->Vertex(7);

        double maxX1 = (myLastVtx).X();
        double maxY1 = (myLastVtx).Y();
        double maxZ1 = (myLastVtx).Z();

        double minX2 = (other.anchor).X();
        double minY2 = (other.anchor).Y();
        double minZ2 = (other.anchor).Z();

        Point othLastVtx = other.Vertex(7);

        double maxX2 = (othLastVtx).X();
        double maxY2 = (othLastVtx).Y();
        double maxZ2 = (othLastVtx).Z();

        if( (minX1 >= minX2) && (minX1 < maxX2) &&
            (minY1 >= minY2) && (minY1 < maxY2) && 
            (minZ1 >= minZ2) && (minZ1 < maxZ2) ) {
          result = true;
        }else if( (minX2 >= minX1) && (minX2 < maxX1) &&
            (minY2 >= minY1) && (minY2 < maxY1) &&
            (minZ2 >= minZ1) && (minZ2 < maxZ1) ) {     
          result = true;
        }

        return result;
      }

      virtual std::ostream & PrintMe(std::ostream & os) const;
      friend std::ostream & operator << (std::ostream & os, Box const & e);
  };

  class Node : public Point {
    private:
      bool isHanging;
    public:
      //constructors
      Node() : Point(), isHanging(false) { }
      Node(bool hn) : Point(), isHanging(hn) {}
      Node(Point pt, bool hn) : Point(pt), isHanging(hn) {}

      //copy constructor	
      Node (const Node & other) : Point(other) {
        this->isHanging = other.isHanging;
      }

      //assignment operator
      Node & operator = (Node const  & other) {
        if(this == (&other)) {return *this;}	
        this->x = other.x;
        this->y = other.y; 
        this->z = other.z; 
        this->isHanging = other.isHanging;
        return *this;
      }//end fn.

      bool  Equals ( Node const  &other) const {
        return ( (this->x == other.x) &&
            (this->y == other.y) &&
            (this->z == other.z) &&
            (this->isHanging == other.isHanging) );
      }

      bool IsHanging() {
        return this->isHanging;
      }

      const Point IgnoreHanging() {
        return (*this);
      }
  };

  class Octant : public Box {
    private:
      unsigned int childNumber, hangingType;
    public:

      Octant() : Box(), childNumber(0), hangingType(0) {
      }

      Octant(unsigned int cNum, unsigned int hnType) : Box() {
        this->childNumber = cNum;
        this->hangingType = hnType;
      }

      Octant(const Box & b, unsigned int cNum,
          unsigned int hnType) : Box(b) {
        this->childNumber = cNum;
        this->hangingType = hnType;
      }

      Octant(Point pt, double sz, unsigned int cNum, 
          unsigned int hnType) : Box(pt, sz) {
        this->childNumber = cNum;
        this->hangingType = hnType;
      }

      //copy constructor	
      Octant (const Octant & other) : Box(other) {
        this->childNumber = other.childNumber;
        this->hangingType = other.hangingType;
      }//end function

      //assignment operator
      Octant & operator = (Octant const  & other) {
        if(this == (&other)) {return *this;}	
        this->anchor = other.anchor;
        this->size = other.size;
        this->childNumber = other.childNumber;
        this->hangingType = other.hangingType;
        return *this;
      }//end fn.

      bool  operator == ( Octant const  &other) const {
        return ( (this->size == other.size) &&
            (this->anchor == other.anchor) &&
            (this->childNumber == other.childNumber) &&
            (this->hangingType == other.hangingType) );
      }

      bool  operator != (Octant const  &other) const {
        return (!((*this) == other));
      }

      bool EqualsIgnoreHanging (Octant const & other) const {
        return ( (this->size == other.size) &&
            (this->anchor == other.anchor) &&
            (this->childNumber == other.childNumber) );
      }

      bool BoxEquals (Octant const & other) const {
        return ( (this->size == other.size) &&
            (this->anchor == other.anchor) );
      }

      void SetChildNumber(unsigned int cNum) {
        this->childNumber = cNum;
      }

      void SetHangingType(unsigned int hnType) {
        this->hangingType = hnType;
      }

      unsigned int ChildNumber() const {
        return this->childNumber;
      }

      unsigned int HangingType() const {
        return this->hangingType;
      }

      const Box Parent() const {
        assert(this->childNumber < 8);
        Box myParent = ((*this)*2.0);
        Point offset = (myParent.Vertex(this->childNumber) - (this->Vertex(this->childNumber)));
        myParent -= offset;
        return myParent;
      }

      const Octant Brother(unsigned int cNumOfBrother) const {
        assert(this->childNumber != cNumOfBrother);
        assert(this->childNumber < 8);
        assert(cNumOfBrother < 8);
        Box myParent = this->Parent();
        Octant brother = (*this);
        brother.childNumber = cNumOfBrother;
        Point offset = ((brother.Vertex(cNumOfBrother)) - (myParent.Vertex(cNumOfBrother)));
        brother -= offset;
        return brother; 
      }

      virtual std::ostream & PrintMe(std::ostream & os) const ;
  };

  template <typename T>
    bool BinarySearch(const T* arr, unsigned int nelem,  T & key, unsigned int *ret_idx) {
      if(!nelem) {*ret_idx = nelem; return false;}
      unsigned int left = 0;
      unsigned int right = (nelem -1);	
      while (left <= right) {
        unsigned int mid = (unsigned int)( left + (unsigned int)(floor((double)(right-left)/2.0)) );
        if (key > arr[mid]) {
          left  = mid+1;
        } else if (key < arr[mid]) {
          if(mid>0) { right = mid-1; }
          else { right = 0; break;}
        } else {
          *ret_idx = mid;
          return true;
        }//end if-else-if
      }//end while
      *ret_idx = nelem;	
      return false;
    }//end function

  //Expects a sorted unique array.
  template <typename T>
    bool maxLowerBound(const std::vector<T>& arr, const T & key, unsigned int &ret_idx,
        unsigned int* leftIdx, unsigned int* rightIdx ) {
      unsigned int nelem = arr.size();
      ret_idx = 0;
      if(!nelem) { return false;}
      if(arr[0] > key) { return false;}
      if(arr[nelem-1] < key) {
        ret_idx = (nelem-1);
        return true;
      }//end if	

      //binary search
      unsigned int left = 0;
      unsigned int right = (nelem -1);	
      unsigned int mid = 0;
      if(leftIdx) {
        left = (*leftIdx);
      }
      if(rightIdx) {
        right = (*rightIdx);
      }
      while (left <= right) {
        mid = (unsigned int)( left + (unsigned int)(floor((double)(right-left)/2.0)) );
        if (key > arr[mid]) {
          left  = mid + ((unsigned int)1);
        } else if (key < arr[mid]){
          if(mid>0) {
            right = mid-1;
          }else {
            right=0;
            break;
          }
        } else {
          ret_idx = mid;
          return true;
        }//end if-else-if
      }//end while

      //If binary search did not find an exact match, it would have
      //stopped one element after or one element before. 

      if( (arr[mid] > key) && (mid > 0) ){ mid--; }	
      if(arr[mid] <= key ) { ret_idx = mid; return true; }
      else { ret_idx = 0; return false;}
    }//end function

  template<typename T> void makeVectorUnique(std::vector<T>& vecT, bool isSorted) {
    if(vecT.size() < 2) { return;}
    if(!isSorted) { sort(vecT.begin(),vecT.end()); }
    std::vector<T> tmp(vecT.size());
    tmp[0] = vecT[0];
    unsigned int tmpSize=1;
    for(unsigned int i=1;i<vecT.size();i++) {	
      if(tmp[tmpSize-1] != vecT[i]) {
        tmp[tmpSize] = vecT[i];
        tmpSize++;
      }
    }//end for
    vecT.clear();
    vecT.insert(vecT.begin(), tmp.begin(), (tmp.begin() + tmpSize));
    tmp.clear();
  }//end function

  void GetCnumBasedHangingMasks(unsigned int cNum, std::vector<unsigned int>& hnMasks);

  void GetNeighboringCnumAndSizes(std::vector<std::vector<
      std::vector<shFn::Octant> > > & nhCnumAndSz);

  int readOctantsFromFile(char* filename, std::vector<shFn::Octant> & octs);
  int readOctantsFromFile(char* filename, std::vector<std::vector<shFn::Octant> >& octs);

  int writeOctantsToFile(char* filename, const std::vector<shFn::Octant> & octs);
  int writeOctantsToFile(char* filename, const std::vector<std::vector<shFn::Octant> >& octs);

  unsigned int ListAllElemsSharing0(const std::vector<std::vector<std::vector<shFn::Octant> > > & nhCnumAndSz);

  unsigned int EliminateOverlaps(unsigned int numElemFilesWritten);

  bool FoundDirectOverlap(const std::vector<shFn::Octant> & list);

  bool FoundInDirectOverlap(const std::vector<shFn::Octant> & list);

  unsigned int IncludeHangingTypes(unsigned int numValidElemFilesWritten);

  void InsertNodesAndOctant(unsigned int hnType, const std::vector<shFn::Node>& tmpNodes,
      std::vector<shFn::Node>& nodesVisited, shFn::Octant& currOct,
      std::vector<shFn::Octant>& currRecord );

  bool IsConfigValid(unsigned int hnMask, const shFn::Octant& elem,
      const shFn::Box& parent, const std::vector<shFn::Node>& nodesVisited,
      std::vector<shFn::Node>& tmpNodes);

  void AddPossibleHnCombinationsPrivate( const std::vector<std::vector<unsigned int> >& potentialHnMasks,
      unsigned int elemId, const std::vector<shFn::Octant> & elements,
      const std::vector<shFn::Node>& nodesVisitedIn,
      const std::vector<shFn::Octant> & currRecordIn,
      std::vector<std::vector<shFn::Octant> >& elemListWithHn );

  void AddPossibleHnCombinations( const std::vector<std::vector<unsigned int> >& potentialHnMasks,
      const std::vector<shFn::Octant> & elements,
      std::vector<std::vector<shFn::Octant> >& elemListWithHn );

  void AddPossibleHnCombinationsForAux( const std::vector<std::vector<unsigned int> >& potentialHnMasks,
      const std::vector<shFn::Octant> & elements,
      std::vector<std::vector<shFn::Octant> >& elemListWithHn );

  unsigned int AddAuxilaryElems(unsigned int numOctFilesWritten);

  unsigned int FinalElemList(unsigned int numAuxFilesWritten);

  void FinalVerticesList(unsigned int numFinalElemFilesWritten, 
      unsigned int & minCnt, unsigned int & maxCnt);

  void RunAllInCore(const std::vector<std::vector<
      std::vector<shFn::Octant> > > & nhCnumAndSz,
      unsigned int & minCnt, unsigned int & maxCnt);

  void EliminateOverlapsInCore(const std::vector<
      std::vector<shFn::Octant> > & elemList,
      unsigned int & minCnt, unsigned int & maxCnt);

  void IncludeHangingTypesInCore(const std::vector<
      std::vector<shFn::Octant> > & validElemList,
      unsigned int & minCnt, unsigned int & maxCnt);

  void AddAuxilaryElemsInCore(const std::vector<
      std::vector<shFn::Octant> > & elemListWithHn,
      unsigned int & minCnt, unsigned int & maxCnt); 

  void FinalElemListInCore(const std::vector<
      std::vector<shFn::Octant> > & elemListWithAux,
      unsigned int & minCnt, unsigned int & maxCnt);

  void FinalVerticesListInCore(const std::vector<
      std::vector<shFn::Octant> > & finalElemList,
      unsigned int & minCnt, unsigned int & maxCnt);

}//end namespace

#endif 


