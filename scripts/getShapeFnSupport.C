#include<mpi.h> 
#include<vector>
#include<iostream>
#include <cstdio>
#include "getShapeFnSupport.h"

namespace shFn {

  std::ostream & operator <<(std::ostream & os, shFn::Point const & p){
    return (os << p.X() <<", "<< p.Y()<<", "<< p.Z());
  }//end fn.

  std::ostream & operator <<(std::ostream & os, shFn::Box const & e){
    return e.PrintMe(os);
  }//end fn.

  std::ostream & Box::PrintMe(std::ostream & os) const{
    return (os << this->anchor <<" @ "<< this->size);
  }//end fn.

  std::ostream & Octant::PrintMe(std::ostream & os) const {
    return (os << this->anchor <<" @ "<< this->size <<
        ": ("<<this->childNumber<<", "<<this->hangingType<<")" );
  }//end fn.

  int readOctantsFromFile(char* filename, std::vector<shFn::Octant> & octs) {
    FILE* infile = fopen(filename,"r");
    unsigned int numOcts;
    fscanf(infile,"%u",&numOcts);
    octs.resize(numOcts);   
    for (unsigned int i =0; i< numOcts; i++) {
      unsigned int cnum,hn;
      double x,y,z,sz;
      fscanf(infile,"%lf",&x);
      fscanf(infile,"%lf",&y);
      fscanf(infile,"%lf",&z);
      fscanf(infile,"%lf",&sz); 
      fscanf(infile,"%u",&cnum);
      fscanf(infile,"%u",&hn);
      octs[i] = shFn::Octant (Point(x,y,z),sz, cnum, hn);
    }//end for i
    fclose(infile);
    return 1;
  }//end function

  //Overloaded for vec of vec
  int readOctantsFromFile(char* filename, std::vector<std::vector<shFn::Octant> >& octs) {
    FILE* infile = fopen(filename,"r");
    unsigned int numOcts;
    fscanf(infile,"%u",&numOcts);
    octs.resize(numOcts);   
    for (unsigned int i = 0; i < numOcts; i++) {
      unsigned int numOcts_i;
      fscanf(infile,"%u",&numOcts_i);
      octs[i].resize(numOcts_i);  
      for(unsigned int j = 0; j < numOcts_i; j++) {
        unsigned int cnum,hn;
        double x,y,z,sz;
        fscanf(infile,"%lf",&x);
        fscanf(infile,"%lf",&y);
        fscanf(infile,"%lf",&z);
        fscanf(infile,"%lf",&sz); 
        fscanf(infile,"%u",&cnum);
        fscanf(infile,"%u",&hn);
        octs[i][j] = shFn::Octant(Point(x,y,z),sz, cnum, hn);
      }//end for j
    }//end for i
    fclose(infile);
    return 1;
  }//end function

  int writeOctantsToFile(char* filename, const std::vector<shFn::Octant> & octs) {
    FILE* outfile = fopen(filename,"w");
    unsigned int numOcts = octs.size();
    fprintf(outfile,"%u\n",numOcts);
    for(unsigned int i = 0; i < numOcts; i++) {
      shFn::Point currAnchor = octs[i].Anchor();
      fprintf(outfile,"%lf %lf %lf %lf %u %u\n",
          currAnchor.X(), currAnchor.Y(), currAnchor.Z(),
          octs[i].Size(),octs[i].ChildNumber(),octs[i].HangingType());
    }//end for i 
    fprintf(outfile,"\n");
    fclose(outfile);
    return 1;
  }//end function

  //Overloaded for vec of vec
  int writeOctantsToFile(char* filename, const std::vector<std::vector<shFn::Octant> >& octs) {
    FILE* outfile = fopen(filename,"w");
    unsigned int numOcts = octs.size();
    fprintf(outfile,"%u\n\n",numOcts);
    for(unsigned int i = 0; i < numOcts; i++) {
      unsigned int numOcts_i = octs[i].size();
      fprintf(outfile,"%u\n",numOcts_i);
      for(unsigned int j = 0; j < numOcts_i; j++) {
        shFn::Point currAnchor = octs[i][j].Anchor();
        fprintf(outfile,"%lf %lf %lf %lf %u %u\n",
            currAnchor.X(), currAnchor.Y(), currAnchor.Z(),
            octs[i][j].Size(),octs[i][j].ChildNumber(),octs[i][j].HangingType());
      }//end for j
      fprintf(outfile,"\n");
    }//end for i 
    fclose(outfile);
    return 1;
  }//end function

  void GetCnumBasedHangingMasks(unsigned int cNum, std::vector<unsigned int>& hnMasks) {
    //There are only 18 valid hanging types depending on the childnumber.
    //Values in the same order as in the C++ code (oda.h), if node i (0-based
    //indexing) is hanging then 1u<<i is set to 1.

    hnMasks.resize(18);
    switch(cNum) {
      case 0: {
                hnMasks[0]  =  0;
                hnMasks[1]  =  4;
                hnMasks[2]  =  2;
                hnMasks[3]  =  6;
                hnMasks[4]  = 16;
                hnMasks[5]  = 20;
                hnMasks[6]  = 18;
                hnMasks[7]  = 22;
                hnMasks[8]  = 14;
                hnMasks[9]  = 30;
                hnMasks[10] = 84;
                hnMasks[11] = 86;
                hnMasks[12] = 94;
                hnMasks[13] = 50;
                hnMasks[14] = 54;
                hnMasks[15] = 62;
                hnMasks[16] = 118;
                hnMasks[17] = 126;
                break;
              }
      case 1: {
                hnMasks[0]  =  0;
                hnMasks[1]  =  1;
                hnMasks[2]  =  8;
                hnMasks[3]  =  9;
                hnMasks[4]  =  32;
                hnMasks[5]  =  33;
                hnMasks[6]  =  40;
                hnMasks[7]  =  41;
                hnMasks[8]  =  13;
                hnMasks[9]  =  45;
                hnMasks[10] =  49;
                hnMasks[11] =  57;
                hnMasks[12] =  61;
                hnMasks[13] =  168;
                hnMasks[14] =  169;
                hnMasks[15] =  173;
                hnMasks[16] =  185;
                hnMasks[17] =  189;
                break;
              }
      case 2: {
                hnMasks[0]  = 0;
                hnMasks[1]  = 8;
                hnMasks[2]  = 1;
                hnMasks[3]  = 9;
                hnMasks[4]  = 64;
                hnMasks[5]  = 72;
                hnMasks[6]  = 65;
                hnMasks[7]  = 73;
                hnMasks[8]  = 11;
                hnMasks[9]  = 75;
                hnMasks[10] = 200;
                hnMasks[11] = 201;
                hnMasks[12] = 203;
                hnMasks[13] = 81;
                hnMasks[14] = 89;
                hnMasks[15] = 91;
                hnMasks[16] = 217;
                hnMasks[17] = 219;
                break;
              }
      case 3: {
                hnMasks[0]  = 0;
                hnMasks[1]  = 2;
                hnMasks[2]  = 4;
                hnMasks[3]  = 6;
                hnMasks[4]  = 128;
                hnMasks[5]  = 130;
                hnMasks[6]  = 132;
                hnMasks[7]  = 134;
                hnMasks[8]  = 7;
                hnMasks[9]  = 135;
                hnMasks[10] = 162;
                hnMasks[11] = 166;
                hnMasks[12] = 167;
                hnMasks[13] = 196;
                hnMasks[14] = 198;
                hnMasks[15] = 199;
                hnMasks[16] = 230;
                hnMasks[17] = 231;
                break;
              }
      case 4: {
                hnMasks[0]  = 0;
                hnMasks[1]  = 1;
                hnMasks[2]  = 32;
                hnMasks[3]  = 33;
                hnMasks[4]  = 64;
                hnMasks[5]  = 65;
                hnMasks[6]  = 96;
                hnMasks[7]  = 97;
                hnMasks[8]  = 35; 
                hnMasks[9]  = 99;
                hnMasks[10] = 69;
                hnMasks[11] = 101;
                hnMasks[12] = 103;
                hnMasks[13] = 224;
                hnMasks[14] = 225;
                hnMasks[15] = 227;
                hnMasks[16] = 229;
                hnMasks[17] = 231;
                break;
              }
      case 5: {
                hnMasks[0]  = 0;
                hnMasks[1]  = 2;
                hnMasks[2]  = 128;
                hnMasks[3]  = 130;
                hnMasks[4]  = 16;
                hnMasks[5]  = 18;
                hnMasks[6]  = 144;
                hnMasks[7]  = 146;
                hnMasks[8]  = 138;
                hnMasks[9]  = 154;
                hnMasks[10] = 19;
                hnMasks[11] = 147;
                hnMasks[12] = 155;
                hnMasks[13] = 208;
                hnMasks[14] = 210;
                hnMasks[15] = 218;
                hnMasks[16] = 211;
                hnMasks[17] = 219;
                break;
              }
      case 6: {
                hnMasks[0]  = 0;
                hnMasks[1]  = 4;
                hnMasks[2]  = 16;
                hnMasks[3]  = 20;
                hnMasks[4]  = 128;
                hnMasks[5]  = 132;
                hnMasks[6]  = 144;
                hnMasks[7]  = 148;
                hnMasks[8]  = 21;
                hnMasks[9]  = 149;
                hnMasks[10] = 140;
                hnMasks[11] = 156;
                hnMasks[12] = 157;
                hnMasks[13] = 176;
                hnMasks[14] = 180;
                hnMasks[15] = 181;
                hnMasks[16] = 188;
                hnMasks[17] = 189;
                break;
              }
      case 7: {
                hnMasks[0]  = 0;
                hnMasks[1]  = 8;
                hnMasks[2]  = 64;
                hnMasks[3]  = 72;
                hnMasks[4]  = 32;
                hnMasks[5]  = 40;
                hnMasks[6]  = 96;
                hnMasks[7]  = 104;
                hnMasks[8]  = 76;
                hnMasks[9]  = 108;
                hnMasks[10] = 42;
                hnMasks[11] = 106;
                hnMasks[12] = 110;
                hnMasks[13] = 112;
                hnMasks[14] = 120;
                hnMasks[15] = 124;
                hnMasks[16] = 122;
                hnMasks[17] = 126;
                break;
              }
      default: {
                 assert(false);
               }
    }
  }//end fn.

  void GetNeighboringCnumAndSizes(std::vector<std::vector<
      std::vector<shFn::Octant> > > & nhCnumAndSz) {

    //Notation: element-i is the element whose i-th vertex is the
    //grid point of interest.
    //For simplicity, element-0 is chosen as the standard reference element
    //[0,1]x[0,1]x[0,1] and the grid point of interest has the coords (0,0,0).

    //A generic element [0,1]x[0,1]x[0,1]

    nhCnumAndSz.resize(8);
    for(unsigned int i = 0; i < 8; i++) {
      nhCnumAndSz[i].resize(7);
    }//end for i

    // Each of 7 surrounding elements could either be of the same size (szOpt =0) or half
    // the size (szOpt = 1) or double the size (szOpt = 2) of the reference element.
    for(unsigned int myCnum = 0; myCnum < 8; myCnum++) { 
      shFn::Octant coordsTemplate(shFn::Point(),1.0,myCnum,0);
      for(unsigned int elemNum = 1; elemNum < 8; elemNum++) {
        for(unsigned int szOpt = 0; szOpt < 3; szOpt++) {
          shFn::Octant elemCoords = coordsTemplate;
          if(szOpt == 1) {
            elemCoords *= 0.5;
          } else if(szOpt == 2) {
            elemCoords *= 2.0;
          }//end if-else
          Point offset = elemCoords.Vertex(elemNum);
          elemCoords -= offset;
          for(unsigned int othCnum = 0; othCnum < 8; othCnum++) {
            elemCoords.SetChildNumber(othCnum);
            bool validBox = true;
            for(unsigned int cNumOfBrother = 0; cNumOfBrother < othCnum; cNumOfBrother++) {
              shFn::Octant brotherCoords = elemCoords.Brother(cNumOfBrother);
              if( coordsTemplate == brotherCoords ) {  
                //Brothers are always valid.
                break;
              }
              bool weIntersect = brotherCoords.Intersects(coordsTemplate);
              if(weIntersect) {
                validBox = false;
                break;
              }
            }//end for cNumBro
            if(validBox) {
              for(unsigned int cNumOfBrother = (othCnum +1);
                  cNumOfBrother < 8; cNumOfBrother++) {
                shFn::Octant brotherCoords = elemCoords.Brother(cNumOfBrother);
                if( coordsTemplate == brotherCoords ) {
                  break;
                }
                bool weIntersect =  brotherCoords.Intersects(coordsTemplate);
                if(weIntersect) {
                  validBox = false;
                  break;
                }
              }//end for cNumBro
              if(validBox) {
                nhCnumAndSz[myCnum][elemNum-1].push_back(elemCoords);
              }
            }
          }//end for othCnum 
        }//end for szOpt
      }//end for elemNum
    }//end for myCnum
  }//end fn.

  unsigned int ListAllElemsSharing0(const std::vector<std::vector<
      std::vector<shFn::Octant> > > & nhCnumAndSz ) {

    shFn::Box coordsTemplate(shFn::Point(),1.0);
    unsigned int elemFileCtr = 0;
    std::vector<std::vector<shFn::Octant> > elemList;

    //First list all possible sets of elements to vist.
    //This is basically selecting one from each elemNum type
    //Each element has the following info: anchor, cNum, sz
    for(unsigned int myCnum = 0; myCnum < 8; myCnum++) {
      std::cout<<"Generating elemList for cNum: "<<myCnum<<std::endl;
      unsigned int lenNh[7];
      for(unsigned int i = 0; i < 7; i++) {
        lenNh[i] = nhCnumAndSz[myCnum][i].size();
        std::cout<<" Len("<<i<<"): "<<lenNh[i]<<std::endl;
      }
      unsigned int elemNumIndices[7];
      for(unsigned int e1Idx = 0; e1Idx < lenNh[0]; e1Idx++) {	
        elemNumIndices[0] = e1Idx;
        for(unsigned int e2Idx = 0; e2Idx < lenNh[1]; e2Idx++) {
          elemNumIndices[1] = e2Idx;
          for(unsigned int e3Idx = 0; e3Idx < lenNh[2]; e3Idx++) {
            elemNumIndices[2] = e3Idx;
            for(unsigned int e4Idx = 0; e4Idx < lenNh[3]; e4Idx++) {
              elemNumIndices[3] = e4Idx;
              for(unsigned int e5Idx = 0; e5Idx < lenNh[4]; e5Idx++) {
                elemNumIndices[4] = e5Idx;
                for(unsigned int e6Idx = 0; e6Idx < lenNh[5]; e6Idx++) {
                  elemNumIndices[5] = e6Idx;
                  for(unsigned int e7Idx = 0; e7Idx < lenNh[6]; e7Idx++) {
                    elemNumIndices[6] = e7Idx;
                    //New Set
                    std::vector<shFn::Octant> record;
                    record.push_back(shFn::Octant(coordsTemplate,myCnum,0));
                    for(unsigned int elemNum = 1; elemNum <8; elemNum++) {
                      shFn::Octant currOct = 
                        nhCnumAndSz[myCnum][elemNum-1][elemNumIndices[elemNum-1]];
                      assert(currOct.Vertex(elemNum) == Point());
                      assert(currOct.HangingType() == 0);
                      record.push_back(currOct);
                    }
                    elemList.push_back(record);
                    if(elemList.size() > 10000) {
                      char fname[256];
                      sprintf(fname,"elemFile_%u.txt",elemFileCtr);
                      writeOctantsToFile(fname, elemList);   
                      elemList.clear();						
                      elemFileCtr++;
                    }
                  }//e7
                }//e6
              }//e5	
            }//e4
          }//e3		
        }//e2
      }//e1       	
    }//myCnum

    if(elemList.size() > 0) {
      char fname[256];
      sprintf(fname,"elemFile_%u.txt",elemFileCtr);
      writeOctantsToFile(fname, elemList);   
      elemList.clear();						
      elemFileCtr++;
    }

    return elemFileCtr;
  }//end fn.

  unsigned int EliminateOverlaps(unsigned int numElemFilesWritten) {
    unsigned int validElemFileCtr = 0;
    std::vector<std::vector<shFn::Octant> > validElemList;

    for(unsigned int elemFileCtr = 0; 
        elemFileCtr < numElemFilesWritten; elemFileCtr++) {
      std::vector<std::vector<shFn::Octant> > elemList;
      char fname[256];
      sprintf(fname,"elemFile_%u.txt",elemFileCtr);
      readOctantsFromFile(fname, elemList);   

      for(unsigned int i = 0; i < elemList.size(); i++) {
        bool foundDirectOverlap = FoundDirectOverlap(elemList[i]);
        if(!foundDirectOverlap) {
          bool foundInDirectOverlap = FoundInDirectOverlap(elemList[i]);
          if(!foundInDirectOverlap) {
            validElemList.push_back(elemList[i]);
          }
        }
      }

      elemList.clear();     

      if(validElemList.size() > 5000) {
        sprintf(fname,"validElemFile_%u.txt",validElemFileCtr);
        writeOctantsToFile(fname, validElemList);   
        validElemList.clear();						
        validElemFileCtr++;
      }
    }

    if(validElemList.size() > 0) {
      char fname[256];
      sprintf(fname,"validElemFile_%u.txt",validElemFileCtr);
      writeOctantsToFile(fname, validElemList);   
      validElemList.clear();						
      validElemFileCtr++;
    }

    return validElemFileCtr;
  }//end fn

  bool FoundDirectOverlap(const std::vector<shFn::Octant> & list) {
    static bool firstOverlapCase = true;

    bool result = false;
    unsigned int len = list.size();
    for(unsigned int i = 0; i < len; i++) {
      for(unsigned int j = i+1; j < len; j++) {
        bool weIntersect = list[i].Intersects(list[j]);
        if(weIntersect) {
          result = true;
          if(firstOverlapCase) {
            std::cout<<"Direct Overlap between: "<<std::endl
              <<list[i]<<std::endl
              <<" and "<<std::endl
              <<list[j]<<std::endl;
            firstOverlapCase = false;
          }
          break;
        }
      }
    }    
    return result;
  }//end fn.

  bool FoundInDirectOverlap(const std::vector<shFn::Octant> & list) {
    static bool firstOverlapCase = true;
    bool result = false;
    unsigned int len = list.size();
    for(unsigned int i = 0; i < len; i++) {
      for(unsigned int j = i+1; j < len; j++) {
        for(unsigned int cNumOfBroOfI = 0;
            cNumOfBroOfI < 8; cNumOfBroOfI++) {
          if(cNumOfBroOfI == list[i].ChildNumber()) {
            continue;
          }
          shFn::Octant broI = list[i].Brother(cNumOfBroOfI);
          if(broI != list[j]) {
            bool weIntersect = broI.Intersects(list[j]);
            if(weIntersect) {
              result = true;
              if(firstOverlapCase) {
                std::cout<<"Indirect Overlap between: "<<std::endl
                  <<list[i]<<std::endl
                  <<" and "<<std::endl
                  <<list[j]<<std::endl;
                firstOverlapCase = false;
              }
              break;
            }
          }
        }
        if(!result) {
          for(unsigned int cNumOfBroOfJ = 0;
              cNumOfBroOfJ < 8; cNumOfBroOfJ++) {
            if(cNumOfBroOfJ == list[j].ChildNumber()) {
              continue;
            }
            shFn::Octant broJ = list[j].Brother(cNumOfBroOfJ);
            if(broJ != list[i]) {
              bool weIntersect = broJ.Intersects(list[i]);
              if(weIntersect) {
                result = true;
                if(firstOverlapCase) {
                  std::cout<<"Indirect Overlap between: "<<std::endl
                    <<list[i]<<std::endl
                    <<" and "<<std::endl
                    <<list[j]<<std::endl;
                  firstOverlapCase = false;
                }
                break;
              }
            }
          }
        }
      }
    }    
    return result;
  }//end fn.

  unsigned int IncludeHangingTypes(unsigned int numValidElemFilesWritten) {
    //Until now all the hanging flags are initialized to 0.
    //Now we create a list of all the possible elements along
    //with their respective hanging types.
    //To do this basically keep track of the vertices you have
    //visited and their hanging flag. [0,0,0] must never be hanging.

    unsigned int octFileCtr = 0;
    std::vector<std::vector<shFn::Octant> > elemListWithHn;

    std::vector<std::vector<unsigned int> > potentialHnMasks(8);
    for(unsigned int i = 0; i< 8; i++) {
      GetCnumBasedHangingMasks( i, potentialHnMasks[i]);  
    }

    for(unsigned int validElemFileCtr = 0; 
        validElemFileCtr < numValidElemFilesWritten;
        validElemFileCtr++) {
      char fname[256];
      std::vector<std::vector<shFn::Octant> > validElemList;
      sprintf(fname,"validElemFile_%u.txt",validElemFileCtr);
      readOctantsFromFile(fname, validElemList);   

      for(unsigned int i = 0; i < validElemList.size(); i++) {
        //Recursion...
        AddPossibleHnCombinations(potentialHnMasks, validElemList[i], elemListWithHn);
        if(elemListWithHn.size() > 5000) {
          sprintf(fname,"octsWithHn_%u.txt",octFileCtr);
          writeOctantsToFile(fname, elemListWithHn);   
          elemListWithHn.clear();						
          octFileCtr++;
        }
      }//end for i

      validElemList.clear();
    }

    if(elemListWithHn.size() > 0) {
      char fname[256];
      sprintf(fname,"octsWithHn_%u.txt",octFileCtr);
      writeOctantsToFile(fname, elemListWithHn);   
      elemListWithHn.clear();						
      octFileCtr++;
    }

    return octFileCtr;    
  }//end fn.

  void AddPossibleHnCombinations( const std::vector<std::vector<unsigned int> >& potentialHnMasks,
      const std::vector<shFn::Octant> & elements,
      std::vector<std::vector<shFn::Octant> >& elemListWithHn ) {

    if(!elements.empty()) {
      shFn::Octant currOct = elements[0];
      shFn::Box parent = currOct.Parent();
      unsigned int currCnum = currOct.ChildNumber();

      for(unsigned int myHnType = 0; myHnType < 18; myHnType++) {  

        std::vector<shFn::Node> nodesVisited;
        std::vector<shFn::Octant> currRecord;
        //[0,0,0] always exists and it is always a node. 
        nodesVisited.push_back(shFn::Node(false));

        unsigned int hnMask = potentialHnMasks[currCnum][myHnType];     

        std::vector<shFn::Node> tmpNodes;		

        bool isConfigValid = IsConfigValid(hnMask, currOct, parent,
            nodesVisited, tmpNodes);

        if(isConfigValid) {	     
          InsertNodesAndOctant(myHnType, tmpNodes, nodesVisited, currOct, currRecord);
          tmpNodes.clear();		
          if(elements.size() > 1) {
            //Recurse
            AddPossibleHnCombinationsPrivate(potentialHnMasks,
                1,elements,nodesVisited,currRecord,elemListWithHn);
          }else {
            elemListWithHn.push_back(currRecord);
          }
        } 

        nodesVisited.clear();
        currRecord.clear();

      }//end loop over all hnTypes
    }

  }//end fn.

  //Assumes nodesVisited to be sorted and unique. This state is preserved.
  void AddPossibleHnCombinationsPrivate( const std::vector<std::vector<unsigned int> >& potentialHnMasks,
      unsigned int elemId, const std::vector<shFn::Octant> & elements,
      const std::vector<shFn::Node>& nodesVisitedIn,
      const std::vector<shFn::Octant> & currRecordIn,
      std::vector<std::vector<shFn::Octant> >& elemListWithHn ) {

    assert(elemId < elements.size());

    shFn::Octant currOct = elements[elemId];
    shFn::Box parent = currOct.Parent();
    unsigned int currCnum = currOct.ChildNumber();

    for(unsigned int myHnType = 0; myHnType < 18; myHnType++) {  

      unsigned int hnMask = potentialHnMasks[currCnum][myHnType];     

      std::vector<shFn::Node> tmpNodes;		

      bool isConfigValid = IsConfigValid(hnMask, currOct, parent, nodesVisitedIn, tmpNodes);

      if(isConfigValid) {	     
        std::vector<shFn::Node> nodesVisitedOut = nodesVisitedIn;
        std::vector<shFn::Octant> currRecordOut = currRecordIn;
        InsertNodesAndOctant(myHnType, tmpNodes, nodesVisitedOut, currOct, currRecordOut);
        tmpNodes.clear();		
        if(elements.size() > (elemId + 1) ) {
          //Recurse
          AddPossibleHnCombinationsPrivate(potentialHnMasks, (elemId+1), elements,
              nodesVisitedOut, currRecordOut, elemListWithHn);
        }else {
          //No more recursion.
          assert(currRecordOut.size() == elements.size());
          elemListWithHn.push_back(currRecordOut);
        }//end if-last element
      }        

    }//end loop over all hnTypes

  }//end fn.

  //Assumes nodesVisited to be sorted and unique. 
  bool IsConfigValid(unsigned int hnMask, const shFn::Octant& elem,
      const shFn::Box& parent, const std::vector<shFn::Node>& nodesVisited,
      std::vector<shFn::Node>& tmpNodes) {

    bool isConfigValid = true; 
    assert(tmpNodes.empty());

    for(unsigned int vtx = 0; vtx < 8; vtx++) {		
      if(hnMask & (1u<<vtx)) {
        shFn::Node vtxElem(elem.Vertex(vtx),true);
        shFn::Node vtxParent(parent.Vertex(vtx),false);
        unsigned int idx;
        bool foundElem = BinarySearch<Node>(&(*nodesVisited.begin()), nodesVisited.size(),  vtxElem, &idx);
        bool foundParent = BinarySearch<Node>(&(*nodesVisited.begin()), nodesVisited.size(),  vtxParent, &idx);
        if(foundElem) {
          if(!(vtxElem.Equals(nodesVisited[idx]))) {
            isConfigValid = false;
            break;		
          }
        } else {
          tmpNodes.push_back(vtxElem);
        }    
        if(foundParent) {
          if(!(vtxParent.Equals(nodesVisited[idx]))) {
            isConfigValid = false;
            break;		
          }
        } else {
          tmpNodes.push_back(vtxParent);
        }
      } else {
        shFn::Node vtxElem(elem.Vertex(vtx),false);
        unsigned int idx;
        bool foundElem = BinarySearch<Node>(&(*nodesVisited.begin()), nodesVisited.size(),  vtxElem, &idx);			
        if(foundElem) {
          if(!(vtxElem.Equals(nodesVisited[idx]))) {
            isConfigValid = false;
            break;		
          }
        } else {
          tmpNodes.push_back(vtxElem);
        }     	  
      }//check if hanging
    }//end for vtx    

    return isConfigValid;
  }//end fn.

  //Assumes nodesVisited to be sorted and unique. 
  void InsertNodesAndOctant(unsigned int myHnType, const std::vector<shFn::Node>& tmpNodes,
      std::vector<shFn::Node>& nodesVisited, shFn::Octant& currOct,
      std::vector<shFn::Octant>& currRecord ) {

    for(unsigned int tmpCtr = 0; tmpCtr < tmpNodes.size(); tmpCtr++) {
      unsigned int idx;
      bool found = maxLowerBound<Node>(nodesVisited, tmpNodes[tmpCtr], idx, NULL, NULL); 
      if(!found) { 			
        nodesVisited.insert(nodesVisited.begin(), tmpNodes[tmpCtr]);
      } else {
        assert(nodesVisited[idx] != tmpNodes[tmpCtr]); 
        nodesVisited.insert(nodesVisited.begin() + idx + 1, tmpNodes[tmpCtr]);
      }	         
    }

    currOct.SetHangingType(myHnType + 1); // 1-based numbering. Since 0 is used for the default case.
    currRecord.push_back(currOct);	
  }//end fn.

  unsigned int AddAuxilaryElems(unsigned int numOctFilesWritten) {
    //loop through the records and add auxilary elements depending
    //on the hanging types of the original elemnts.
    //If your i-th node is hanging the auxilary element will be the
    //i-th child of your parent.
    //These auxilary elements will be initialized with a 0 hanging flag.

    bool foundSecondaryConflict = false;

    unsigned int auxFileCtr = 0;
    std::vector<std::vector<shFn::Octant> > elemListWithAux;

    std::vector<std::vector<unsigned int> > potentialHnMasks(8);
    for(unsigned int i = 0; i< 8; i++) {
      GetCnumBasedHangingMasks( i, potentialHnMasks[i]);  
    }

    for(unsigned int octFileCtr = 0; octFileCtr < numOctFilesWritten; octFileCtr++) {
      char fname[256];
      std::vector<std::vector<shFn::Octant> > elemListWithHn;
      sprintf(fname,"octsWithHn_%u.txt",octFileCtr);
      readOctantsFromFile(fname, elemListWithHn);   

      for(unsigned int i = 0; i < elemListWithHn.size(); i++) {
        std::vector<shFn::Octant> extraOcts;

        for(unsigned int j = 0; j < elemListWithHn[i].size(); j++) {
          shFn::Octant currOct = elemListWithHn[i][j];
          assert(currOct.HangingType() > 0);
          unsigned int hnMask = potentialHnMasks[currOct.ChildNumber()][(currOct.HangingType())-1];
          for(unsigned int vtx = 0; vtx < 8; vtx++) {
            if(hnMask & (1u<<vtx)) {
              shFn::Octant auxElem = currOct.Brother(vtx);
              auxElem.SetHangingType(0);
              extraOcts.push_back(auxElem);
            }
          }
        }//end for j

        bool canAdd = true;
        std::vector<shFn::Octant> newRecord = elemListWithHn[i];

        //compare extraOcts with newRecord and then insert into newRecord.
        for(unsigned int j = 0; j < extraOcts.size(); j++) {
          bool isNew = true; 
          for(unsigned int k = 0; k < newRecord.size(); k++) {
            if(extraOcts[j].EqualsIgnoreHanging(newRecord[k])) {
              isNew = false;
              break;
            }
          }//end for k
          if(isNew) {
            //Check direct overlap
            for(unsigned int k = 0; k < newRecord.size(); k++) {
              bool weIntersect = extraOcts[j].Intersects(newRecord[k]);
              if(weIntersect) {
                canAdd = false;
                break;
              }
            }//end for k
            //check indirect overlap
            if(canAdd) {
              for(unsigned int k = 0; k < newRecord.size(); k++) {
                for(unsigned int cNumOfBroOfK = 0; cNumOfBroOfK < 8; cNumOfBroOfK++) {
                  if(cNumOfBroOfK == newRecord[k].ChildNumber()) {
                    continue;
                  }
                  shFn::Octant broOfK = newRecord[k].Brother(cNumOfBroOfK);
                  if(!broOfK.EqualsIgnoreHanging(extraOcts[j])) {
                    bool weIntersect = extraOcts[j].Intersects(broOfK);
                    if(weIntersect) {
                      canAdd = false;
                      break;
                    }
                  }
                }
                if(canAdd) {
                  for(unsigned int cNumOfBroOfJ = 0; cNumOfBroOfJ < 8; cNumOfBroOfJ++) {
                    if(cNumOfBroOfJ == extraOcts[j].ChildNumber()) {
                      continue;
                    }
                    shFn::Octant broOfJ = extraOcts[j].Brother(cNumOfBroOfJ);
                    if(!broOfJ.EqualsIgnoreHanging(newRecord[k])) {
                      bool weIntersect = newRecord[k].Intersects(broOfJ);
                      if(weIntersect) {
                        canAdd = false;
                        break;
                      }
                    }
                  }
                }
              }//end for k
              if(canAdd) {
                newRecord.push_back(extraOcts[j]);
              }
            }
          }//end if new
        }//end for j

        if(canAdd) {
          elemListWithAux.push_back(newRecord);
          if(elemListWithAux.size() > 5000) {
            sprintf(fname,"aux_%u.txt",auxFileCtr);
            writeOctantsToFile(fname, elemListWithAux);   
            elemListWithAux.clear();						
            auxFileCtr++;
          }
        }else {
          //I am not sure if this is possible or not.  
          foundSecondaryConflict = true;
        }

      }//end for i

      elemListWithHn.clear();
    }//end for octFileCtr

    if(elemListWithAux.size() > 0) {
      char fname[256];
      sprintf(fname,"aux_%u.txt",auxFileCtr);
      writeOctantsToFile(fname, elemListWithAux);   
      elemListWithAux.clear();						
      auxFileCtr++;
    }

    if(foundSecondaryConflict) {
      std::cout<<"Found a Secondary Conflict."<<std::endl;
    }

    return auxFileCtr;

  }//end fn.

  unsigned int FinalElemList(unsigned int numAuxFilesWritten) {
    //Assign all possible hanging flags to the auxilary elements
    unsigned int finalElemCtr = 0;
    std::vector<std::vector<shFn::Octant> > finalElemList;

    std::vector<std::vector<unsigned int> > potentialHnMasks(8);
    for(unsigned int i = 0; i< 8; i++) {
      GetCnumBasedHangingMasks( i, potentialHnMasks[i]);  
    }

    for(unsigned int auxFileCtr = 0;  auxFileCtr < numAuxFilesWritten; auxFileCtr++) {
      char fname[256];
      std::vector<std::vector<shFn::Octant> > elemListWithAux;
      sprintf(fname,"aux_%u.txt",auxFileCtr);
      readOctantsFromFile(fname, elemListWithAux);   

      for(unsigned int i = 0; i < elemListWithAux.size(); i++) {
        //Recursion...
        AddPossibleHnCombinationsForAux(potentialHnMasks, elemListWithAux[i], finalElemList);
        if(finalElemList.size() > 1000) {
          sprintf(fname,"finalElem_%u.txt",finalElemCtr);
          writeOctantsToFile(fname, finalElemList);   
          finalElemList.clear();						
          finalElemCtr++;
        }
      }//end for i

      elemListWithAux.clear();						
    }//end for auxFileCtr

    if(finalElemList.size() > 0) {
      char fname[256];
      sprintf(fname,"finalElem_%u.txt",finalElemCtr);
      writeOctantsToFile(fname, finalElemList);   
      finalElemList.clear();						
      finalElemCtr++;
    }

    return finalElemCtr;
  }//end fn.

  //Assumes nodesVisited to be sorted and unique. 
  void AddPossibleHnCombinationsForAux( const std::vector<std::vector<unsigned int> >& potentialHnMasks,
      const std::vector<shFn::Octant> & elements,
      std::vector<std::vector<shFn::Octant> >& elemListWithHn ) {

    unsigned int elemId = 0;

    std::vector<shFn::Node> nodesVisited;
    std::vector<shFn::Octant> currRecord;

    while( (elemId < elements.size()) && (elements[elemId].HangingType() > 0) ) {
      shFn::Octant currOct = elements[elemId];
      shFn::Box parent = currOct.Parent();
      unsigned int currCnum = currOct.ChildNumber();
      unsigned int myHnType = currOct.HangingType();
      unsigned int hnMask = potentialHnMasks[currCnum][myHnType-1];     
      currRecord.push_back(currOct);
      for(unsigned int vtx = 0; vtx < 8; vtx++) {
        if( hnMask & (1u<<vtx) ) {
          nodesVisited.push_back(shFn::Node(currOct.Vertex(vtx),true));
          nodesVisited.push_back(shFn::Node(parent.Vertex(vtx),false));
        }else {
          nodesVisited.push_back(shFn::Node(currOct.Vertex(vtx),false));
        }
      }//end for vtx
      elemId++;
    }

    //Test all the elements in the tailing end have 0 hangingTypes. These must
    //be the aux elements.
    for(unsigned int i = elemId; i < elements.size(); i++) {
      assert(elements[i].HangingType() == 0);
    }

    if(elemId < elements.size()) {
      //The comparator here ignores the hanging type. However, we know
      //that there are no conflicts. Hence, whenever two points are equal, then
      //the hanging types must also be equal.
      shFn::makeVectorUnique<shFn::Node>(nodesVisited,false);
      //Recurse
      AddPossibleHnCombinationsPrivate(potentialHnMasks, elemId, elements,
          nodesVisited, currRecord, elemListWithHn);
    }else {
      elemListWithHn.push_back(currRecord);
    }

    nodesVisited.clear();
    currRecord.clear();
  }//end fn.

  void FinalVerticesList(unsigned int numFinalElemFilesWritten, 
      unsigned int & minCnt, unsigned int & maxCnt) {
    //Loop through each set and store the non-hanging vertices. Make this set of vertices unique.
    //Count the number of unique vertices in each set. Find the max and min among all the sets.
    bool isFirstSet = true;

    std::vector<std::vector<unsigned int> > potentialHnMasks(8);
    for(unsigned int i = 0; i< 8; i++) {
      GetCnumBasedHangingMasks( i, potentialHnMasks[i]);  
    }

    for(unsigned int finalElemCtr = 0;
        finalElemCtr < numFinalElemFilesWritten; finalElemCtr++) {
      char fname[256];
      std::vector<std::vector<shFn::Octant> > finalElemList;
      sprintf(fname,"finalElem_%u.txt",finalElemCtr);
      readOctantsFromFile(fname, finalElemList);   

      for(unsigned int i = 0; i < finalElemList.size(); i++) {

        std::vector<shFn::Point> vertices;

        for(unsigned int j = 0; j < finalElemList[i].size(); j++) {
          shFn::Octant currOct = finalElemList[i][j];
          shFn::Box parent = currOct.Parent();
          unsigned int currCnum = currOct.ChildNumber();
          unsigned int myHnType = currOct.HangingType();
          assert(myHnType > 0);
          unsigned int hnMask = potentialHnMasks[currCnum][myHnType-1];     
          for(unsigned int vtx = 0; vtx < 8; vtx++) {
            if(!(hnMask & (1u<<vtx))) {
              vertices.push_back(currOct.Vertex(vtx));
            }else {
              vertices.push_back(parent.Vertex(vtx));
            }
          }//end for vtx
        }//end for j

        shFn::makeVectorUnique<shFn::Point>(vertices,false);
        unsigned int numVtxInCurrSet = vertices.size();

        if(isFirstSet) {
          minCnt = numVtxInCurrSet;
          maxCnt = numVtxInCurrSet;
          isFirstSet = false;
        }else {
          if(minCnt > numVtxInCurrSet) {
            minCnt = numVtxInCurrSet;
          }
          if(maxCnt < numVtxInCurrSet) {
            maxCnt = numVtxInCurrSet;
          }
        }

        vertices.clear();
      }//end for i

      finalElemList.clear();						
    }//end for finalElemCtr

  }//end fn.

  void RunAllInCore(const std::vector<std::vector<
      std::vector<shFn::Octant> > > & nhCnumAndSz,
      unsigned int & minCnt, unsigned int & maxCnt) {

    shFn::Box coordsTemplate(shFn::Point(),1.0);
    std::vector<std::vector<shFn::Octant> > elemList;

    //First list all possible sets of elements to vist.
    //This is basically selecting one from each elemNum type
    //Each element has the following info: anchor, cNum, sz
    for(unsigned int myCnum = 0; myCnum < 8; myCnum++) {
      std::cout<<"Generating elemList for cNum: "<<myCnum<<std::endl;
      unsigned int lenNh[7];
      for(unsigned int i = 0; i < 7; i++) {
        lenNh[i] = nhCnumAndSz[myCnum][i].size();
        std::cout<<" Len("<<i<<"): "<<lenNh[i]<<std::endl;
      }
      unsigned int elemNumIndices[7];
      for(unsigned int e1Idx = 0; e1Idx < lenNh[0]; e1Idx++) {	
        elemNumIndices[0] = e1Idx;
        for(unsigned int e2Idx = 0; e2Idx < lenNh[1]; e2Idx++) {
          elemNumIndices[1] = e2Idx;
          for(unsigned int e3Idx = 0; e3Idx < lenNh[2]; e3Idx++) {
            elemNumIndices[2] = e3Idx;
            for(unsigned int e4Idx = 0; e4Idx < lenNh[3]; e4Idx++) {
              elemNumIndices[3] = e4Idx;
              for(unsigned int e5Idx = 0; e5Idx < lenNh[4]; e5Idx++) {
                elemNumIndices[4] = e5Idx;
                for(unsigned int e6Idx = 0; e6Idx < lenNh[5]; e6Idx++) {
                  elemNumIndices[5] = e6Idx;
                  for(unsigned int e7Idx = 0; e7Idx < lenNh[6]; e7Idx++) {
                    elemNumIndices[6] = e7Idx;
                    //New Set
                    std::vector<shFn::Octant> record;
                    record.push_back(shFn::Octant(coordsTemplate,myCnum,0));
                    for(unsigned int elemNum = 1; elemNum <8; elemNum++) {
                      shFn::Octant currOct = 
                        nhCnumAndSz[myCnum][elemNum-1][elemNumIndices[elemNum-1]];
                      assert(currOct.Vertex(elemNum) == Point());
                      assert(currOct.HangingType() == 0);
                      record.push_back(currOct);
                    }
                    elemList.push_back(record);
                    if(elemList.size() > 5000) {
                      shFn::EliminateOverlapsInCore(elemList, minCnt, maxCnt);
                      elemList.clear();						
                    }
                  }//e7
                }//e6
              }//e5	
            }//e4
          }//e3		
        }//e2
      }//e1       	
    }//myCnum

    if(elemList.size() > 0) {
      shFn::EliminateOverlapsInCore(elemList, minCnt, maxCnt);
      elemList.clear();						
    }
  }//end fn.

  void EliminateOverlapsInCore(const std::vector<
      std::vector<shFn::Octant> > &elemList,
      unsigned int & minCnt, unsigned int & maxCnt) {

    static unsigned int overlapInCoreCtr = 0;
    if( (overlapInCoreCtr % 100) == 0) {
      std::cout<<"Eliminating Overlaps."<<std::endl;
    }
    overlapInCoreCtr++;

    std::vector<std::vector<shFn::Octant> > validElemList;

    for(unsigned int i = 0; i < elemList.size(); i++) {
      bool foundDirectOverlap = FoundDirectOverlap(elemList[i]);
      if(!foundDirectOverlap) {
        bool foundInDirectOverlap = FoundInDirectOverlap(elemList[i]);
        if(!foundInDirectOverlap) {
          validElemList.push_back(elemList[i]);
          if(validElemList.size() > 5000) {
            shFn::IncludeHangingTypesInCore(validElemList, minCnt, maxCnt);
            validElemList.clear();						
          }
        }
      }
    }//end for i

    if(validElemList.size() > 0) {
      shFn::IncludeHangingTypesInCore(validElemList, minCnt, maxCnt);
      validElemList.clear();						
    }
  }//end fn

  void IncludeHangingTypesInCore(const std::vector<
      std::vector<shFn::Octant> > &validElemList,
      unsigned int & minCnt, unsigned int & maxCnt) {
    //Until now all the hanging flags are initialized to 0.
    //Now we create a list of all the possible elements along
    //with their respective hanging types.
    //To do this basically keep track of the vertices you have
    //visited and their hanging flag. [0,0,0] must never be hanging.

    static unsigned int hangingInCoreCtr = 0;
    if( (hangingInCoreCtr % 100) == 0) {
      std::cout<<"Including Hanging Types."<<std::endl;
    }
    hangingInCoreCtr++;

    std::vector<std::vector<shFn::Octant> > elemListWithHn;

    std::vector<std::vector<unsigned int> > potentialHnMasks(8);
    for(unsigned int i = 0; i< 8; i++) {
      GetCnumBasedHangingMasks( i, potentialHnMasks[i]);  
    }

    for(unsigned int i = 0; i < validElemList.size(); i++) {
      //Recursion...
      AddPossibleHnCombinations(potentialHnMasks, validElemList[i], elemListWithHn);
      if(elemListWithHn.size() > 5000) {
        shFn::AddAuxilaryElemsInCore(elemListWithHn, minCnt, maxCnt);
        elemListWithHn.clear();						
      }
    }//end for i

    if(elemListWithHn.size() > 0) {
      shFn::AddAuxilaryElemsInCore(elemListWithHn, minCnt, maxCnt);
      elemListWithHn.clear();						
    }
  }//end fn.

  void AddAuxilaryElemsInCore(const std::vector<
      std::vector<shFn::Octant> > &elemListWithHn,
      unsigned int & minCnt, unsigned int & maxCnt) {
    //loop through the records and add auxilary elements depending
    //on the hanging types of the original elemnts.
    //If your i-th node is hanging the auxilary element will be the
    //i-th child of your parent.
    //These auxilary elements will be initialized with a 0 hanging flag.

    static unsigned int auxInCoreCtr = 0;
    if( (auxInCoreCtr % 100) == 0) {
      std::cout<<"Adding Auxilary Elements."<<std::endl;
    }
    auxInCoreCtr++;

    std::vector<std::vector<shFn::Octant> > elemListWithAux;

    std::vector<std::vector<unsigned int> > potentialHnMasks(8);
    for(unsigned int i = 0; i< 8; i++) {
      GetCnumBasedHangingMasks( i, potentialHnMasks[i]);  
    }

    for(unsigned int i = 0; i < elemListWithHn.size(); i++) {
      std::vector<shFn::Octant> extraOcts;

      for(unsigned int j = 0; j < elemListWithHn[i].size(); j++) {
        shFn::Octant currOct = elemListWithHn[i][j];
        assert(currOct.HangingType() > 0);
        unsigned int hnMask = 
          potentialHnMasks[currOct.ChildNumber()][(currOct.HangingType())-1];
        for(unsigned int vtx = 0; vtx < 8; vtx++) {
          if(hnMask & (1u<<vtx)) {
            shFn::Octant auxElem = currOct.Brother(vtx);
            auxElem.SetHangingType(0);
            extraOcts.push_back(auxElem);
          }
        }
      }//end for j

      bool canAdd = true;
      std::vector<shFn::Octant> newRecord = elemListWithHn[i];

      //compare extraOcts with newRecord and then insert into newRecord.
      for(unsigned int j = 0; j < extraOcts.size(); j++) {
        bool isNew = true; 
        for(unsigned int k = 0; k < newRecord.size(); k++) {
          if(extraOcts[j].EqualsIgnoreHanging(newRecord[k])) {
            isNew = false;
            break;
          }
        }//end for k
        if(isNew) {
          //Check direct overlap
          for(unsigned int k = 0; k < newRecord.size(); k++) {
            bool weIntersect = extraOcts[j].Intersects(newRecord[k]);
            if(weIntersect) {
              canAdd = false;
              break;
            }
          }//end for k
          //check indirect overlap
          if(canAdd) {
            for(unsigned int k = 0; k < newRecord.size(); k++) {
              for(unsigned int cNumOfBroOfK = 0; cNumOfBroOfK < 8; cNumOfBroOfK++) {
                if(cNumOfBroOfK == newRecord[k].ChildNumber()) {
                  continue;
                }
                shFn::Octant broOfK = newRecord[k].Brother(cNumOfBroOfK);
                if(!broOfK.EqualsIgnoreHanging(extraOcts[j])) {
                  bool weIntersect = extraOcts[j].Intersects(broOfK);
                  if(weIntersect) {
                    canAdd = false;
                    break;
                  }
                }
              }
              if(canAdd) {
                for(unsigned int cNumOfBroOfJ = 0; cNumOfBroOfJ < 8; cNumOfBroOfJ++) {
                  if(cNumOfBroOfJ == extraOcts[j].ChildNumber()) {
                    continue;
                  }
                  shFn::Octant broOfJ = extraOcts[j].Brother(cNumOfBroOfJ);
                  if(!broOfJ.EqualsIgnoreHanging(newRecord[k])) {
                    bool weIntersect = newRecord[k].Intersects(broOfJ);
                    if(weIntersect) {
                      canAdd = false;
                      break;
                    }
                  }
                }
              }
            }//end for k
            if(canAdd) {
              newRecord.push_back(extraOcts[j]);
            }
          }
        }//end if new
      }//end for j

      if(canAdd) {
        elemListWithAux.push_back(newRecord);
        if(elemListWithAux.size() > 5000) {
          shFn::FinalElemListInCore(elemListWithAux, minCnt, maxCnt);
          elemListWithAux.clear();						
        }
      }

    }//end for i

    if(elemListWithAux.size() > 0) {
      shFn::FinalElemListInCore(elemListWithAux, minCnt, maxCnt);
      elemListWithAux.clear();						
    }
  }//end fn.

  void FinalElemListInCore(const std::vector<
      std::vector<shFn::Octant> > &elemListWithAux,
      unsigned int & minCnt, unsigned int & maxCnt) {
    //Assign all possible hanging flags to the auxilary elements

    static unsigned int finalElemInCoreCtr = 0;
    if( (finalElemInCoreCtr % 100) == 0) {
      std::cout<<"Building Final Element Lists."<<std::endl;
    }
    finalElemInCoreCtr++;

    std::vector<std::vector<shFn::Octant> > finalElemList;

    std::vector<std::vector<unsigned int> > potentialHnMasks(8);
    for(unsigned int i = 0; i< 8; i++) {
      GetCnumBasedHangingMasks( i, potentialHnMasks[i]);  
    }

    for(unsigned int i = 0; i < elemListWithAux.size(); i++) {
      //Recursion...
      AddPossibleHnCombinationsForAux(potentialHnMasks, elemListWithAux[i], finalElemList);
      if(finalElemList.size() > 5000) {
        shFn::FinalVerticesListInCore(finalElemList, minCnt, maxCnt);
        finalElemList.clear();						
      }
    }//end for i

    if(finalElemList.size() > 0) {
      shFn::FinalVerticesListInCore(finalElemList, minCnt, maxCnt);
      finalElemList.clear();						
    }
  }//end fn.

  void FinalVerticesListInCore(const std::vector<
      std::vector<shFn::Octant> > & finalElemList,
      unsigned int & minCnt, unsigned int & maxCnt) {
    //Loop through each set and store the non-hanging vertices.
    //Make this set of vertices unique.
    //Count the number of unique vertices in each set.
    //Find the max and min among all the sets.

    static unsigned int finalVtxInCoreCtr = 0;
    if( (finalVtxInCoreCtr % 100) == 0) {
      std::cout<<"Building Final Vertices Lists."<<std::endl;
    }
    finalVtxInCoreCtr++;

    static bool isFirstSet = true;

    std::vector<std::vector<unsigned int> > potentialHnMasks(8);
    for(unsigned int i = 0; i< 8; i++) {
      GetCnumBasedHangingMasks( i, potentialHnMasks[i]);  
    }

    for(unsigned int i = 0; i < finalElemList.size(); i++) {
      std::vector<shFn::Point> vertices;

      for(unsigned int j = 0; j < finalElemList[i].size(); j++) {
        shFn::Octant currOct = finalElemList[i][j];
        shFn::Box parent = currOct.Parent();
        unsigned int currCnum = currOct.ChildNumber();
        unsigned int myHnType = currOct.HangingType();
        assert(myHnType > 0);
        unsigned int hnMask = potentialHnMasks[currCnum][myHnType-1];     
        for(unsigned int vtx = 0; vtx < 8; vtx++) {
          if(!(hnMask & (1u<<vtx))) {
            vertices.push_back(currOct.Vertex(vtx));
          }else {
            vertices.push_back(parent.Vertex(vtx));
          }
        }//end for vtx
      }//end for j

      shFn::makeVectorUnique<shFn::Point>(vertices,false);
      unsigned int numVtxInCurrSet = vertices.size();

      if(isFirstSet) {
        minCnt = numVtxInCurrSet;
        maxCnt = numVtxInCurrSet;
        isFirstSet = false;
        std::cout<<"Initial minCnt = maxCnt = "<<numVtxInCurrSet<<std::endl;
      }else {
        if(minCnt > numVtxInCurrSet) {
          minCnt = numVtxInCurrSet;
          std::cout<<"minCnt updated: "<<minCnt<<std::endl;
        }
        if(maxCnt < numVtxInCurrSet) {
          maxCnt = numVtxInCurrSet;
          std::cout<<"maxCnt updated: "<<maxCnt<<std::endl;
        }
      }

      vertices.clear();
    }//end for i

  }//end fn.


}//end namespace

int main(int argc, char**argv) {

  MPI_Init(&argc, &argv);

  std::cout<<"Starting..."<<std::endl;

  bool outOfCore = (bool)atoi(argv[1]);

  std::vector<std::vector<std::vector<shFn::Octant> > > nhCnumAndSz;
  shFn::GetNeighboringCnumAndSizes(nhCnumAndSz); 
  std::cout<<"Got Neighboring Cnums and Sizes"<<std::endl;

  unsigned int minNumVtx = 0;
  unsigned int maxNumVtx = 0;
  if(outOfCore) {
    std::cout<<" Running out of core..."<<std::endl;

    //First list all possible sets of elements to vist.
    //This is basically selecting one from each elemNum type
    //Each element has the following info: anchor, cNum, sz
    unsigned int numElemFilesWritten = shFn::ListAllElemsSharing0(nhCnumAndSz);
    nhCnumAndSz.clear();
    std::cout<<"Got List of all elems sharing 0."<<std::endl;

    //Now we have a list of all elements which share the node [0,0,0]

    //Loop through the above list and eliminate cases which are not possible
    //Eliminate direct overlaps (2 octants in the list intersect),
    //Eliminate indirect overlaps(2 octants do not directly intersect,
    //but they do not have comptatible sizes and cNums)
    //Unlike direct overlaps, indirect overlaps are not symmetric. i.e.
    // A.Intersects(B) and B.Intersects(A) are always the same.
    // However, for indirect overlaps we must check if any of A's brothers
    // intersect with B and if any of B's brothers intersect with A. Both the
    // tests are necessary. Just like in getNhCnumAndSz, here too valid brothers
    // must not be treated as intersects. 

    unsigned int numValidElemFilesWritten = shFn::EliminateOverlaps(numElemFilesWritten);
    std::cout<<"Eliminated overlaps."<<std::endl;

    //Now we have a list of valid elements which share the node [0,0,0]

    //Until now all the hanging flags are initialized to 0.
    //Now we create a list of all the possible elements along
    //with their respective hanging types.
    //To do this basically keep track of the vertices you have
    //visited and their hanging flag. [0,0,0] must never be hanging.

    unsigned int numOctFilesWritten = shFn::IncludeHangingTypes(numValidElemFilesWritten);
    std::cout<<"Included Hanging Types."<<std::endl;

    //Now we have a list of valid elements along with their hanging
    //types which share the node [0,0,0]

    //loop through the records and add auxilary elements depending
    //on the hanging types of the original elemnts.
    //If your i-th node is hanging the auxilary element will be the
    //i-th child of your parent.
    //These auxilary elements will be initialized with a 0 hanging flag.

    unsigned int numAuxFilesWritten = shFn::AddAuxilaryElems(numOctFilesWritten);
    std::cout<<"Included Auxilary Elements."<<std::endl;

    //Now we have a list of all elements within the support of the
    //shape function rooted at [0,0,0] 

    //Assign all possible hanging flags to the auxilary elements
    unsigned int numFinalElemFilesWritten = shFn::FinalElemList(numAuxFilesWritten);
    std::cout<<"Final Element list ready."<<std::endl;

    //Loop through each set and store the non-hanging vertices. Make this set of vertices unique.
    //Count the number of unique vertices in each set. Find the max among all the sets.
    shFn::FinalVerticesList(numFinalElemFilesWritten, minNumVtx, maxNumVtx);
  } else {
    std::cout<<" Running in core..."<<std::endl;

    shFn::RunAllInCore(nhCnumAndSz, minNumVtx, maxNumVtx);

    nhCnumAndSz.clear();
  }

  std::cout<<"minVtxCount: "<<minNumVtx<<", maxVtxCount: "<<maxNumVtx<<std::endl;

  std::cout<<"Finished Successfully."<<std::endl;

  MPI_Finalize();

}


