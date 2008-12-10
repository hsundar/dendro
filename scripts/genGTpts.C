
#include <iostream>
#include <cstdio>
#include <vector>
#include <cassert>

void writePtsToFile( const double* pts, const unsigned int ptsLen, char* filename) {
  FILE*  outfile = fopen(filename,"w");
  if(ptsLen >0) {
    unsigned int numPts = ptsLen/3;
    fwrite(&numPts, sizeof(unsigned int),1,outfile);
    fwrite(pts, sizeof(double),ptsLen,outfile);
  }
  fclose(outfile);
}//end function

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "<exe> outFile h1 h2" << std::endl;
    return -1;
  }

  double h1 = atof(argv[2]);
  double h2 = atof(argv[3]);

  assert(h1 <= 0.1);
  assert(h2 <= 0.05);

  std::vector<double> outPts;

  for(double ty = 0.3; ty <= 0.7; ty += h2) {
    //G top pts
    for(double tx = 0.05; tx <= 0.45; tx += h1) {
      for(double tz = 0.7; tz <= 0.8; tz += h2) {
        outPts.push_back(tx);
        outPts.push_back(ty);
        outPts.push_back(tz);
      }
    }

    //G column pts
    for(double tx = 0.05; tx <= 0.15; tx += h1) {
      for(double tz = 0.2; tz <= 0.7; tz += h2) {
        outPts.push_back(tx);
        outPts.push_back(ty);
        outPts.push_back(tz);
      }
    }

    //G bottom part 1 pts
    for(double tx = 0.15; tx <= 0.40; tx += h1) {
      for(double tz = 0.2; tz <= 0.3; tz += h2) {
        outPts.push_back(tx);
        outPts.push_back(ty);
        outPts.push_back(tz);
      }
    }

    //G bottom part 2 pts
    for(double tx = 0.35; tx <= 0.40; tx += h1) {
      for(double tz = 0.30; tz <= 0.50; tz += h2) {
        outPts.push_back(tx);
        outPts.push_back(ty);
        outPts.push_back(tz);
      }
    }

    //G bottom part 3 pts
    for(double tx = 0.40; tx <= 0.475; tx += h1) {
      for(double tz = 0.40; tz <= 0.50; tz += h2) {
        outPts.push_back(tx);
        outPts.push_back(ty);
        outPts.push_back(tz);
      }
    }

    //G bottom part 4 pts
    for(double tx = 0.425; tx <= 0.475; tx += h1) {
      for(double tz = 0.20; tz <= 0.40; tz += h2) {
        outPts.push_back(tx);
        outPts.push_back(ty);
        outPts.push_back(tz);
      }
    }

    //T top pts
    for(double tx = 0.55; tx <= 0.95; tx += h1) {
      for(double tz = 0.7; tz <= 0.8; tz += h2) {
        outPts.push_back(tx);
        outPts.push_back(ty);
        outPts.push_back(tz);
      }
    }

    //T leg pts
    for(double tx = 0.7; tx <= 0.8; tx += h1) {
      for (double tz = 0.2; tz <= 0.7; tz += h2) {
        outPts.push_back(tx);
        outPts.push_back(ty);
        outPts.push_back(tz);
      }
    }

  }

  std::cout<<"Writing "<<(outPts.size()/3)<<" points."<<std::endl;

  writePtsToFile(&(*outPts.begin()), outPts.size(), argv[1]);

  outPts.clear();

  return 0;
}



