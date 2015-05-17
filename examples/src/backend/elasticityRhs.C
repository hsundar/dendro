
/**
  @file elasticityRhs.C
  @author Rahul Sampath, rahul.sampath
  */

#include "petsc.h"
#include "petscvec.h"
#include "omg.h"
#include "oda.h"
#include "elasticityJac.h"
#include "vecMass.h"

PetscErrorCode ComputeElasticityRHS(ot::DAMG damg, Vec rhs) {
  PetscFunctionBegin;	 	 
  ot::DA* da = damg->da;
  Vec tmp;
  VecDuplicate(rhs, &tmp);
  PetscScalar *inarray;
  VecZeroEntries(tmp);
  da->vecGetBuffer(tmp, inarray, false, false, false, 3);

  ElasticityData* data = (static_cast<ElasticityData*>(damg->user));
  unsigned char* bdyArr = data->bdyArr;
  unsigned int maxD;
  unsigned int balOctmaxD;
  if(da->iAmActive()) {
    maxD = da->getMaxDepth();
    balOctmaxD = maxD - 1;
    for(da->init<ot::DA_FLAGS::ALL>(); da->curr() < da->end<ot::DA_FLAGS::ALL>(); da->next<ot::DA_FLAGS::ALL>())  
    {
      Point pt;
      pt = da->getCurrentOffset();
      unsigned levelhere = da->getLevel(da->curr()) - 1;
      double hxOct = (double)((double)(1u << (balOctmaxD - levelhere))/(double)(1u << balOctmaxD));
      double x = (double)(pt.xint())/((double)(1u << (maxD-1)));
      double y = (double)(pt.yint())/((double)(1u << (maxD-1)));
      double z = (double)(pt.zint())/((double)(1u << (maxD-1)));

      unsigned int indices[8];
      da->getNodeIndices(indices); 
      double coord[8][3] = {
        {0.0,0.0,0.0},
        {1.0,0.0,0.0},
        {0.0,1.0,0.0},
        {1.0,1.0,0.0},
        {0.0,0.0,1.0},
        {1.0,0.0,1.0},
        {0.0,1.0,1.0},
        {1.0,1.0,1.0}
      };
      unsigned char hn = da->getHangingNodeIndex(da->curr());

      for(int i = 0; i < 8; i++)
      {
        if (!(hn & (1 << i))){
          //Dirichlet nodes will be set to 0
          if(!(bdyArr[indices[i]])) {
            //Some dummy loading with Fx = Fy = Fz at all pts.
            double xhere, yhere, zhere;
            xhere = x + coord[i][0]*hxOct ;
            yhere = y + coord[i][1]*hxOct;
            zhere = z + coord[i][2]*hxOct; 
            double rhsSum = 0.0;
            for(int freqCnt = 1; freqCnt < 10; freqCnt++) {
              double facsol = freqCnt;
              rhsSum  += (1.0 + 3*facsol*facsol*M_PI*M_PI)*cos(facsol*M_PI*xhere)*
                cos(facsol*M_PI*yhere)*cos(facsol*M_PI*zhere);				               
            }
            inarray[3*indices[i]] = rhsSum;
            inarray[(3*indices[i])+1] = rhsSum;
            inarray[(3*indices[i])+2] = rhsSum;
          }
        }
      }
    }
  }

  da->vecRestoreBuffer(tmp,inarray,false,false,false,3); 

  Mat vecMassMat;
  CreateConstVecMass(damg, &vecMassMat);
  ComputeConstVecMass(damg, vecMassMat, vecMassMat);

  MatMult(vecMassMat,tmp,rhs);
  MatDestroy(&vecMassMat);
  VecDestroy(&tmp);

  VecScale(rhs,-1.0);

  PetscFunctionReturn(0);
}


