clc;clear;
%X right, Y Back Z top
%1-based indexing

%Type-2 Coarse and Fine are the same. The integration is on the domain
%[-1,1] in the local coordinate system (master element).

minNonZero = 0;
maxNonZero = 0;
for cNum=1:8
    HnMasks = cNumBasedHangingMasks(cNum);
    ParentCoords = parentOfCoarse3d(cNum);
    for coarseType=1:18
        phiCoarse = coarseShapeFunction3DCoeffs(HnMasks(coarseType,:),ParentCoords);
        Rmat = zeros(8,8);
        for i=1:8
            cSh = phiCoarse(i,:);
            Rmat(i,1) = evalPhi(cSh,-1,-1,-1);
            Rmat(i,2) = evalPhi(cSh, 1,-1,-1);
            Rmat(i,3) = evalPhi(cSh,-1, 1,-1);
            Rmat(i,4) = evalPhi(cSh, 1, 1,-1);
            Rmat(i,5) = evalPhi(cSh,-1,-1, 1);
            Rmat(i,6) = evalPhi(cSh, 1,-1, 1);
            Rmat(i,7) = evalPhi(cSh,-1, 1, 1);
            Rmat(i,8) = evalPhi(cSh, 1, 1, 1);
        end
        thisNz = nonzeros(Rmat);
        thisMin = min(thisNz);
        thisMax = max(thisNz);
        if(minNonZero) 
           if(thisMin < minNonZero)
              minNonZero = thisMin;
           end
        else
           minNonZero = thisMin;
        end
        if(thisMax > maxNonZero)
           maxNonZero = thisMax;
        end
        fname = ['RmatType2_',int2str(cNum),'_',int2str(coarseType),'.mat'];
        save(fname,'Rmat','thisMin','thisMax');
    end
end

display(minNonZero)
display(maxNonZero)

