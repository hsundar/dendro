clc;clear;
%X right, Y Back Z top
%1-based indexing

xlb = [-1 0 -1 0 -1 0 -1 0];
xub = [0 1 0 1 0 1 0 1];
ylb = [-1 -1 0 0 -1 -1 0 0];
yub = [0 0 1 1 0 0 1 1];
zlb = [-1 -1 -1 -1 0 0 0 0];
zub = [0 0 0 0 1 1 1 1];

%Type-1 Coarse and Fine are not the same.
minNonZero = 0;
maxNonZero = 0;
for cNumCoarse=1:8
    cHnMasks = cNumBasedHangingMasks(cNumCoarse);
    ParentCoords = parentOfCoarse3d(cNumCoarse);
    for coarseType=1:18
        phiCoarse = coarseShapeFunction3DCoeffs(cHnMasks(coarseType,:),ParentCoords);
        for cNumFine=1:8
            Rmat = zeros(8,8);
            xl = xlb(cNumFine);
            yl = ylb(cNumFine);
            zl = zlb(cNumFine);
            xu = xub(cNumFine);
            yu = yub(cNumFine);
            zu = zub(cNumFine);
            for i=1:8
                    cSh = phiCoarse(i,:);
                    Rmat(i,1) = evalPhi(cSh,xl,yl,zl);
                    Rmat(i,2) = evalPhi(cSh,xu,yl,zl);
                    Rmat(i,3) = evalPhi(cSh,xl,yu,zl);
                    Rmat(i,4) = evalPhi(cSh,xu,yu,zl);
                    Rmat(i,5) = evalPhi(cSh,xl,yl,zu);
                    Rmat(i,6) = evalPhi(cSh,xu,yl,zu);
                    Rmat(i,7) = evalPhi(cSh,xl,yu,zu);
                    Rmat(i,8) = evalPhi(cSh,xu,yu,zu);
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
            fname = ['RmatType1_',int2str(coarseType),'_',int2str(cNumCoarse),...
                  '_',int2str(cNumFine),'.mat'];
            save(fname,'Rmat','thisMin','thisMax');
        end
    end
end

display(minNonZero)
display(maxNonZero)


