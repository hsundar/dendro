clc;clear;
%X right, Y Back Z top
%1-based indexing

%The global coordinate system and the local coordinate system are aligned. So
%the jacobian is only a diagonal matrix with entries h/2.
%Type-2 Coarse and Fine are the same. The integration is on the domain
%[-1,1] in the local coordinate system (master element).

%Gauss-point (3-pt rule)
%w1 =8/9, x1 =0, w2=w3 =5/9, x2 = +sqrt(3/5), x3 = -sqrt(3/5)
gpt = [0, sqrt(3/5), -sqrt(3/5)];

for cNum=1:8
    HnMasks = cNumBasedHangingMasks(cNum);
    ParentCoords = parentOfCoarse3d(cNum);
    for eType=1:18
        phiCoarse = coarseShapeFunction3DCoeffs(HnMasks(eType,:),ParentCoords);
        shFnMat = zeros(8,3,3,3);
        for j = 1:8
            cSh_j = phiCoarse(j,:);
            for m=1:3
                for n=1:3
                    for p=1:3
                        shFnMat(j,m,n,p) = evalPhi(cSh_j,gpt(m),gpt(n),gpt(p));
                    end
                end
            end
        end
        fname = ['shFnMat_',int2str(cNum),'_',int2str(eType),'.mat'];
        save(fname,'shFnMat');
    end
end

