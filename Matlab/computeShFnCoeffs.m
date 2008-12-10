clc;clear;
%X right, Y Back Z top
%1-based indexing

%The global coordinate system and the local coordinate system are aligned. So
%the jacobian is only a diagonal matrix with entries h/2.
%Type-2 Coarse and Fine are the same. The integration is on the domain
%[-1,1] in the local coordinate system (master element).

for cNum=1:8
    HnMasks = cNumBasedHangingMasks(cNum);
    ParentCoords = parentOfCoarse3d(cNum);
    for eType=1:18
        ShFnCoeffs = coarseShapeFunction3DCoeffs(HnMasks(eType,:),ParentCoords);
        fname = ['ShFnCoeffs_',int2str(cNum),'_',int2str(eType),'.mat'];
        save(fname,'ShFnCoeffs');
    end
end

