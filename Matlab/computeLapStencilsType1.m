clc;clear;
%X right, Y Back Z top
%1-based indexing

%The global coordinate system and the local coordinate system are aligned. So
%the jacobian is only a diagonal matrix with entries h/2.

xlb = [-1 0 -1 0 -1 0 -1 0];
xub = [0 1 0 1 0 1 0 1];
ylb = [-1 -1 0 0 -1 -1 0 0];
yub = [0 0 1 1 0 0 1 1];
zlb = [-1 -1 -1 -1 0 0 0 0];
zub = [0 0 0 0 1 1 1 1];

%Type-1 Coarse and Fine are not the same.
%[-1,1] in the local coordinate system (master element).
%Integration is carried out on each of the 8 children of the master element.
%Gauss-point (2-pt rule, wts=1)
gpt = [(-1/sqrt(3)), (1/sqrt(3))];

for cNumCoarse=1:8
    cHnMasks = cNumBasedHangingMasks(cNumCoarse);
    ParentCoords = parentOfCoarse3d(cNumCoarse);
    for eType=1:18
        phiCoarse = coarseShapeFunction3DCoeffs(cHnMasks(eType,:),ParentCoords);
        for cNumFine=1:8
            Lmat = zeros(8,8);
            Xfac1 = ((xub(cNumFine) - xlb(cNumFine))/2.0);
            Xfac2 = ((xub(cNumFine) + xlb(cNumFine))/2.0);
            Yfac1 = ((yub(cNumFine) - ylb(cNumFine))/2.0);
            Yfac2 = ((yub(cNumFine) + ylb(cNumFine))/2.0);
            Zfac1 = ((zub(cNumFine) - zlb(cNumFine))/2.0);
            Zfac2 = ((zub(cNumFine) + zlb(cNumFine))/2.0);
            for i=1:8
                cSh_i = phiCoarse(i,:);
                for j=1:8
                    cSh_j = phiCoarse(j,:);
                    %3-D Gauss Quadrature
                    for m=1:2
                        for n=1:2
                            for p=1:2
                                for xIdx=1:3
                                    cShVal_i = evalGradPhi(cSh_i,xIdx,((Xfac1*gpt(m))+Xfac2),((Yfac1*gpt(n))+Yfac2),((Zfac1*gpt(p))+Zfac2));
                                    cShVal_j = evalGradPhi(cSh_j,xIdx,((Xfac1*gpt(m))+Xfac2),((Yfac1*gpt(n))+Yfac2),((Zfac1*gpt(p))+Zfac2));
                                    Lmat(i,j) = Lmat(i,j) + (cShVal_i*cShVal_j);
                                end
                            end
                        end
                    end
                end
            end
            Lmat = (Xfac1*Yfac1*Zfac1)*Lmat;
            fname = ['LmatType1_',int2str(eType),'_',int2str(cNumCoarse),'_',int2str(cNumFine),'.mat'];
            save(fname,'Lmat');
        end
    end
end


%When using the stencil you must multiply by h/2. This factor was removed from
%the integrand. h is the size of the coarse element.


