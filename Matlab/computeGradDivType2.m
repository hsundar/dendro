clc;clear;
%X right, Y Back Z top
%1-based indexing

%The global coordinate system and the local coordinate system are aligned. So
%the jacobian is only a diagonal matrix with entries h/2.
%Type-2 Coarse and Fine are the same. The integration is on the domain
%[-1,1] in the local coordinate system (master element).

%Gauss-point (2-pt rule, wts=1)
gpt = [(-1/sqrt(3)), (1/sqrt(3))];

for cNum=1:8
    HnMasks = cNumBasedHangingMasks(cNum);
    ParentCoords = parentOfCoarse3d(cNum);
    for eType=1:18
        phiCoarse = coarseShapeFunction3DCoeffs(HnMasks(eType,:),ParentCoords);
        GDmat = zeros(24,24);
        for i = 1:8
            cSh_i = phiCoarse(i,:);
            for j = 1:8
                cSh_j = phiCoarse(j,:);
                %Integrate using quadratures
                %3-D Gauss Quadrature
                for m=1:2
                    for n=1:2
                        for p=1:2
                            for dof1 = 1:3
                                for dof2 = 1:3
                                    cShVal_i = evalGradPhi(cSh_i,dof1,gpt(m),gpt(n),gpt(p));
                                    cShVal_j = evalGradPhi(cSh_j,dof2,gpt(m),gpt(n),gpt(p));
                                    GDmat((3*(i-1)+dof1),(3*(j-1)+dof2)) = ...
                                    GDmat((3*(i-1)+dof1),(3*(j-1)+dof2)) + (cShVal_i*cShVal_j);
                                end
                            end
                        end
                    end
                end
            end
        end
        fname = ['GDmatType2_',int2str(cNum),'_',int2str(eType),'.mat'];
        save(fname,'GDmat');
    end
end

%When using the stencil you must multiply by h/2. This factor was removed from
%the integrand


