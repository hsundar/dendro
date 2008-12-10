clc;clear;
%X right, Y Back Z top
%1-based indexing

% syms x1 y1 z1 hx hy hz x y z hc real;
% %Morton Order: X first Y second Z last
% %Global Coordinates
% XYZ = [x1 y1 z1;
%     (x1+hx) y1 z1;
%     x1 (y1+hy) z1;
%     (x1+hx) (y1+hy) z1;
%     x1 y1 (z1+hz);
%     (x1+hx) y1 (z1+hz);
%     x1 (y1+hy) (z1+hz)
%     (x1+hx) (y1+hy) (z1+hz) ];
%
% jac = elem3DJacTransform(XYZ(:,1),XYZ(:,2),XYZ(:,3),1);
% detJac = simplify(det(jac));

%detJac is a constant = hc^3/8

%Type-2 Coarse and Fine are the same. The integration is on the domain
%[-1,1] in the local coordinate system (master element).

%Gauss-point (2-pt rule, wts=1)
gpt = [(-1/sqrt(3)), (1/sqrt(3))];

for cNum=1:8
    HnMasks = cNumBasedHangingMasks(cNum);
    ParentCoords = parentOfCoarse3d(cNum);
    for coarseType=1:18
        phiCoarse = coarseShapeFunction3DCoeffs(HnMasks(coarseType,:),ParentCoords);
        for fineType=1:18
            phiFine = coarseShapeFunction3DCoeffs(HnMasks(fineType,:),ParentCoords);

            Rmat = zeros(8,8);
            for i=1:8
                for j=1:8
                    %The factor of detJac has been removed since it is
                    %a constant for this integration and can be
                    %introduced while using the matrices
                    cSh = phiCoarse(i,:);
                    fSh = phiFine(j,:);

                    %3-D Gauss Quadrature
                    Rmat(i,j) = 0;
                    for m=1:2
                        for n=1:2
                            for p=1:2
                                cShVal = evalPhi(cSh,gpt(m),gpt(n),gpt(p));
                                fShVal = evalPhi(fSh,gpt(m),gpt(n),gpt(p));
                                Rmat(i,j) = Rmat(i,j) + (cShVal*fShVal);
                            end
                        end
                    end
                end
            end

            fname = ['L2basedRmatType2_',int2str(cNum),'_',int2str(coarseType),'_',int2str(fineType),'.mat'];
            save(fname,'Rmat');
        end
    end
end

%When using the stencil it must be multiplied by hc^3
