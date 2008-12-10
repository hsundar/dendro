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

%The integration is on the domain
%[-1,1] in the local coordinate system (master element).

%Gauss-point (2-pt rule, wts=1)
gpt = [(-1/sqrt(3)), (1/sqrt(3))];

for cNum=1:8
    HnMasks = cNumBasedHangingMasks(cNum);
    ParentCoords = parentOfCoarse3d(cNum);
    for eType=1:18
        phi = coarseShapeFunction3DCoeffs(HnMasks(eType,:),ParentCoords);
        Mmat = zeros(8,8);
        for i=1:8
            for j=1:8
                %The factor of detJac has been removed since it is
                %a constant for this integration and can be
                %introduced while using the matrices
                sh_i = phi(i,:);
                sh_j = phi(j,:);

                %3-D Gauss Quadrature
                Mmat(i,j) = 0;
                for m=1:2
                    for n=1:2
                        for p=1:2
                            ShVal_i = evalPhi(sh_i,gpt(m),gpt(n),gpt(p));
                            ShVal_j = evalPhi(sh_j,gpt(m),gpt(n),gpt(p));
                            Mmat(i,j) = Mmat(i,j) + (ShVal_i*ShVal_j);
                        end
                    end
                end
            end
        end

        fname = ['MassMatrix_',int2str(cNum),'_',int2str(eType),'.mat'];
        save(fname,'Mmat');
    end
end

%When using the stencil it must be multiplied by hc^3
