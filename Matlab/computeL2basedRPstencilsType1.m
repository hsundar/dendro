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

xlb = [-1 0 -1 0 -1 0 -1 0];
xub = [0 1 0 1 0 1 0 1];
ylb = [-1 -1 0 0 -1 -1 0 0];
yub = [0 0 1 1 0 0 1 1];
zlb = [-1 -1 -1 -1 0 0 0 0];
zub = [0 0 0 0 1 1 1 1];

%Gauss-point (2-pt rule, wts=1)
gpt = [(-1/sqrt(3)), (1/sqrt(3))];

%Type-1 Coarse and Fine are not the same.
for cNumCoarse=1:8
    cHnMasks = cNumBasedHangingMasks(cNumCoarse);
    ParentCoords = parentOfCoarse3d(cNumCoarse);
    for coarseType=1:18
        phiCoarse = coarseShapeFunction3DCoeffs(cHnMasks(coarseType,:),ParentCoords);
        for cNumFine=1:8
            fHnMasks = cNumBasedHangingMasks(cNumFine);
            Xfac1 = ((xub(cNumFine) - xlb(cNumFine))/2);
            Xfac2 = ((xub(cNumFine) + xlb(cNumFine))/2);
            Yfac1 = ((yub(cNumFine) - ylb(cNumFine))/2);
            Yfac2 = ((yub(cNumFine) + ylb(cNumFine))/2);
            Zfac1 = ((zub(cNumFine) - zlb(cNumFine))/2);
            Zfac2 = ((zub(cNumFine) + zlb(cNumFine))/2);
            for fineType=1:18
                phiFine = fineOnCoarseShapeFunction3DCoeffs(cNumFine,fHnMasks(fineType,:));
                
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
                                     cShVal = evalPhi(cSh,((Xfac1*gpt(m))+Xfac2),((Yfac1*gpt(n))+Yfac2),((Zfac1*gpt(p))+Zfac2));
                                     fShVal = evalPhi(fSh,((Xfac1*gpt(m))+Xfac2),((Yfac1*gpt(n))+Yfac2),((Zfac1*gpt(p))+Zfac2));
                                     Rmat(i,j) = Rmat(i,j) + (cShVal*fShVal);
                                 end
                             end
                         end                         
                    end
                end         
                Rmat = (Xfac1*Yfac1*Zfac1)*Rmat;
                fname = ['L2basedRmatType1_',int2str(coarseType),'_',int2str(cNumCoarse),...
                    '_',int2str(fineType),'_',int2str(cNumFine),'.mat'];
                save(fname,'Rmat');
            end
        end
    end
end

%When using the stencil it must be multiplied by hc^3
