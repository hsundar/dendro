clc;clear;
%Morton Ordered Element. Child Num 1, X grows first Y second Z last
%X right, Y Back Z top
%1-based indexing
%J1 = gradPhi_i*gradPhi_j
%J2 = Phi_i*Phi_j
%dzdx dydz dxdy dz dx dy
%6 7 4 5 2 3
v= [0 0 0 0 0 0;
    0 0 0 0 0 1;
    0 0 0 0 1 0;
    0 0 0 0 1 1;
    0 0 0 1 0 0;
    0 0 0 1 0 1;
    0 0 0 1 1 0;
    0 0 0 1 1 1;
    0 0 1 0 1 1;
    0 0 1 1 1 1;
    0 1 0 1 0 1;
    0 1 0 1 1 1;
    0 1 1 1 1 1;
    1 0 0 1 1 0;
    1 0 0 1 1 1;
    1 0 1 1 1 1;
    1 1 0 1 1 1;
    1 1 1 1 1 1];

[numCases,dum] = size(v);

syms x1 y1 z1 hx hy hz x y z real
%Morton Order: X first Y second Z last
%Global Coordinates
XYZ = [x1 y1 z1;
    (x1+hx) y1 z1;
    x1 (y1+hy) z1;
    (x1+hx) (y1+hy) z1;
    x1 y1 (z1+hz);
    (x1+hx) y1 (z1+hz);
    x1 (y1+hy) (z1+hz)
    (x1+hx) (y1+hy) (z1+hz) ];

% for(cNum=1:8)
%      jac = elem3DJacTransform(XYZ(:,1),XYZ(:,2),XYZ(:,3),cNum);
%      detJac = simplify(det(jac))
% end


for(cNum=1:8)
    %The only thing that depends on the childnumber is the global-to-local
    %transformations...
    jac = elem3DJacTransform(XYZ(:,1),XYZ(:,2),XYZ(:,3),cNum);
    detJac = simplify(det(jac));
    invJac = inv(jac);
    syms h real;
    for(i=1:numCases)
        %This is computed in local numbering
        [phi, gradPhi] = shapeFunction3D(v(i,1), v(i,2), v(i,3), v(i,4), v(i,5), v(i,6),invJac);
        %J1 and J2 will be stored in local numbering.
        for(j=1:8)
            for(k=1:8)
                integrand1 = dotProduct(gradPhi(j,:),gradPhi(k,:))*detJac;
                integrand2 = phi(j)*phi(k)*detJac;
                J1(j,k,i,cNum) = int(int(int(integrand1,x,-1,1),y,-1,1),z,-1,1);
                J2(j,k,i,cNum) = int(int(int(integrand2,x,-1,1),y,-1,1),z,-1,1);
                J1(j,k,i,cNum) = subs(J1(j,k,i,cNum),{hx,hy,hz},{h,h,h});
                J2(j,k,i,cNum) = subs(J2(j,k,i,cNum),{hx,hy,hz},{h,h,h});
                J1(j,k,i,cNum) = simplify(J1(j,k,i,cNum));
                J2(j,k,i,cNum) = simplify(J2(j,k,i,cNum));
            end
        end
    end
end

