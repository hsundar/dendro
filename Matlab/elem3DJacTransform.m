function jac = elem3DJacTransform(X,Y,Z,cNum)

%1-indexing for cNum
%Local coordinates
syms x y z real;

%Note: It is not symmetric.
%XYZ - Morton Order, C-0 configuration
T = [-1 -1 -1;
    1 -1 -1;
    -1 1 -1;
    1 1 -1;
    -1 -1 1;
    1 -1 1;
    -1 1 1;
    1 1 1;];

%ChildNum Global-Local Transformations
% Global(j) = Local(C(cNum,j)

C = [1 2 3;
     -2 1 3;
     2 -1 3;
     -1 -2 3;
     1 3 -2;
     -3 1 -2;
     3 -1 -2;
     -1 -3 -2];
 
     for i =1:8
         for j =1:3
             if(C(cNum,j) > 0)
                 M(i,C(cNum,j)) = T(i,j);
             else
                 M(i,-C(cNum,j)) = -T(i,j);
             end
         end
     end

     M = [ones(8,1),M];

for(i=1:8)
    M(i,5) = M(i,2)*M(i,3);
    M(i,6) = M(i,3)*M(i,4);
    M(i,7) = M(i,4)*M(i,2);
    M(i,8) = M(i,2)*M(i,3)*M(i,4);
end

xc = inv(M)*X;
yc = inv(M)*Y;
zc = inv(M)*Z;

xPhi = xc(1) + xc(2)*x + xc(3)*y + xc(4)*z + xc(5)*x*y + xc(6)*y*z + xc(7)*z*x + xc(8)*x*y*z;
yPhi = yc(1) + yc(2)*x + yc(3)*y + yc(4)*z + yc(5)*x*y + yc(6)*y*z + yc(7)*z*x + yc(8)*x*y*z;
zPhi = zc(1) + zc(2)*x + zc(3)*y + zc(4)*z + zc(5)*x*y + zc(6)*y*z + zc(7)*z*x + zc(8)*x*y*z;

jac(1,1) = diff(xPhi,x);
jac(1,2) = diff(xPhi,y);
jac(1,3) = diff(xPhi,z);

jac(2,1) = diff(yPhi,x);
jac(2,2) = diff(yPhi,y);
jac(2,3) = diff(yPhi,z);

jac(3,1) = diff(zPhi,x);
jac(3,2) = diff(zPhi,y);
jac(3,3) = diff(zPhi,z);

