function jac = elem2DJacTransform(X,Y)

syms x y real;

M = [1 -1 -1 1;
     1 1 -1 -1;
     1 1 1 1;
     1 -1 1 -1];

xc = inv(M)*X;
yc = inv(M)*Y;

xPhi = xc(1) + xc(2)*x + xc(3)*y + xc(4)*x*y;
yPhi = yc(1) + yc(2)*x + yc(3)*y + yc(4)*x*y;

jac(1,1) = diff(xPhi,x);
jac(1,2) = diff(xPhi,y);
jac(2,1) = diff(yPhi,x);
jac(2,2) = diff(yPhi,y);


