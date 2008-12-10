function [phi, gradPhi] = shapeFunction2D(dxOn,dyOn,invJac)

syms x y real;

if(dxOn)
    dx = 2;
else
    dx = 0;
end

if(dyOn)
    dy = 2;
else
    dy = 0;
end

M=[-1 -1;
    (1+dx) -1;
    1 1;
    -1 (1+dy)];

M = [ones(4,1),M];

for(i=1:4)
    M(i,4) = M(i,2)*M(i,3);
end


for(i=1:4)
    v =zeros(4,1);
    v(i)=1;
    c = inv(M)*v;
    phi(i,1) = c(1) + c(2)*x + c(3)*y + c(4)*x*y;
    gradPhi(i,1) = diff(phi(i),x)*invJac(1,1) + diff(phi(i),y)*invJac(2,1);
    gradPhi(i,2) = diff(phi(i),x)*invJac(1,2) + diff(phi(i),y)*invJac(2,2);
end
