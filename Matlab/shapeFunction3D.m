function [phi, gradPhi] = shapeFunction3D(dzdxOn, dydzOn, dxdyOn, dzOn, dxOn, dyOn,varargin)

%Only for C-0 configuration

syms x y z real;

if(dzdxOn)
    dxp1=2;
    dzp1=2;
    dxOn=1;
    dzOn=1;
else
    dxp1=0;
    dzp1=0;
end

if(dydzOn)
    dyp2=2;
    dzp2=2;
    dyOn=1;
    dzOn=1;
else
    dyp2=0;
    dzp2=0;
end

if(dxdyOn)
    dxp3=2;
    dyp3=2;
    dxOn=1;
    dyOn=1;
else
    dxp3=0;
    dyp3=0;
end

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

if(dzOn)
    dz = 2;
else
    dz = 0;
end
%X Y Z
M = [-1 -1 -1;
    (1 + dx) -1 -1;
    -1 (1+dy) -1;
    (1+dxp3) (1+dyp3) -1;
    -1 -1 (1+dz);
    (1+dxp1) -1 (1+dzp1);
    -1 (1+dyp2) (1+dzp2);
    1 1 1];

M = [ones(8,1),M];

for(i=1:8)
    M(i,5) = M(i,2)*M(i,3);
    M(i,6) = M(i,3)*M(i,4);
    M(i,7) = M(i,4)*M(i,2);
    M(i,8) = M(i,2)*M(i,3)*M(i,4);
end

if(nargin == 7)
    display('Accepted invJac')
    invJac = varargin{1};
end

for(i=1:8)
    v =zeros(8,1);
    v(i)=1;
    c = inv(M)*v;
    phi(i,1) = c(1) + c(2)*x + c(3)*y + c(4)*z + c(5)*x*y + c(6)*y*z + c(7)*z*x + c(8)*x*y*z;
    if(nargout == 2)
        gradPhi(i,1) = diff(phi(i),x)*invJac(1,1) + diff(phi(i),y)*invJac(2,1) + diff(phi(i),z)*invJac(3,1);
        gradPhi(i,2) = diff(phi(i),x)*invJac(1,2) + diff(phi(i),y)*invJac(2,2) + diff(phi(i),z)*invJac(3,2);
        gradPhi(i,3) = diff(phi(i),x)*invJac(1,3) + diff(phi(i),y)*invJac(2,3) + diff(phi(i),z)*invJac(3,3);
    end
end

if(nargout == 2)
    display('Computed gradPhi')
end
