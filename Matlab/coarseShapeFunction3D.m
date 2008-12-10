function phi = coarseShapeFunction3D(hnMask,ParentCoords)

%cNum is the childNumber of the coarse element wr.t its parent.
%The coarse element will always be defined on [-1,1] and its
%parent will be moved around it

syms x y z real;

%X Y Z
%The coords of the coarse
M = [-1 -1 -1
    1 -1 -1
    -1 1 -1
    1 1 -1
    -1 -1 1
    1 -1 1
    -1 1 1
    1 1 1];

for(i=1:8)
    if hnMask(i)
        M(i,:) = ParentCoords(i,:);
    end
end

M = [ones(8,1),M];

for(i=1:8)
    M(i,5) = M(i,2)*M(i,3);
    M(i,6) = M(i,3)*M(i,4);
    M(i,7) = M(i,4)*M(i,2);
    M(i,8) = M(i,2)*M(i,3)*M(i,4);
end

for(i=1:8)
    v =zeros(8,1);
    v(i)=1;
    c = inv(M)*v;
    phi(i,1) = c(1) + c(2)*x + c(3)*y + c(4)*z + c(5)*x*y + c(6)*y*z + c(7)*z*x + c(8)*x*y*z;
end

