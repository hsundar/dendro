function phiCoeffs = fineOnCoarseShapeFunction3DCoeffs(cNum,hnMask)

%cNum is the childNumber of the fine element wr.t the coarse element.
%The corresponding coarse element will be defined on [-1,1]

%X Y Z
%The coords of the coarse(parent of the fine)
P = [-1 -1 -1
    1 -1 -1
    -1 1 -1
    1 1 -1
    -1 -1 1
    1 -1 1
    -1 1 1
    1 1 1];

%epsilon, nu, gamma
epnuga = [-1, -1, -1;
    1, -1, -1;
    -1, 1, -1;
    1, 1, -1;
    -1, -1, 1;
    1, -1, 1;
    -1, 1, 1;
    1, 1, 1;
    0, -1, -1;
    0, 1, -1;
    -1, 0, -1;
    1, 0, -1;
    0, -1, 1;
    0, 1, 1;
    -1, 0, 1;
    1, 0, 1;
    -1, -1, 0;
    1, -1, 0;
    -1, 1, 0;
    1, 1, 0;
    0, 0, -1;
    0, 0, 1;
    -1, 0, 0;
    1, 0, 0;
    0, -1, 0;
    0, 1, 0;
    0, 0, 0 ];

if(cNum == 1)
    MtIdx = [0,8,10,20,16,24,22,26];
elseif(cNum == 2)
    MtIdx = [8,1,20,11,24,17,26,23];
elseif(cNum==3)
    MtIdx = [10,20,2,9,22,26,18,25];
elseif(cNum==4)
    MtIdx = [20,11,9,3,26,23,25,19];
elseif(cNum==5)
    MtIdx = [16,24,22,26,4,12,14,21];
elseif(cNum==6)
    MtIdx = [24,17,26,23,12,5,21,15];
elseif(cNum==7)
    MtIdx = [22,26,18,25,14,21,6,13];
else
    MtIdx = [26,23,25,19,21,15,13,7];
end


for i=1:8
    if hnMask(i)
        M(i,:) = P(i,:);
    else
        M(i,:) = epnuga(MtIdx(i)+1,:);
    end
end

M = [ones(8,1),M];

for i=1:8
    M(i,5) = M(i,2)*M(i,3);
    M(i,6) = M(i,3)*M(i,4);
    M(i,7) = M(i,4)*M(i,2);
    M(i,8) = M(i,2)*M(i,3)*M(i,4);
end

phiCoeffs = zeros(8,8);

for i=1:8
    v =zeros(8,1);
    v(i)=1;
    c = inv(M)*v;
    for j =1:8
        phiCoeffs(i,j) = c(j);
    end
end

