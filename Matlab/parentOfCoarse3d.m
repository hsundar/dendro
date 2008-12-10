function P = parentOfCoarse3d(cNum)

%cNum is the childNumber of the coarse element wr.t its parent.
%The coarse element will always be defined on [-1,1] and its
%parent will be moved around it

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
    P(i,:) = 4*epnuga(MtIdx(i)+1,:);
end

offSet = P(cNum,:) - M(cNum,:);
P = P-repmat(offSet,8,1);



