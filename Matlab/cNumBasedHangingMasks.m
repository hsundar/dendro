function hangingFlags = cNumBasedHangingMasks(cnum)

%1-based indexing for cNums
%There are only 18 valid hanging types depending on the childnumber.
%Values in the same order as in the C++ code (oda.h), if node i (0-based
%indexing) is hanging then (1 << i) is set to 1.
if(cnum == 1)
    hnMasks = [ 0, 4, 2, 6, 16, 20, 18, 22, 14, 30, 84, 86, 94, 50, 54, 62, 118, 126 ];
elseif(cnum == 2)
    hnMasks = [0, 1, 8, 9, 32, 33, 40, 41, 13, 45, 49, 57, 61, 168, 169, 173, 185, 189 ];
elseif(cnum == 3)
    hnMasks = [0, 8, 1, 9, 64, 72, 65, 73, 11, 75, 200, 201, 203, 81, 89, 91, 217, 219 ];
elseif(cnum == 4)
    hnMasks = [0, 2, 4, 6, 128, 130, 132, 134, 7, 135, 162, 166, 167, 196, 198, 199, 230, 231 ];
elseif(cnum == 5)
    hnMasks = [0, 1, 32, 33, 64, 65, 96, 97, 35, 99, 69, 101, 103, 224, 225, 227, 229, 231 ];
elseif(cnum == 6)
    hnMasks = [0, 2, 128, 130, 16, 18, 144, 146, 138, 154, 19, 147, 155, 208, 210, 218, 211, 219 ];
elseif(cnum == 7)
    hnMasks = [0, 4, 16, 20, 128, 132, 144, 148, 21, 149, 140, 156, 157, 176, 180, 181, 188, 189 ];
else
    hnMasks = [0, 8, 64, 72, 32, 40, 96, 104, 76, 108, 42, 106, 110, 112, 120, 124, 122, 126 ];
end

%Now convert to a format suitable for use within matlab
hangingFlags = zeros(18,8);
for i=2:18    
    tmp = dec2bin(hnMasks(i),8);
    for j=1:8                        
        hangingFlags(i,j) = str2double(tmp(9-j));
    end
end
