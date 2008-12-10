
function nhCnumAndSz = getNeighboringCnumAndSizes()

%1-based numbering

%Notation: element-i is the element whose i-th vertex is the
%grid point of interest.

%For simplicity, element-1 is chosen as the standard reference element
%[0,1]x[0,1]x[0,1] and the grid point of interest has the coords (0,0,0).

%A generic element [0,1]x[0,1]x[0,1]
coordsTemplate = [0 0 0;
    1 0 0;
    0 1 0;
    1 1 0;
    0 0 1;
    1 0 1;
    0 1 1;
    1 1 1];

nhCnumAndSz = cell(8,7);

% Each of 7 surrounding elements could either be of the same size (szOpt =1) or half
% the size (szOpt = 2) or double the size (szOpt = 3) of the reference element.
for myCnum = 1:8
    for elemNum = 2:8
        for szOpt = 1:3
            elemCoords = coordsTemplate;
            if(szOpt == 2)
                elemCoords = 0.5*elemCoords;
            end
            if(szOpt == 3)
                elemCoords = 2*elemCoords;
            end
            offset = elemCoords(elemNum,:);
            elemCoords = elemCoords - repmat(offset,8,1);
            for othCnum = 1:8
                validElem = true;
                for cNumOfBrother = 1: (othCnum -1)
                    brotherCoords = getBrother(othCnum,cNumOfBrother,elemCoords);
                    if( (szOpt == 1) && (cNumOfBrother == myCnum) )
                        if(length(find(coordsTemplate == brotherCoords)) == 24)
                            break;
                        end
                    end
                    weIntersect =  elemsIntersect(brotherCoords, coordsTemplate);
                    if(weIntersect)
                        validElem = false;
                        break;
                    end
                end
                if(validElem)
                    for cNumOfBrother = (othCnum +1) : 8
                        brotherCoords = getBrother(othCnum,cNumOfBrother,elemCoords);
                        if( (szOpt == 1) && (cNumOfBrother == myCnum) )
                            if(length(find(coordsTemplate == brotherCoords)) == 24)
                                break;
                            end
                        end
                        weIntersect =  elemsIntersect(brotherCoords, coordsTemplate);
                        if(weIntersect)
                            validElem = false;
                            break;
                        end
                    end
                    if(validElem)
                       nhCnumAndSz{myCnum,(elemNum-1)}(length(nhCnumAndSz{myCnum,(elemNum-1)})+1) = struct('sz',szOpt,'cNum',othCnum);             
                    end
                end
            end
        end
    end
end


