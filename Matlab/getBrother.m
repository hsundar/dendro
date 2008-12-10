function brother = getBrother(myCnum, cNumOfBrother, myVertices)

if(myCnum == cNumOfBrother)
    error('myCnum and cNumOfBrother cannot be the same.');
end

if(myCnum < 1 || myCnum > 8)
    error('myCnum must be >=1 and <= 8.');
end

if(cNumOfBrother < 1 || cNumOfBrother > 8)
    error('cNumOfBrother must be >=1 and <= 8.');
end

myParent = getParent(myCnum,myVertices);

brother = myVertices;
offset = brother(cNumOfBrother,:) - myParent(cNumOfBrother,:);
brother = brother - repmat(offset,8,1);
