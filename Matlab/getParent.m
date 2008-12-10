function myParent = getParent(myCnum, myVertices)

if(myCnum < 1 || myCnum > 8)
    error('myCnum must be >=1 and <= 8.');
end

myParent = 2*myVertices;
offset = myParent(myCnum,:) - myVertices(myCnum,:);
myParent = myParent - repmat(offset,8,1);
