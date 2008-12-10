function result = elemsIntersect(elem1, elem2)

result = false;

minX1 = min(elem1(:,1));
minY1 = min(elem1(:,2));
minZ1 = min(elem1(:,3));

maxX1 = max(elem1(:,1));
maxY1 = max(elem1(:,2));
maxZ1 = max(elem1(:,3));

minX2 = min(elem2(:,1));
minY2 = min(elem2(:,2));
minZ2 = min(elem2(:,3));

maxX2 = max(elem2(:,1));
maxY2 = max(elem2(:,2));
maxZ2 = max(elem2(:,3));

if ( (minX1 >= minX2) && (minX1 < maxX2) && (minY1 >= minY2) && (minY1 < maxY2) && (minZ1 >= minZ2) && (minZ1 < maxZ2) ) 
    %anchor of elem1 is inside elem2
    result = true;
end

if ( (minX2 >= minX1) && (minX2 < maxX1) && (minY2 >= minY1) && (minY2 < maxY1) && (minZ2 >= minZ1) && (minZ2 < maxZ1) )      
    %anchor of elem2 is inside elem1
    result = true;
end
