function res = dotProduct(a,b)

res = (a(1)*b(1));
for i=2:length(a)
 res = res + (a(i)*b(i));
end
