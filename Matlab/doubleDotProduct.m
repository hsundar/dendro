%Tensor Double Dot Product = AijBji (sum on i and j)
function res = doubleDotProduct(A,B)

[rA, cA] = size(A);
[rB, cB] = size(B);


if(rA ~= cB)
    error('Wrong Dimensions')
end

if(rB ~= cA)
    error('Wrong Dimensions')
end

res = 0;
for(i=1:rA)
    for(j=1:cA)
        res = res + (A(i,j)*B(j,i));
    end
end