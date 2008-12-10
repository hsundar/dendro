function gradPhiVal = evalGradPhi(c,xIdx,x,y,z)

if(xIdx == 1)
    gradPhiVal = c(2) + c(5)*y + c(7)*z + c(8)*y*z;
end

if(xIdx == 2)
    gradPhiVal = c(3) + c(5)*x + c(6)*z + c(8)*x*z;
end

if(xIdx == 3)
    gradPhiVal = c(4) + c(6)*y + c(7)*x + c(8)*x*y;
end
