function res = gradPsi(num)

%x1,y1 is the anchor of the element and hx and hy are its dimensions.

syms x y x1 y1 hx hy real

a = 2/hx;
c = 2/hy;
b = -(1+ (a*x1));
d = -(1 + (c*y1));

ep = ((a*x) + b);
eta = ((c*y) + d);

if(num == 1)
    psi = (1/4)*(1-ep)*(1-eta);
else if(num == 2)
        psi = (1/4)*(1+ep)*(1-eta);
    else if(num == 3)
            psi = (1/4)*(1+ep)*(1+eta);
        else
            psi = (1/4)*(1-ep)*(1+eta);
        end
    end
end

res(1) = diff(psi,x);
res(2) = diff(psi,y);
