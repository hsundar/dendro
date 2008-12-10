function res = psi(num)
%x1,y1 is the anchor of the element and hx and hy are its dimensions.
syms x y x1 y1 hx hy real

a = 2/hx;
c = 2/hy;
b = -(1+ (a*x1));
d = -(1 + (c*y1));

ep = ((a*x) + b);
eta = ((c*y) + d);

if(num == 1)
    res = (1/4)*(1-ep)*(1-eta);
else if(num == 2)
        res = (1/4)*(1+ep)*(1-eta);
    else if(num == 3)
            res = (1/4)*(1+ep)*(1+eta);
        else
            res = (1/4)*(1-ep)*(1+eta);
        end
    end
end


