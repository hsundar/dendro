function res = psiCoarseOnFine(coarseNum,elemNum)

%x1,y1 is the anchor of the fine element and hx and hy are the dimensions
%of the fine element.
%The coarse elements dimensions are (2*hx) and (2*hy)

syms x y x1 y1 hx hy real

hxc  = (2*hx);
hyc = (2*hy);

if(elemNum == 1)
    x1c = x1;
    y1c = y1;
else if(elemNum == 2)
        x1c = (x1 - hx);
        y1c = y1;
    else if(elemNum == 3)
            x1c = (x1 - hx);
            y1c = (y1 - hy);
        else
            x1c = x1;
            y1c = (y1 - hy);
        end
    end
end

a = 2/hxc;
c = 2/hyc;
b = -(1 + (a*x1c));
d = -(1 + (c*y1c));

ep = ((a*x) + b);
eta = ((c*y) + d);

if(coarseNum == 1)
    res = ((1/4)*(1-ep)*(1-eta));
else if(coarseNum == 2)
        res = ((1/4)*(1+ep)*(1-eta));
    else if(coarseNum == 3)
            res = ((1/4)*(1+ep)*(1+eta));
        else
            res = ((1/4)*(1-ep)*(1+eta));
        end
    end
end

