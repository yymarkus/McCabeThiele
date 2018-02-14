function NTUstr = NTUstr(xd, xf, xtrans, xb, R, F, D)

x = xb;
%if xtrans xtrans < xb;
%    xtrans;
%end

range = xtrans - xb;
stepsize = 0.0001;
steps = range/stepsize;
sum = 0;
for i=1:steps
    if i == 1 || i == steps
        multip = 1;
    else
        multip = 2;
    end
    denominator =  (R+1)*fstar(x) - (R+F/D)*x + F/D*xf - xd;
    %denominator =  (R+1)*fstar(x) - (R+F/D)*x + xb*(F/D-1);
    Integrand = (R + F/D)/denominator;
    sum = sum + Integrand*multip;
    x = x + stepsize;
end
NTUstr = sum*stepsize/2;
end