function NTUenr = NTUenr(xd, xtrans, R)

x = xtrans;

range = xd - xtrans;
stepsize = 0.0001;
steps = range/stepsize;
sum = 0;
%figure
for i=1:steps
    if i == 1 || i == steps
        multip = 1;
    else
        multip = 2;
    end
    Denominator = (R+1)*fstar(x) - R*x - xd;
    Integrand = R/Denominator;
    sum = sum + Integrand*multip;
    x = x + stepsize;
    %plot([i], [Denominator], 'og')
    %hold on
end
NTUenr = sum*stepsize/2;
end
