function y = yStripping(x, xf, xd, R, F, D, xb, B)

%y = ((R+(F/D))/(R+1))*x - ((F/D)*xf - xd)/(R+1);
%L = R*D
y = ((R+(F/D))/(R+1))*x - xb*(F/D-1)/(R+1);
%y = L/(L-B)*x - B/(L-B)*xb


end
