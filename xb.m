function xb_ = xb(xf, xd, F, D)

xb_ = ((F/D)*xf - xd)/(F/D-1);
%xb_ = (xf*F - xd*D)/(F-D);

end