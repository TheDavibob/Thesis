function [C,dC]=VeitchContour(c,D,W)

% Using the contour given by Veitch and Peake 2008. c is the interesection
% of the real axis, D and W shape parameters. Gives z=C(t)

B=3*W^4; A=4*1i*D*W^3;

C=@(t) t+c-A*t./(B+t.^4);
dC=@(t) 1-A*(B-3*t.^4)./(B+t.^4).^2;

end