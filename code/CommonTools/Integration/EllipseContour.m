function [C,dC]=EllipseContour(H,W)
% an ellipse from 0 to W of height h, otherwise a flat line.

b=2*H/W;

ge=@(t) sqrt(H^2-b^2*(t-W/2).^2);
dge=@(t) -b^2*(t-W/2)./ge(t);

g=@(t) ge(t).*heaviside(W-t);
dg=@(t) dge(t).*heaviside(W-t);

Cp=@(t) t+1i*g(t);
dCp=@(t) 1+dg(t);

Cm=@(t) t;
dCm=@(t) 1;

C=@(t) Cp(t).*heaviside(t)+Cm(t).*heaviside(-t);
dC=@(t) dCp(t).*heaviside(t)+dCm(t).*heaviside(-t);

end