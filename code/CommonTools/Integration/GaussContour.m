function [C,dC]=GaussContour(Cinf,H,W)

% a Gaussian contour passing through 0, with C(-inf)=C(inf)=Cinf, maximum H
% and width between zeros W (second zero at W+0i).

% Cinf<0

ginf=Cinf;
g0=H;
t0=W/2;
sigma2=log((ginf-g0)/ginf)/(t0^2);

g = @(t) ginf+(g0-ginf)*exp(-sigma2*(t-t0).^2);
dg = @(t) -2*sigma2*(g0-ginf)*(t-t0).*exp(-sigma2*(t-t0).^2);

C=@(t) t+1i*g(t);
dC=@(t) 1+1i*dg(t);

end