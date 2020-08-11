function [kappa1,dkappa1,y2] = kappa1_Quad(omega,U,delta,N);
% Given boundary-layer profile, computes N quadrature points for kappa1
% (and weights). Then computes the relevant y2 for this setup.

a = omega./U.f(0);
b = omega./U.f(delta);

[kappa1,dkappa1] = lgwt(N,a,b);
kappa1 = kappa1(end:-1:1).';
dkappa1 = dkappa1(end:-1:1).';

% Does a quick interpolation to compute the inverse.
y2_prior = linspace(0,delta,10*N);
kappa1_prior = omega./U.f(y2_prior);

y2 = interp1(kappa1_prior,y2_prior,kappa1);

end