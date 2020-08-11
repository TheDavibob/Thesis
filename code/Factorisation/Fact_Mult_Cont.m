function [Gplus,Gminus] = Fact_Mult_Cont(G,contour,prec,branch_cut)
% Multiplicatively factorises G (as per Fact_Mult) and analytically
% continues to the other half plane.

% Requires the contour to be such that Re(C(t)) = t. This could be
% extended, but I'm not going to.

[Gp,Gm] = Fact_Mult(G,contour,prec,branch_cut);

Ip = @(k1) heaviside(imag(k1) - imag(contour.Cp(real(k1))));
Im = @(k1) heaviside(imag(k1) - imag(contour.Cm(real(k1))));

Gplus = @(k1) Gp(k1).*Ip(k1) + (G(k1)./Gm(k1)).*(1-Ip(k1));
Gminus = @(k1) (G(k1)./Gp(k1)).*Im(k1) + Gm(k1).*(1-Im(k1));

% i.e. Gplus and Gminus are Gp and Gm on the respective well-defined bits,
% and extended for the less well-defined bits.