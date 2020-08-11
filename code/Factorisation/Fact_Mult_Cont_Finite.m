function [Gplus,Gminus] = Fact_Mult_Cont_Finite(G,contour,prec,branch_cut)
% Multiplicatively factorises G (as per Fact_Mult) and analytically
% continues to the other half plane.

% Requires the contour to be such that Re(C(t)) = t + const. This could be
% extended, but I'm not going to.

shift_p = contour.Cp(0); % This should be real.
shift_m = contour.Cm(0);

[Gp,Gm] = Fact_Mult_Finite(G,contour,prec,branch_cut);

Ip = @(k1) heaviside(imag(k1) - imag(contour.Cp(real(k1)-shift_p)));
Im = @(k1) heaviside(imag(k1) - imag(contour.Cm(real(k1)-shift_m)));

Gplus = @(k1) Gp(k1).*Ip(k1) + (G(k1)./Gm(k1)).*(1-Ip(k1));
Gminus = @(k1) (G(k1)./Gp(k1)).*Im(k1) + Gm(k1).*(1-Im(k1));

% i.e. Gplus and Gminus are Gp and Gm on the respective well-defined bits,
% and extended for the less well-defined bits.