function [Fp,Fm] = Fact_Additive_Finite(F,contour,prec)
% Additively factorised function F, which decays at infinity, using a
% substitution. Otherwise exactly as before.

% contour contains Cp,m and dCp,m fields, as per usual (combined)
% prec contains an N field (number of steps) and max and min fields (limit of
% integration).

% Substitution:
T = @(s) s./(1-s.^2) ; % Not quite the Veitch/Reinstra choice but valid in both directions
dT = @(s) (1+s.^2)./((1-s.^2).^2);
% This maps (-inf,inf) to the finite range (-1,1) on which can do normal
% quadrature. Thereby don't use most of prec.


[S,Ws]=lgwt(prec.N,-1,1);

S = S(:).';
s = S(end:-1:1);
Ws = Ws(:).';
ws = Ws(end:-1:1);

t = T(s);
wt = ws.*dT(s);

zp = contour.Cp(t);
dzp = wt.*contour.dCp(t);

zm = contour.Cm(t);
dzm = wt.*contour.dCm(t);

Cauchy_Kernelp = @(k) 1/(2*pi*1i) .* 1./(zp(:) - k(:).') ;
Cauchy_Kernelm = @(k) 1/(2*pi*1i) .* 1./(zm(:) - k(:).') ;

F_nan = @(k1) NoNan(F,k1) ; % Removes NaN singularities

Fp = @(k) F_nan(zp).*dzp * Cauchy_Kernelp(k) ;
Fm = @(k) -F_nan(zm).*dzm * Cauchy_Kernelm(k) ;

end

function F_nan = NoNan(F,k1)
    F_nan = F(k1);
    F_nan(isnan(F_nan)) = 0;
end