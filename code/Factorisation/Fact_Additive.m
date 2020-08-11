function [Fp,Fm] = Fact_Additive(F,contour,prec)
% Additively factorised function F, which decays at infinity, using a
% Gauss-Legendre-type thing

% contour contains Cp,m and dCp,m fields, as per usual (combined)
% prec contains an N field (number of steps) and max and min fields (limit of
% integration).

[T,Wt]=lgwt(prec.N,prec.min,prec.max);

T = T(:).';
t = T(end:-1:1);
Wt = Wt(:).';
wt = Wt(end:-1:1);

zp = contour.Cp(t);
dzp = wt.*contour.dCp(t);

zm = contour.Cm(t);
dzm = wt.*contour.dCm(t);

Cauchy_Kernelp = @(k) 1/(2*pi*1i) .* 1./(zp(:) - k(:).') ;
Cauchy_Kernelm = @(k) 1/(2*pi*1i) .* 1./(zm(:) - k(:).') ;

Fp = @(k) F(zp).*dzp * Cauchy_Kernelp(k) ;
Fm = @(k) -F(zm).*dzm * Cauchy_Kernelm(k) ;

end