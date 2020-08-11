function [y] = Two_Solutions_Complex(alpha,beta,param)

% Given complex anonymous functions of z, alpha and beta, solves the
% differential equation
% z^2 y'' + alpha z y' + beta y = 0
% giving two solutions.

% The two solutions are structured functions, with a variety of fields.

% param governs the series expansion. It has
% .N fields (the number of terms of the series),
% .r field (the radius of convergence chosen for expanding  alpha and beta,
% ideally greater than one if possible)
% .m field (the number of points for computing coefficients, m >> N).
% .cut is the direction of the branch cut when logs and fractional powers
% are involved.

[a] = TaylorCoeff_Complex(alpha,param.N,param.r,param.m);
[b] = TaylorCoeff_Complex(beta,param.N,param.r,param.m);

[c,d,C,s1,s2,I] = Two_Series(param.N,a,b);

% Clean up the functions a little bit.
%c(abs(c)<1e-5) = 0;
%d(abs(d)<1e-5) = 0;
C(abs(C)<1e-5) = 0;
s1(abs(s1)<1e-5) = 0;
s2(abs(s2)<1e-5) = 0;

f1 = Taylor_Series(c);
f2 = Taylor_Series(d);

nf1 = Taylor_Series(c.*(0:param.N));
nf2 = Taylor_Series(d.*(0:param.N));


y.y1 = @(z) PowerA(z,s1,param.cut).*f1(z);
y.dy1 = @(z) PowerA(z,s1,param.cut)./z.*(s1*f1(z)+nf1(z));

y.y2 = @(z) C*LogA(z,param.cut) .*y.y1(z) + PowerA(z,s2,param.cut).*f2(z);
y.dy2 = @(z) C*y.y1(z)./z + C*y.dy1(z).*LogA(z,param.cut) + PowerA(z,s2,param.cut)./z.*(s2*f2(z)+nf2(z));

y.s1 = s1; % The larger root of indicial equation
y.s2 = s2; % The smaller
y.C = C; % The constant preceding the log term, possibly zero
y.c = c; % The coefficients of y1
y.d = d; % The coefficients of y2
y.I = I; % The type of solution (behaviour of roots of indicial equation)

y.a = a;
y.b = b;

end