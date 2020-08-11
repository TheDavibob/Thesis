function [U] = Struct_U_Cubic(U0,sigma0)
% Produces a structure containing U and its derivatives for a simple
% cubic profile, with slip velocity and shear prescribed, and smooth at
% edge of boundary layer (delta = 1, Uinf = 1).

% U0 = U(0), sigma0 = U'(0);

% Profile monotonic provided 3U0 + sigma0 < 3 (and sigma0 > 0)
% Inflection free provided 3*U0 + 2*sigma0 > 3, and the above

A = U0;
B = sigma0;
C = 3 - 3*U0-2*sigma0;
D = -2 + 2*U0 + sigma0;

U.f = @(x2) A + B*x2 + C*x2.^2 + D*x2.^3;
U.df = @(x2) B + 2*C*x2 + 3*D*x2.^2;
U.d2f = @(x2) 2*C + 6*D*x2;
U.d3f = @(x2) 6*D;

end