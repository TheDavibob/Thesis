function [U] = Struct_U_Cubic_B(U_points,delta_points)
% Produces a structure containing U and its derivatives for a simple
% cubic profile. The profile passes through U0, U1 and Uinf, with vanishing
% shear at deltainf.

U0 = U_points(1);
U1 = U_points(2);
Uinf = U_points(3);

delta0 = delta_points(1);
delta1 = delta_points(2);
deltainf = delta_points(3);

B = (U1 - Uinf)./((delta1 - deltainf).^2*(delta1-delta0)) - ...
    (U0 - Uinf)./((delta0 - deltainf).^2*(delta1-delta0));
A = (U1 - Uinf)./((delta1 - deltainf).^2) - ...
    B*(delta1 - deltainf);

U.f = @(x2) Uinf + A*(x2 -deltainf).^2 + B*(x2 - deltainf).^3;
U.df = @(x2) 2*A*(x2 -deltainf) + 3*B*(x2 - deltainf).^2;
U.d2f = @(x2) 2*A + 6*B*(x2 - deltainf);

end