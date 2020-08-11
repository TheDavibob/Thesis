function [U] = Struct_U_Quartic(U_points,delta_points)
% Produces a structure containing U and its derivatives for a simple
% quartic profile. The profile passes through U0, U1 and Uinf, with vanishing
% shear at deltainf, and should always be decreasing?

U0 = U_points(1);
U1 = U_points(2);
Uinf = U_points(3);

delta0 = delta_points(1);
delta1 = delta_points(2);
deltainf = delta_points(3);

B = (U1 - Uinf)./((delta1 - deltainf).^2*((delta1-deltainf).^2-(delta0 - deltainf).^2)) - ...
    (U0 - Uinf)./((delta0 - deltainf).^2*((delta1-deltainf).^2-(delta0 - deltainf).^2));
A = (U1 - Uinf)./((delta1 - deltainf).^2) - ...
    B*(delta1 - deltainf).^2;

U.f = @(x2) Uinf + A*(x2 -deltainf).^2 + B*(x2 - deltainf).^4;
U.df = @(x2) 2*A*(x2 -deltainf) + 4*B*(x2 - deltainf).^3;
U.d2f = @(x2) 2*A + 12*B*(x2 - deltainf).^2;

end