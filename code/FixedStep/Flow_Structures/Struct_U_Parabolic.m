function [U] = Struct_U_Parabolic(U0,Uinf,deltainf)
% Produces a structure containing U and its derivatives for a simple
% parabolic profile, with continuous shear at the edge of the boundary
% layer.

% U(0) = U0, U(deltainf)=Uinf, etc.

U.f = @(x2) Uinf - (Uinf-U0)/(deltainf^2)*(x2-deltainf).^2 ;
U.df = @(x2) -2 * (Uinf-U0)/(deltainf^2)*(x2-deltainf);
U.d2f = @(x2) - 2 * (Uinf-U0)/deltainf^2*ones(size(x2));

end