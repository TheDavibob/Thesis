function [U] = Struct_U_Const(U0)
% Produces a structure containing U and its derivatives for a constant
% flow profile

% U(0) = U0, U(deltainf)=Uinf, etc.

U.f = @(x2) U0*ones(size(x2)) ;
U.df = @(x2) zeros(size(x2));
U.d2f = @(x2) zeros(size(x2));

end