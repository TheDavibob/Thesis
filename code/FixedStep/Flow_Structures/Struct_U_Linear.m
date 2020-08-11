function [U] = Struct_U_Linear(U0,Uinf,deltainf)
% Produces a structure containing U and its derivatives for a simple linear
% profile.

% U(0) = U0, U(deltainf)=Uinf, etc.

U.f = @(x2) U0+(Uinf-U0).*x2/deltainf ;
U.df = @(x2) (Uinf-U0)/deltainf * ones(size(x2));
U.d2f = @(x2) zeros(size(x2));

end