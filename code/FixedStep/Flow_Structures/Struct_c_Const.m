function [c] = Struct_c_Const(c0)
% Produces a structure containing c and its first derivative in the trivial
% case it is constant

c.f = @(x2) c0*ones(size(x2)) ;
c.df = @(x2) zeros(size(x2)) ;

end