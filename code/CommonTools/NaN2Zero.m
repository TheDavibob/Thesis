function g=NaN2Zero(f,x)
% Makes NaNs into zeros -- run as g=@(x) NaN2Zero(f,x);

% Created 2018-01-24

g=f(x);
g(isnan(g))=0;

end