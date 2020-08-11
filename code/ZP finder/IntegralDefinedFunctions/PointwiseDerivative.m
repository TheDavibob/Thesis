function [df] = PointwiseDerivative(f,z)
%PointwiseDerivative Simple two-point derivative for a function defined
%pointwise.

% h is real and fixed wlog.

h=1e-10;

df=(f(z+h)-f(z-h))/(2*h);
end

