function [x] = SqrtA(x2,cut)
%SqrtA Modifies the sqrt function to allow for arbitrary cut.
%   By shifting and unshifting the argument of x2, the cut can be moved to
%   any straight line from the origin. cut gives the radian angle of the
%   cut

alpha=1i*(cut-pi);
x=exp(alpha/2)*sqrt(exp(-alpha)*x2);


end

