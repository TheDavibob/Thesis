function [x] = LogA(x2,cut)
%LogA Modifies the log function to allow for arbitrary cut.
%   By shifting and unshifting the argument of x2, the cut can be moved to
%   any straight line from the origin. cut gives the radian angle of the
%   cut

% Note care is needed with the angle of this cut - -pi/2 and 3pi/2 are not
% the same (they will differ by 2ipi).

alpha=1i*(cut-pi);
x=log(exp(-alpha)*x2)+alpha;


end

