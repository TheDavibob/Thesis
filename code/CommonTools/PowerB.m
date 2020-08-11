function [x] = PowerB(x2,power,cut)
%PowerB Complex power with arbitrary cut, with power a function of x2.
%   Essentially straightforward use of LogA

% Note care is needed with the angle of this cut - -pi/2 and 3pi/2 are not
% the same (they will differ by 2ipi).

x3=LogA(x2,cut);
x=exp(power(x2).*x3);

end