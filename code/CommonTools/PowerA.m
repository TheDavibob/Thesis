function [x] = PowerA(x2,power,cut)
%PowerA Complex power with arbitrary cut.
%   Essentially straightforward use of LogA

% Note care is needed with the angle of this cut - -pi/2 and 3pi/2 are not
% the same (they will differ by 2ipi).

x3=LogA(x2,cut);
x=exp(power*x3);

end