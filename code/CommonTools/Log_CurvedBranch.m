function [Log] = Log_CurvedBranch(z,f)
% Computes log, with jump across branch cut defined by theta = f(|z|).

% f an anonymous function which takes values from 0 to infinity and from
% this spits out values of theta, which can be anything.

r = abs(z);
theta = pi+f(r)+angle(z.*exp(-1i*(f(r)+pi))); % Between f(r) and f(r) + 2pi

Log = log(r) + 1i*theta;

end

