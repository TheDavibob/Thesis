function [I] = PVTrapz(f,a,b,p,epsilon)
%PVTRAPZ Computes a definite principal valued integral using the trapezium
%rule
%   Integrates f from a to b, leaving a region of epsilon on either side of
%   each pole p(j). Uses a brute force trapezium rule -- quadrature would
%   probably be better.

% Created 2018-01-24

% doesn't work if epsilon is too small, particularly if N isn't
% sufficiently large

% Constructs n+1 intervals, where n is the number of poles
n=numel(p);

Pp=[a,sort(p)+epsilon]; % start points of each interval
Pm=[sort(p)-epsilon,b]; % end points of each interval

intervals=[Pp(:),Pm(:)]; % each row is a start and end point of an integral, ish.

N=1000;

integral=zeros(1,n+1);
for j=1:n+1
    range=linspace(intervals(j,1),intervals(j,2),N);
    integral(j) = trapz(range,f(range));
end

I=sum(integral);

end

