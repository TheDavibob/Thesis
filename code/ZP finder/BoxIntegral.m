function [Int]=BoxIntegral(fun,x,y)

% Integrates function fun around a closed box -- should indicate the
% presence (or not) of poles.

f=@(z) fun(z)./(2*pi*1i); % obtain the sum of the residues

Int= integral(f,x(1)+1i*y(1),x(1)+1i*y(1),'Waypoints',[x(2)+1i*y(1),x(2)+1i*y(2),x(1)+1i*y(2)]);

end