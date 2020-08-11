function [ZP]=ZerosPoles(fun,dfun,x,y)

% computes Z-P for function fun in the box x1<Re(z)<x2, y1<Im(z)<y2;

% That is, zeros-poles (with multiplicity)
% Analytic input of fun and its derivative (currently)

f=@(z) dfun(z)./(2*pi*1i*fun(z));

ZP= integral(f,x(1)+1i*y(1),x(1)+1i*y(1),'Waypoints',[x(2)+1i*y(1),x(2)+1i*y(2),x(1)+1i*y(2)]);

ZP=real(ZP); % ZP is pretty much always real unless something's gone wrong.
%Also mostly integral, though branch points etc can make that skewy.

end