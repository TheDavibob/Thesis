function [sumzj]=SumZerosj(fun,dfun,x,y,j)

% Computes the sum (z_n^j) of the zeros of fun. This allows finding the
% zeros.

f=@(z) z.^j.*dfun(z)./(2*pi*1i*fun(z));

sumzj= integral(f,x(1)+1i*y(1),x(1)+1i*y(1),'Waypoints',[x(2)+1i*y(1),x(2)+1i*y(2),x(1)+1i*y(2)]);

end
