function [C] = CL_Ellipse_Contour(kinf,k0,a,b)
%CL_ELLIPSE_CONTOUR Elliptical contour looping the critical layer

% Critical layer from kinf to k0, with the contour passing -a from either
% end in the real direction and b in the imaginary direction.

if kinf > k0
    error('End of critical layer must be larger than start')
end

dk = 0.5*(k0 - kinf);
sk = 0.5*(k0 + kinf);
scale = 1./sqrt(1 - dk^2/(dk+a)^2);

% Real part of the function, t running from -1 to 1
x.f = @(t) ( a +dk )*t + sk;
x.df = @(t) ( a + dk )*ones(size(t));

% In the upper half plane
C.U = @(t) x.f(t) + 1i*b*scale*sqrt(1 - ( x.f(t)-sk).^2/((dk+a).^2) ) ;
C.dU = @(t) x.df(t) + 1i*b*scale*x.df(t).*( - (x.f(t)-sk)/(dk+a)^2./sqrt(1 - ( x.f(t)-sk).^2/((dk+a).^2) ) ) ;

% In the lower half-plane
C.L = @(t)  x.f(t) - 1i*b*scale*sqrt(1 - ( x.f(t)-sk).^2/((dk+a).^2) ) ;
C.dL = @(t) x.df(t) - 1i*b*scale*x.df(t).*( - (x.f(t)-sk)/(dk+a)^2./sqrt(1 - ( x.f(t)-sk).^2/((dk+a).^2) ) ) ;

end

