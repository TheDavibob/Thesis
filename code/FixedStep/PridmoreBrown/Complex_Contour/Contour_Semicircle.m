function [ z ] = Contour_Semicircle( delta ,HP)
%Contour_Semicircle: an integration path going from delta to 0 via a
% semicircular contour

% if HP > 0, goes via the upper half-plane
% if HP < 0, goes via the lower half-plane

% In all cases, z(0) = delta and z(1) = 0

% z is a structured function

if HP > 0
    I = 1;
elseif HP < 0
    I = -1;
else
    error('Half-plane of integration unclear')
end

z.f = @(s) (delta/2)*(1+ exp(I*1i*pi*s)) ;
z.df = @(s) (delta/2)*I*1i*pi*exp(I*1i*pi*s);

%z.s = @(x) 1/(1i*pi)*LogA( I*(2*x/delta - 1) ,pi);

end

