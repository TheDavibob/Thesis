function I = ComplexIntegral_Line(fun,z1,z2,N)
% LinearIntegral - Integrates fun from z1 to z2, where zj can be complex.

% Uses standard Gaussian quadrature methods to compute int fun dz from z1
% to z2. The quadrature precision is given by N.

[s,w] = GaussLegendre(N) ; % Points and weights for the quadrature

z = (z1+z2)/2 + s*(z2-z1)/2 ; % Contour paramterisation

I2 = fun(z).*w.*(z2-z1)/2;
I = sum(I2);

end