function Kff = FFKernel_HardSoft_NoSlip(omega,sigma,c0,k3,W)
% The far field kernel, factorised as per Singh, for the no-slip case.
% W allows messing with the winding number, if required.

% sigma here is U'(0)

A = 1i*(sigma - omega);

% The tricky exponent
L = 1/(2*pi*1i) * log( (sigma + omega)/(omega - sigma) );

alpha = -1/2 + L + W;
beta = -1/2 - L - W;

kbp = sqrt((omega/c0)^2 - k3^2);
kbm = -sqrt((omega/c0)^2 - k3^2);

Kff.p = @(k1) A*PowerA(k1 - kbp,alpha,3*pi/2);
Kff.m = @(k1) PowerA(k1 - kbm,beta,pi/2);

Kff.f = @(k1) Kff.m(k1).*Kff.p(k1);

end