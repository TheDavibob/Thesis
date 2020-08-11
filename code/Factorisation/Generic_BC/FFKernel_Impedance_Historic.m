function [Kff] = FFKernel_Impedance_Historic(Zu,Zd,omega,U0,c0,k3,W)
%FFKERNEL_IMPEDANCE factorises far-field K for two impedance boundary
%conditions

% THIS IS WRONG WRONG WRONG DO NOT USE: Had slightly buggered up the
% boundary conditions

% NOT EXACTLY CORRECT FOR U0 = 0: see the hard/soft case, and use the other
% code, hopefully.

beta = sqrt(1-U0^2/c0^2);
redU = U0./(beta*omega);
K_inftyplus = (Zd - redU) ./ (Zu - redU) ;
K_inftyminus = (Zd + redU) ./ (Zu + redU) ;


A = K_inftyplus ;

% The tricky exponent
L = 1/(2*pi*1i) * log( K_inftyminus/K_inftyplus );

alpha = L + W;
beta = - L - W;

kbp = sqrt((omega/c0)^2 - k3^2);
kbm = -sqrt((omega/c0)^2 - k3^2);

Kff.p = @(k1) A*PowerA(k1 - kbp,alpha,3*pi/2);
Kff.m = @(k1) PowerA(k1 - kbm,beta,pi/2);

Kff.f = @(k1) Kff.m(k1).*Kff.p(k1);

end

