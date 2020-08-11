function [Kff] = FFWake_Kernel(k3,omega,U0,c0,sigma1,sigma2,Z,W)

% There are various different cases of the wake kernel to be considered.
% This function attempts to factorise all of them. Good luck.

% Winding number variation permisable, as usual.

sigma0 = (sigma1 - sigma2) /2;
% sigma2 = -U2.f(0) if defined as positive

% Branch points for clarity reasons
kbp = sqrt((omega/c0)^2 - k3^2);
kbm = -sqrt((omega/c0)^2 - k3^2);

% Common functions
beta0 = sqrt(1-U0^2/c0^2);

% K_inftyplus: the coefficient of k1 at +infty
% K_inftyminus: the coefficient of k1 at -infty
% k1_dep the power of k1 at either end.

if Z == inf
    if U0 ~= 0
        K_inftyplus = 2*1i*U0/beta0 ;
        K_inftyminus = -2*1i*U0/beta0 ;
        k1_dep = 0;
    else
        K_inftyplus = 2*1i*(sigma0 - omega) ;
        K_inftyminus = 2*1i*(sigma0 + omega) ;
        k1_dep = -1;
    end
elseif Z == 0
    if U0 ~= 0
        K_inftyplus = -2*1i*beta0/U0 ;
        K_inftyminus = 2*1i*beta0/U0 ;
        k1_dep = 0;
    else
        K_inftyplus = -1i*( (1/(sigma1 - omega)) + (1/(-sigma2 - omega)) ) ;
        K_inftyminus = 1i*( (1/(sigma1 + omega)) + (1/(-sigma2 + omega)) ) ;
        k1_dep = 1;
    end
else
    Zed = 1i*omega*Z;
    if U0 ~= 0
        K_inftyplus = -2*1i*beta0*Zed/(U0^3);
        K_inftyminus = 2*1i*beta0*Zed/(U0^3);
        k1_dep = -2;
    else
        K_inftyplus = 2*1i*(sigma0 - omega)/Zed;
        K_inftyminus = -2*1i*(sigma0 + omega)/Zed;
        k1_dep = -1;
    end
end

A = K_inftyplus ;
L = (1/(2*pi*1i) * log( K_inftyminus/K_inftyplus ));
alpha = L + k1_dep + W;
beta = - L - W;

Kff.p = @(k1) A*PowerA(k1 - kbp,alpha,3*pi/2);
Kff.m = @(k1) PowerA(k1 - kbm,beta,pi/2);

Kff.f = @(k1) Kff.m(k1).*Kff.p(k1);