function [gamma,gammap,gammam]=Gamma_FF(k3,omega,cinf,Uinf)
% Computes gamma, the coefficient in the exponent for large x2 beyond a
% boundary-layer, as a function of k1.

% Also computes +- factorisation, as usual.

% gamma^2 = (k1^2+k3^2) - (omega - Uk1)^2/c^2, branch cuts chosen to avoid
% the real axis when Im(omega)<0. Not scaled with respect to c, unlike some
% other code I've previously written.

% Works for all k3, k3 a single parameter.

% Created 2018-02-01
% Fixed 2018-04-27

% Error checks
if isreal(Uinf) ~= 1
    error('Uinf must be real')
elseif isreal(cinf) ~=1
    error('cinf must be real')
elseif abs(Uinf)>cinf
    error('Flow must be subsonic') ;
end

% The scaling can be done earlier by just setting cinf=1
M = Uinf/cinf;
k0 = omega/cinf; 

% Cp is the square root of C2, either in the UHP or on the negative real
% axis. Cm is the other one.
C2 = k0^2/((1-M^2)^2) - k3^2/(1-M^2) ;
Cpm = [sqrt(C2),-sqrt(C2)];
if imag(Cpm(1)) == 0
    Cp = min(Cpm);
    Cm = max(Cpm);
elseif imag(Cpm(1)) > 0
    Cp = Cpm(1);
    Cm = Cpm(2);
else
    Cp=Cpm(2);
    Cm=Cpm(1);
end

beta = sqrt(1-M^2);

gammap = @(k1) sqrt(beta)*SqrtA(  ( k1 + M*k0/(1-M^2) ) - Cm , 3*pi/2) ; % analytic in the UHP, branch point at Cm
gammam = @(k1) sqrt(beta)*SqrtA(  ( k1 + M*k0/(1-M^2) ) - Cp , pi/2) ; % analytic in the LHP, branch point at Cp

gamma = @(k1) gammap(k1).*gammam(k1);

end