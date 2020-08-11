function Kff = FFKernel_HardSoft_Slip(omega,U0,c0,k3,W)
% The far field kernel, factorised reasonably trivially, in the slipping
% case.

% W the "winding number" of the normal thing around 0.

% All flow variables evaluated on the wall.

[gamma,gammap,gammam] = Gamma_FF(k3,omega,c0,U0);

C = @(k1) 1i*(omega  - U0*k1);

% Kff.f = @(k1) -gamma(k1)./(c0^2*C(k1));

% Kff.p = @(k1) -gammap(k1)./(c0^2*C(k1));
% Kff.m = @(k1) gammam(k1);

M = U0/c0;
k0 = omega/c0;

kb1 = -M/(1-M^2)*k0 + sqrt(k0^2 + k3^2*(1-M^2))./(1-M^2);
kb2 = -M/(1-M^2)*k0 - sqrt(k0^2 + k3^2*(1-M^2))./(1-M^2);

kp = kb2;
km = kb1;


% If W = 0, then have the trivial factorisation.

Log_Continuity_Factor = @(k1) (k1-kp).^W./(k1-km).^W; % There are probably better choices.
% This satisfies problems with where the poles end up.

Kff.f = @(k1) -(C(k1)./gamma(k1)).*Log_Continuity_Factor(k1);
Kff.p = @(k1) -(C(k1)./gammap(k1))./(k1 - km).^W;
Kff.m = @(k1) (1./gammam(k1)).*(k1 - kp).^W;



end