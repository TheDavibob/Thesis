function [phiinf]=IC_Decaying_PB_Aug18A(k1,k3,omega,U,c,delta)
% Initial conditions for the decaying PB problem, evaluated at x2 = delta

Cinf = 1i*(omega - U.f(delta)*k1);
dCinf =-1i*U.df(delta)*k1;
Uinf = U.f(delta);
cinf=c.f(delta);

[gamma,~,~]=Gamma_FF(k3,omega,cinf,Uinf) ;

% pinf = -Cinf.^3*exp(-delta*gamma(k1)); % Pressure at edge of BL
% vinf = -Cinf.^2*exp(-delta*gamma(k1)); % Velocity at edge of BL


phiinf = zeros(2,numel(k1));
phiinf(1,:) = exp(-delta*gamma(k1)) ; % phi is continuous
phiinf(2,:) = ( -gamma(k1)-3*dCinf./Cinf ).*exp(-delta*gamma(k1)); % so v is continuous, not quite phi'
% so p and v are continuous -- jump in shear adds a little.

% Might be quite a lot neater if simply have phi ~ exp(-gamma x2) outside
% the boundary layer.
% phiinf(1,:) = exp(-delta*gamma(k1)) ;
% phiinf(2,:) = -gamma(k1).*(1i*cinf^2*(omega-Uinf*k1)).*phiinf(1,:);
% unchanged other than scaling due to phi


end