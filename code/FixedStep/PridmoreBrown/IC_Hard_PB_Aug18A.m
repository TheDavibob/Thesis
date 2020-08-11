function [Phi0]=IC_Hard_PB_Aug18A(k1,omega,U,c)
% Initial conditions for the hard-wall PB problem in the simplest
% formulation


%smol=1e-10;

%c0 = c.f(0);
%dc0 = c.df(0);
u0=U.f(0);
du0=U.df(0);

C0 = 1i*(omega - u0*k1);
dC0 = -1i*du0*k1;

Phi0 = zeros(2,numel(k1));
Phi0(1,:) = -1./(C0.^3) ; % so that p(0)=1
Phi0(2,:) = 3*dC0./(C0.^4);  % so that v(0)=0


% TEST

end
