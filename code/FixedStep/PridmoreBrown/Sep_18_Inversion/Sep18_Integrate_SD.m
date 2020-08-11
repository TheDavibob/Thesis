function [I,Integrand,Hard,Dec,Wron,k1] = Sep18_Integrate_SD(x1,x2,y2,k3,omega,U,c,delta)
% Integrand for the Steepest Descent integration: note that this doesn't
% require particularly careful integration across the boundary-layer, as
% we're remote from the critical-layer.

% k1 will be a big 3-array, so be wary of crippling things.

% UNFINISHED: A HECK OF A LOT OF AWKWARD WORK IS NEEDED

disp('Setting up integration');

% Integration across the boundary-layer
N_SD = 100;

% Big Grid of evaluation points
[k1] = Steepest_Descent_Contour(x1,x2,omega,U.f(delta),c.f(delta),k3,N_SD);

X1 = k1.X1;
X2 = k1.X2;

X1b = repmat(X1,1,1,N_SD); % Barely used until the very end
X2b = repmat(X2,1,1,N_SD); % So both the same size as K1, below
K1 = k1.SD;

Dec.phi = zeros(size(K1));
Dec.p = zeros(size(K1));
Dec.v = zeros(size(K1));
Dec.dphi = zeros(size(K1));
Hard = Dec;

% Integration across the boundary-layer
disp('Computing auxiliary functions within the boundary-layer');
N_BL = 100;
[dec] = Dec_Aug18A_GridX_skip(X2b( (X2b > y2) & (X2b < delta) ),K1( (X2b > y2) & (X2b < delta) ),k3,omega,U,c,delta,N_BL);
[hard] = Hard_Aug18A_GridX_skip(X2b( (X2b < y2) & (X2b < delta) ),K1( (X2b < y2) & (X2b < delta) ),k3,omega,U,c,delta,N_BL);

Dec.phi( (X2b > y2) & (X2b < delta) ) = dec.phi;
Dec.dphi( (X2b > y2) & (X2b < delta) ) = dec.dphi;
Dec.p( (X2b > y2) & (X2b < delta) ) = dec.p;
Dec.v( (X2b > y2) & (X2b < delta) ) = dec.v;

Hard.phi( (X2b < y2) & (X2b < delta) ) = hard.phi;
Hard.dphi( (X2b < y2) & (X2b < delta) ) = hard.dphi;
Hard.p( (X2b < y2) & (X2b < delta) ) = hard.p;
Hard.v( (X2b < y2) & (X2b < delta) ) = hard.v;

% Outwith the boundary-layer: "known" results
disp('Computing auxiliary functions outside the boundary-layer');
gammainf = Gamma_FF(k3,omega,c.f(delta),U.f(delta));
G_inf = gammainf(K1);
C_inf = 1i*(omega - U.f(delta)*K1);
dC_inf = -1i*U.df(delta)*K1;

Dec.phi(X2b >= delta) = exp(-G_inf(X2b>=delta).*X2b(X2b>=delta));
Dec.dphi(X2b >= delta) = -G_inf(X2b>= delta).*exp(-G_inf(X2b>=delta).*X2b(X2b>=delta));

Dec.p(X2b >= delta) = -C_inf(X2b>= delta).^3.*Dec.phi(X2b >= delta);
Dec.v(X2b >= delta) = -3*dC_inf(X2b>= delta).*C_inf(X2b>= delta).*Dec.phi(X2b >= delta)-C_inf(X2b>= delta).^2.*Dec.dphi(X2b >= delta);

% Can also do the same for the hard-wall solution, though if the source is
% within the boundary-layer this is not necessary. FUTURE WORK

% Internal solution: y2 < delta: note that we always use .phi as this isn't
% really differentiated at any point: things end up being discontinuous by
% construction. This allows us to just multiply for PlusMinus.
disp('Computing auxiliary functions at y2');
[dec_y2] = Dec_Aug18A_GridX_skip(y2*ones(size(K1((X2b <= y2) & (X2b < delta)))),K1( (X2b <= y2) & (X2b < delta) ),k3,omega,U,c,delta,N_BL);
[hard_y2] = Hard_Aug18A_GridX_skip(y2*ones(size(K1(X2b >= y2))),K1(X2b >=y2),k3,omega,U,c,delta,N_BL);

Dec.phi( (X2b <= y2) & (X2b < delta) ) = dec_y2.phi;
Dec.dphi( (X2b <= y2) & (X2b < delta) ) = dec_y2.phi;
Dec.p( (X2b <= y2) & (X2b < delta) ) = dec_y2.phi;
Dec.v( (X2b <= y2) & (X2b < delta) ) = dec_y2.phi;

Hard.phi( (X2b >= y2) ) = hard_y2.phi;
Hard.dphi( (X2b >= y2) ) = hard_y2.phi;
Hard.p( (X2b >= y2) ) = hard_y2.phi;
Hard.v( (X2b >= y2) ) = hard_y2.phi;

% Dispersion relation on the wall (hard-wall, though this should be
% reasonably easy to substitute to whatever the heck we want)

disp('Computing dispersion relationship');
%[dec_0] = Dec_Aug18A_GridX_skip(y2*ones(size(K1)),K1,k3,omega,U,c,delta,N_BL);
[dec_0] = Dec_Aug18A_GridX_skip(zeros(size(K1)),K1,k3,omega,U,c,delta,N_BL);
Disp_hard = dec_0.v;
C_0 = 1i*(omega - U.f(0)*K1);
C_y2 = 1i*(omega - U.f(y2)*K1);
Wron = (c.f(0).^2.*C_0.^4)./(c.f(y2).^2.*C_y2.^4).*Disp_hard;
% This ensures the jump in dphi is unity, which is a useful check for
% validity of the solution. This isn't yet correct: need to fudge around
% some things due to the normalisation.

disp('Extra components of the integration');
Integrand.Fourier = exp(-1i*K1.*X1b);
Integrand.Weight = k1.wdt; % Both weights and coordinate transformation contribution
Integrand.Cauchy = 1./(2*pi*1i*(K1- omega./U.f(y2)));
Integrand.Disp_rec = 1./Wron;

Fourier = Integrand.Fourier;
Weight = Integrand.Weight;
Cauchy = Integrand.Cauchy;
Disp_rec = Integrand.Disp_rec;

Integrand.phi = Fourier.*Dec.phi.*Hard.phi.*Weight.*Cauchy.*Disp_rec;
Integrand.dphi = Fourier.*Dec.dphi.*Hard.dphi.*Weight.*Cauchy.*Disp_rec;
Integrand.p = Fourier.*Dec.p.*Hard.p.*Weight.*Cauchy.*Disp_rec;
Integrand.v = Fourier.*Dec.v.*Hard.v.*Weight.*Cauchy.*Disp_rec;


disp('Integrating');
I.phi = sum(Integrand.phi,3);
I.dphi = sum(Integrand.dphi,3);
I.p = sum(Integrand.p,3);
I.v = sum(Integrand.v,3);

end