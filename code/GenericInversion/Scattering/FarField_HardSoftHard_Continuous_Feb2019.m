%% FarField_HardSoftHard_Continuous_Feb2019
% Computes the far-field noise due to the double junction setup as per
% previous HSH scattering code.

%% Parameter setup [Fork of HSH_Cts_Feb2010]
delta = 1;
%U = Struct_U_Linear(0,1,delta);
U = Struct_U_Parabolic(0.3,1,delta);
c = Struct_c_Const(5);
omega = 4;
k3 = 0;
N_bl = 1000;

% Vortex sheet location
y2 = 0.4;
kappa1 = omega./U.f(y2);

% Junction locations
x1_u = -7.5;
x1_d = 0;

l_junction = x1_d - x1_u;

disp('Parameters set up');

%% Far-field bit: setting up theta, locating saddle points
N_theta = 1000;
theta = linspace(0,pi,N_theta);

x1 = cos(theta);
x2 = sin(theta);

[k1] = Steepest_Descent_Contour(x1,x2,omega,U.f(delta),c.f(delta),k3,1);
Theta = diag(k1.Theta).';
Saddle = k1.K0.*cos(Theta) - (k1.M).*k1.k0/(1-k1.M^2);
% The above is the steepest descent point as a function of Mach-modified
% angle, Theta.

%% Boundary conditions

C0 = @(k1) 1i*(omega - U.f(0)*k1);
dC0 = @(k1) -1i*U.df(0)*k1;

% Pressure release
l1_p = @(k1) zeros(size(k1));
l2_p = @(k1) -C0(k1).^3;

% Hard-wall
l1_h = @(k1) -C0(k1).^2;
l2_h = @(k1) -3*C0(k1).*dC0(k1);

%% Zero finding stuff
% An assumption: there's a conjugate pair of pressure-release zeros. The
% maths and therefore code assumes this, so the stable case isn't going to work.

smol = 1e-12; % For differentiation of dispersion relationship.
D_l = @(k1) Dispersion_Generic_Nov2018(k1,k3,omega,U,c,l1_p,l2_p,delta,N_bl,1);
dD_l = @(k1) (D_l(k1 + smol) - D_l(k1 - smol)) ./(2*smol);

zeds(1) = NewtonRaphson(D_l,dD_l,omega/U.f(delta)+0.5i);
zeds(2) = conj(zeds(1));

kappap = zeds(1); % The critical, unstable wavenumber.
if imag(kappap) <= 0
    disp('No pressure-release mode not found, ignoring it');
    kappap = 0;
end

disp('Unstable modes determined, computing upstream behaviour.');

% This code does *not* extend the critical-layer integrals, for reasonably
% obvious reasons. It therefore picks up residual contributions directly.
% Given these are critical to the analysis, this isn't much of a suprise.

dDp_kappap = dD_l(kappap); % Appears directly in the mathematics
dDp_zeds = dD_l(zeds); % Appears in the first scattering problem


%% Upstream terms (which, in turn, will fix the scalings
Q0 = 1e-10; % Depends on y2 in a variety of reasonably subtle ways, but scales out

% Useful previous work: easy to replace with Generic, too. Note:
% deformation of critical layer into lower HP to "preserve causality" in
% some bollox sense.
N_bly = N_bl;

% Some scalings for the scattered thing:
Wall_Upstream = Y_GreensHard(0,y2,kappa1,k3,omega,U,c,delta,N_bly,-1);
Pressure_Upstream = Q0*Wall_Upstream.p; % i.e. phid(y2)*ph(0): easy enough to genericise

%% Y-dependence - generic for all scattering
% Useful as contains the derivatives
[Y] = Y_Continuous_Scattering_Nov2018_Mod(delta,Saddle,k3,omega,U,c,delta,10);

%% Dispersion relationships
disp('Finding and factorising the kernel');
% Contours: z0 along real axis, zp and zm semicircles
z0.f = @(s) delta*(1-s);
z0.df = @(s) -delta;
% Dispersion relationships without deformation
dec0 = Comp_Dec_Aug18A_skip(Saddle,k3,omega,U,c,delta,N_bl,z0);
Dh0 = dec0.v;
Dp0 = dec0.p;

% Specific (other) points of interest (note already have dDp_kappap from
% above). Note: all evaluation with a straight contour.
dec_kappap = Comp_Dec_Aug18A_skip(kappap,k3,omega,U,c,delta,N_bl,z0);
Dh_kappap = dec_kappap.v;

%% Kernel and factorisations
% Kernel is the "same" for both scattering problems, which helps somewhat
z_kernel = Contour_Semicircle(delta,1); % So branch cut in LHP, as with zp earlier

clear prec
prec.N = 5000;
% Number of integration points for kernel factorisation

shift = 0.5*(0.5.*(Saddle(1) - Saddle(end)));
[contour.Cp,contour.dCp] = TanhContour(0,3,-1+shift,1+shift);
[contour.Cm,contour.dCm] = TanhContour(0,3,-1-shift,1-shift);
% Contour strictly not symmetric in this case

% Branch cut for multiplicative factorisation, requiring logs. For the
% hard-soft factorisation this is normally trivial.
branch_cut = pi;

% "Exact" integrated kernel
Kernel = @(k1) Kernel_HardSoft(k1,k3,omega,U,c,delta,N_bl,z_kernel);

% Far-field kernel, with allowance for non-slip.
if U.f(0) ~= 0
    Kernel_ff = FFKernel_HardSoft_Slip(omega,U.f(0),c.f(0),k3,0);
else
    Kernel_ff = FFKernel_HardSoft_NoSlip(omega,U.df(0),c.f(0),k3,0);
end

[K_factorised] = Fact_Mult_Cont_Scaled_Finite(Kernel,Kernel_ff,contour,prec,branch_cut);
Kernel_p_0 = K_factorised.p(Saddle);
Kernel_m_0 = K_factorised.m(Saddle);

% Some points require specific evaluation. All are in the generalised LHP,
% so all use Km.
Km_kappa1 = K_factorised.m(kappa1);
Kp_kappa1 = Kernel(kappa1)./Km_kappa1;
Km_kappap = K_factorised.m(kappap);

Km_zeds = K_factorised.m(zeds);
disp('Kernel found and factorised, setting up integrals')

%% Cauchy-like terms
Cauchy_kappa1_0 = 1./(Saddle - kappa1);
Cauchy_kappap_0 = 1./(Saddle - kappap);

Cauchy_kappap_kappa1 = 1./(kappap - kappa1);

%% Weighting functions
A_s1 = -Pressure_Upstream./(1i*Km_kappa1); 
if kappap ~= 0
    A_s2_1 = ( Pressure_Upstream.*Dh_kappap.*Km_kappap.^2.*exp(-1i*kappap*l_junction) ) ...
    ./ (1i*Km_kappa1.*(kappap - kappa1).*dDp_kappap);
else
    A_s2_1 = 0;
end
A_s2_2 = ( Pressure_Upstream*exp(-1i*kappa1*l_junction) )./(1i*Kp_kappa1);

%% Saddle evaluation
% Upstream junction
Kp_s1 = A_s1*Cauchy_kappa1_0./(Kernel_p_0.*Dh0);
Km_s1 = A_s1*Cauchy_kappa1_0.*Kernel_m_0./Dp0;

% Downstream junction: modal unscattering
Kp_s2_1 = A_s2_1*Cauchy_kappap_0.*Kernel_p_0./(Dp0);
Km_s2_1 = A_s2_1*Cauchy_kappap_0./(Kernel_m_0.*Dh0);

% Downstream junction: forcing term
Kp_s2_2 = A_s2_2*Cauchy_kappa1_0.*Kernel_p_0./(Dp0);
Km_s2_2 = A_s2_2*Cauchy_kappa1_0./(Kernel_m_0.*Dh0);

Kp = Kp_s1 + Kp_s2_1 + Kp_s2_2; % Saddle point evaluation of integrand, save Fourier and decaying bit.
Km = Km_s1 + Km_s2_1 + Km_s2_2;

%% Far-field sound

FF_generic = struct('phi',[],'dphi',[],'p',[],'v',[]);
physical_variables = fieldnames(FF_generic);

FFp = FF_generic;
FFm = FF_generic;
FF = FF_generic;

FFp_s1 = FF_generic;
FFm_s1 = FF_generic;
FF_s1 = FF_generic;

FFp_s2_1 = FF_generic;
FFm_s2_1 = FF_generic;
FF_s2_1 = FF_generic;

FFp_s2_2 = FF_generic;
FFm_s2_2 = FF_generic;
FF_s2_2 = FF_generic;

for j = 1:numel(physical_variables)
    s = physical_variables{j};
    FFp.(s) = Kp.*sin(Theta).*Y.(s);
    FFm.(s) = Km.*sin(Theta).*Y.(s);
    FF.(s) = FFp.(s);
    FF.(s)(theta < pi/2) = FFm.(s)(theta<pi/2);
    
    FFp_s1.(s) = Kp_s1.*sin(Theta).*Y.(s);
    FFm_s1.(s) = Km_s1.*sin(Theta).*Y.(s);
    FF_s1.(s) = FFp_s1.(s);
    FF_s1.(s)(theta < pi/2) = FFm_s1.(s)(theta<pi/2);
    
    FFp_s2_1.(s) = Kp_s2_1.*sin(Theta).*Y.(s);
    FFm_s2_1.(s) = Km_s2_1.*sin(Theta).*Y.(s);
    FF_s2_1.(s) = FFp_s2_1.(s);
    FF_s2_1.(s)(theta < pi/2) = FFm_s2_1.(s)(theta<pi/2);
    
    FFp_s2_2.(s) = Kp_s2_2.*sin(Theta).*Y.(s);
    FFm_s2_2.(s) = Km_s2_2.*sin(Theta).*Y.(s);
    FF_s2_2.(s) = FFp_s2_2.(s);
    FF_s2_2.(s)(theta < pi/2) = FFm_s2_2.(s)(theta<pi/2);
end