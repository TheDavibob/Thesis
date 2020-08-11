%% FarField_ImpImp_Continuous_Feb2019
% Computes the far-field noise due to a generic junction between a pair of
% lined solutions. Note: this is forked from HSH_Continuous, which is
% slightly snazzier code than earlier.

%% Parameter setup
delta = 1;
%U = Struct_U_Linear(0,1,delta);
U = Struct_U_Parabolic(0,1,delta);
c = Struct_c_Const(5);
omega = 1;
k3 = 0;
N_bl = 1000;

Z = 1-1i;
Zu = R^(-1)*Z; % Can be infinity
Zd = R*Z; % As above
[Lu,Ld] = ImpedanceUpDown(Zu,Zd,omega,U.f(0),U.df(0));
% Generates boundary condition operators, as Lu.l1, Lu.l2 etc.

y2 = 0.5;
kappa1 = omega./U.f(y2);

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

%% Y-dependence - generic for all scattering
% Useful as contains the derivatives
[Y] = Y_Continuous_Scattering_Nov2018_Mod(delta,Saddle,k3,omega,U,c,delta,10);

%% Upstream terms (which, in turn, will fix the scalings
Q0 = 1e-10; % Depends on y2 in a variety of reasonably subtle ways, but scales out

% Useful previous work: easy to replace with Generic, too. Note:
% deformation of critical layer into lower HP to "preserve causality" in
% some bollox sense.
N_bly = N_bl;

% Some scalings for the scattered thing:
Wall_Upstream = Y_GreensGeneric(0,y2,kappa1,k3,omega,U,c,Lu.l1,Lu.l2,delta,N_bly,-1);
Incident_Upstream = Q0*( Ld.l1(kappa1).*Wall_Upstream.dphi + Ld.l2(kappa1).*Wall_Upstream.phi );
% i.e. phid(y2)*Ld(phi)(0)

disp('Upstream terms computed');

%% Dispersion relationships
z0.f = @(s) delta*(1-s);
z0.df = @(s) -delta;

decp = Comp_Dec_Aug18A_skip(Saddle,k3,omega,U,c,delta,N_bl,z0);
Du_p = Lu.l1(Saddle).*decp.dphi + Lu.l2(Saddle).*decp.phi;

decm = Comp_Dec_Aug18A_skip(Saddle,k3,omega,U,c,delta,N_bl,z0);
Dd_m = Ld.l1(Saddle).*decm.dphi + Ld.l2(Saddle).*decm.phi;


%% Kernel and factorisations
z_kernel = Contour_Semicircle(delta,1); % So branch cut in LHP

clear prec
prec.N = 5000;

shift = 0.5*(0.5.*(Saddle(1) - Saddle(end)));
[contour.Cp,contour.dCp] = TanhContour(0,3,-1+shift,1+shift);
[contour.Cm,contour.dCm] = TanhContour(0,3,-1-shift,1-shift);
% Contour strictly not symmetric in this case

branch_cut = pi;

Kernel = @(k1) Kernel_Generic(k1,k3,omega,U,c,delta,Lu,Ld,N_bl,z_kernel);

if U.f(0) ~= 0
    Kernel_ff = FFKernel_Impedance(Zu,Zd,omega,U.f(0),c.f(0),k3,0);
else
    Kernel_ff = FFKernel_Impedance(Zu,Zd,omega,U.f(0),c.f(0),k3,0); % This should be fine.
end

[K_factorised] = Fact_Mult_Cont_Scaled_Finite(Kernel,Kernel_ff,contour,prec,branch_cut);

Kernel_p_0 = K_factorised.p(Saddle);
Kernel_m_0 = K_factorised.m(Saddle);

Km_kappa1 = K_factorised.m(kappa1);

disp('Kernel factorised');
%% Cauchy-like terms
Cauchy_kappa1_0 = 1./(Saddle - kappa1);

%% Weighting functions
A = -Incident_Upstream./(1i*Km_kappa1); 

%% Saddle evaluation
Kp = A*Cauchy_kappa1_0./(Kernel_p_0.*Du_p);
Km = A*Cauchy_kappa1_0.*Kernel_m_0./Dd_m;

%% Far-field sound

FF_generic = struct('phi',[],'dphi',[],'p',[],'v',[]);
physical_variables = fieldnames(FF_generic);

FFp = FF_generic;
FFm = FF_generic;
FF = FF_generic;


for j = 1:numel(physical_variables)
    s = physical_variables{j};
    FFp.(s) = Kp.*sin(Theta).*Y.(s);
    FFm.(s) = Km.*sin(Theta).*Y.(s);
    FF.(s) = FFp.(s);
    FF.(s)(theta < pi/2) = FFm.(s)(theta<pi/2);
end

disp('All done');