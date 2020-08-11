%% The far-field noise for the trailing-edge problem
% Only really optimised for hard-wall problem at the moment
% To do: make work for other Z.

%% Parameter setup
c0 = 5;

% Flow field above plate (x2 > 0)
deltap = 1;
%U0 = 0.3;
Uinfp = 1;
% U = Struct_U_Linear(0.3,1,delta);
Up = Struct_U_Parabolic(U0,Uinfp,deltap);
cp = Struct_c_Const(c0);

% Flow field below plate (function of -x2 > 0)
deltam = 1;
Uinfm = 1;
Um = Struct_U_Parabolic(U0,Uinfm,deltam);
cm = Struct_c_Const(c0);

% Other physical parameters
omega = 5; % Frequency
Z = inf; % Plate impedance, can be 0 or inf
k3 = 0; % Spanwise wavenumber
q0 = 1; % Vortex sheet strength (we'll piss around with this)
y2 = 0.1; % Vortex sheet location

% Integration parameters
N_bl = 1000;

% Derived paramters (do not edit)
kappa1 = omega./Up.f(y2);
% Note: assuming y2 is positive. This could be easily relaxed, but shouldn't be

disp('Parameters set up');


%% Far-field bit
N_theta_p = 300;
N_theta_m = 300;
theta1 = linspace(0,pi,N_theta_p);
theta2 = linspace(pi,2*pi,N_theta_m);

theta = [theta1,theta2];

x1_1 = cos(theta1);
x2_1 = sin(theta1);

x1_2 = cos(theta2);
x2_2 = sin(theta2);

[k1_1] = Steepest_Descent_Contour(x1_1,x2_1,omega,Up.f(deltap),cp.f(deltap),k3,1);
[k1_2] = Steepest_Descent_Contour(x1_2,-x2_2,omega,Um.f(deltam),cm.f(deltam),k3,1);

Theta_1 = diag(k1_1.Theta).';
Saddle_1 = k1_1.K0.*cos(Theta_1) - (k1_1.M).*k1_1.k0/(1-k1_1.M^2);

Theta_2 = diag(k1_2.Theta).';
Saddle_2 = k1_2.K0.*cos(Theta_2) - (k1_2.M).*k1_2.k0/(1-k1_2.M^2);
% The above is the steepest descent point as a function of Mach-modified
% angle, Theta.

Theta = [Theta_1,Theta_2];

%% Boundary condition setup
% Upstream, lined boundary (sameish above and below the plate)
Zed = 1i*omega*Z;

U0 = Up.f(0);
sigma1 = Up.df(0);
sigma2 = Um.df(0);

C0 = @(k1) 1i*(omega - U0*k1);
dC1 = @(k1) -1i*sigma1*k1;
dC2 = @(k1) -1i*sigma2*k1;

% Velocity operator
V1.l1 = @(k1) -C0(k1).^2;
V1.l2 = @(k1) -3*C0(k1).*dC1(k1);
V2.l1 = @(k1) -C0(k1).^2;
V2.l2 = @(k1) -3*C0(k1).*dC2(k1);

% Pressure operator (lazy, since continuous)
P.l1 = @(k1) zeros(size(k1));
P.l2 = @(k1) -C0(k1).^3;

if Z == inf
    L1 = V1;
    L2.l1 = @(k1) -V2.l1(k1);
    L2.l2 = @(k1) -V2.l2(k1);
else
    L1.l1 = @(k1) Zed*V1.l1(k1) - C0(k1).*P.l1(k1);
    L1.l2 = @(k1) Zed*V1.l2(k1) - C0(k1).*P.l2(k1);
    L2.l1 = @(k1) -Zed*V2.l1(k1) - C0(k1).*P.l1(k1);
    L2.l2 = @(k1) -Zed*V2.l2(k1) - C0(k1).*P.l2(k1);
end   

disp('Boundary conditions set up');

%% Two alternative factorisations
disp('Finding and factorising Wiener-Hopf kernel');
% Kernel is the same across all problems, predictably, which helps somewhat
z1 = Contour_Semicircle(deltap,1); % So branch cut in LHP, as with zp earlier
z2 = Contour_Semicircle(deltam,1);

clear prec
% Number of integration points for kernel factorisation
prec.N = 10000;

shift = 0.5*(0.5.*(Saddle_1(1) - Saddle_1(end)));
scale = 1; % Vertical stretching

clear contour
%[contour.Cp,contour.dCp] = TanhContour(0,scale*omega,-omega+shift,omega+shift);
%[contour.Cm,contour.dCm] = TanhContour(0,scale*omega,-omega-shift,omega-shift);
[contour.Cp,contour.dCp]=VeitchContour(Saddle_1(1) - shift,-scale*omega,0.5*omega./Up.f(0));
[contour.Cm,contour.dCm]=VeitchContour(Saddle_1(end) + shift,-scale*omega,0.5*omega./Up.f(0));

% Other contours are available, but the symmetric one does the job nicely.

% Branch cut for multiplicative factorisation, requiring logs. For the
% hard-soft factorisation this is normally trivial.
branch_cut = 3*pi/2;

% "Exact" integrated kernel
Kernel = @(k1) Wake_Kernel(k1,k3,omega,Up,cp,deltap,Um,cm,deltam,Z,N_bl,z1,z2);

% Far-field kernel
Kernel_ff = FFWake_Kernel(k3,omega,Up.f(0),cp.f(0),Up.df(0),-Um.df(0),Z,0);

[K_factorised] = Fact_Mult_Cont_Scaled_Finite(Kernel,Kernel_ff,contour,prec,branch_cut);
Kernel_p_1 = K_factorised.p(Saddle_1);
Kernel_m_1 = K_factorised.m(Saddle_1);
Kernel_p_2 = K_factorised.p(Saddle_2);
Kernel_m_2 = K_factorised.m(Saddle_2);

% Some points require specific evaluation, which will depend on their
% respective generalised half-planes

disp('Kernel found and factorised, setting up integrals')
shift_y2 = 0.001;
Km_kappa1 = K_factorised.m(kappa1);
% Wall_Upstream = Y_GreensGeneric(0,y2,kappa1,k3,omega,Up,cp,L1.l1,L1.l2,deltap,N_bl,-1);
[WallPressure] = Distributed_WallPressure(y2,omega./Up.f(y2),shift_y2,omega,k3,Up,cp,deltap,N_bl,L1.l1,L1.l2,-1);
Pressure_Upstream = -q0*WallPressure.*cp.f(y2).^2./Up.df(y2);
%Velocity_Upstream = Q0*( V1.l1(kappa1).*Wall_Upstream.dphi + V1.l2(kappa1).*Wall_Upstream.phi );
% A = -1i*Pressure_Upstream./Km_kappa1;


%% Dispersion relationship setup
disp('Computing dispersion relationships along contours');
DZ1_1 = Dispersion_Generic_Nov2018(Saddle_1,k3,omega,Up,cp,L1.l1,L1.l2,deltap,N_bl,0);
% DZ2_1_temp = Dispersion_Generic_Nov2018(Saddle_1,k3,omega,Um,cm,L2.l1,L2.l2,deltam,N_bl,0);
DZ2_1 = Dispersion_Generic_Nov2018(Saddle_1,k3,omega,Um,cm,L2.l1,L2.l2,deltam,N_bl,0);

Dh1_1 = Dispersion_Generic_Nov2018(Saddle_1,k3,omega,Up,cp,V1.l1,V1.l2,deltap,N_bl,0);
Dp1_1 = Dispersion_Generic_Nov2018(Saddle_1,k3,omega,Up,cp,P.l1,P.l2,deltap,N_bl,0);
Dh2_1 = Dispersion_Generic_Nov2018(Saddle_1,k3,omega,Um,cm,V2.l1,V2.l2,deltam,N_bl,0);
Dp2_1 = Dispersion_Generic_Nov2018(Saddle_1,k3,omega,Um,cm,P.l1,P.l2,deltam,N_bl,0);

Dw_1 = -Dp1_1.*Dh2_1 - Dp2_1.*Dh1_1;

DZ1_2 = Dispersion_Generic_Nov2018(Saddle_2,k3,omega,Up,cp,L1.l1,L1.l2,deltap,N_bl,0);
DZ2_2 = Dispersion_Generic_Nov2018(Saddle_2,k3,omega,Um,cm,L2.l1,L2.l2,deltam,N_bl,0);
% DZ2_2 = -DZ2_2_temp; % Not sure this is robust.

Dh1_2 = Dispersion_Generic_Nov2018(Saddle_2,k3,omega,Up,cp,V1.l1,V1.l2,deltap,N_bl,0);
Dp1_2 = Dispersion_Generic_Nov2018(Saddle_2,k3,omega,Up,cp,P.l1,P.l2,deltap,N_bl,0);
Dh2_2 = Dispersion_Generic_Nov2018(Saddle_2,k3,omega,Um,cm,V2.l1,V2.l2,deltam,N_bl,0);
Dp2_2 = Dispersion_Generic_Nov2018(Saddle_2,k3,omega,Um,cm,P.l1,P.l2,deltam,N_bl,0);

Dw_2 = -Dp1_2.*Dh2_2 - Dp2_2.*Dh1_2;
% This is more intensive process than necessarily required.

%% Leading coefficient
if Z ~= 0
    I0 = Pressure_Upstream;
else
    I0 = Velocity_Upstream;
end

if Z == inf
    A = I0/(1i*Km_kappa1);
    A_upstream = A;
else
    A = I0/(1i*Zed*Km_kappa1);
    A_upstream = I0/(1i*Km_kappa1);
end

%% Functional dependence
Cauchy_1 = 1./(Saddle_1 - kappa1);
Cauchy_2 = 1./(Saddle_2 - kappa1);

% The "steepest descent" contours
K1_p = A_upstream*Cauchy_1./(DZ1_1.*Kernel_p_1);
K2_p = A_upstream*Cauchy_2./(DZ2_2.*Kernel_p_2);
K1_m = A*Cauchy_1.*DZ2_1.*Kernel_m_1./(Dw_1);
K2_m = A*Cauchy_2.*DZ1_2.*Kernel_m_2./(Dw_2);

% Only because it contains the derivatives
[Y_1] = Y_Continuous_Scattering_Nov2018_Mod(deltap,Saddle_1,k3,omega,Up,cp,deltap,10);
[Y_2] = Y_Continuous_Scattering_Nov2018_Mod(deltam,Saddle_2,k3,omega,Um,cm,deltam,10);


%% Far-field stuff (putting it all together)
FF1p.p = K1_p.*sin(Theta_1).*Y_1.p;
FF1p.v = K1_p.*sin(Theta_1).*Y_1.v;
FF1p.phi = K1_p.*sin(Theta_1).*Y_1.phi;
FF1m.p = K1_m.*sin(Theta_1).*Y_1.p;
FF1m.v = K1_m.*sin(Theta_1).*Y_1.v;
FF1m.phi = K1_m.*sin(Theta_1).*Y_1.phi;

FF2p.p = K2_p.*sin(Theta_2).*Y_2.p;
FF2p.v = K2_p.*sin(Theta_2).*Y_2.v;
FF2p.phi = K2_p.*sin(Theta_2).*Y_2.phi;
FF2m.p = K2_m.*sin(Theta_2).*Y_2.p;
FF2m.v = K2_m.*sin(Theta_2).*Y_2.v;
FF2m.phi = K2_m.*sin(Theta_2).*Y_2.phi;

FFp.phi = [FF1p.phi,FF2p.phi];
FFp.p = [FF1p.p,FF2p.p];
FFp.v = [FF1p.v,FF2p.v];
FFm.phi = [FF1m.phi,FF2m.phi];
FFm.p = [FF1m.p,FF2m.p];
FFm.v = [FF1m.v,FF2m.v];

FF = FFp;
FF.p( (theta < pi/2) | (theta > 3*pi/2)) = FFm.p((theta < pi/2) | (theta > 3*pi/2));
FF.v((theta < pi/2) | (theta > 3*pi/2)) = FFm.v((theta < pi/2) | (theta > 3*pi/2));
FF.phi((theta < pi/2) | (theta > 3*pi/2)) = FFm.phi((theta < pi/2) | (theta > 3*pi/2));
% polarplot(theta,abs(FF.p.'));
