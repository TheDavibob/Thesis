%% Scattering from a hard plate with distributed vorticity distribution
% A direct fork of lined trailing-edge, which is optimised only for
% hard-wall stuff at time of creation (March 2019) anyway.

%% Parameter setup
% Flow field above plate (x2 > 0)
deltap = 1;
U0 = 0.3;
Uinfp = 1;
% U = Struct_U_Linear(0.3,1,delta);
Up = Struct_U_Parabolic(U0,Uinfp,deltap);
cp = Struct_c_Const(5);

% Flow field below plate (function of -x2 > 0)
deltam = 1;
Uinfm = 1.5;
Um = Struct_U_Parabolic(U0,Uinfm,deltam);
cm = Struct_c_Const(5);

% Other physical parameters
omega = 1; % Frequency
Z = inf; % Plate impedance, can be 0 or inf
k3 = 0; % Spanwise wavenumber
q0 = 1; % Vortex sheet strength

% Grid
x1 = -10:0.05:10;
x2 = -2:0.02:2;
% x2 = 0:0.005:2;

% Integration parameters
N_bl = 2000;
N_distribution = 80;

% Vorticity parameters
[kappa1,dkappa1,y2] = kappa1_Quad(omega,Up,deltap,N_distribution);
shift_y2 = 0.01;
omega_c = ones(size(y2));
y1 = -10;
%omega_c = exp(+1i*y1*kappa1);
%omega_c = y2.*(deltap-y2).*exp(+1i*y1*kappa1);

%omega_c = zeros(size(y2));
%omega_c(50) = 100;

disp('Parameters set up');

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
    L2 = V2;
else
    L1.l1 = @(k1) Zed*V1.l1(k1) - C0(k1).*P.l1(k1);
    L1.l2 = @(k1) Zed*V1.l2(k1) - C0(k1).*P.l2(k1);
    L2.l1 = @(k1) Zed*V2.l1(k1) - C0(k1).*P.l1(k1);
    L2.l2 = @(k1) Zed*V2.l2(k1) - C0(k1).*P.l2(k1);
end   

disp('Boundary conditions set up');


%% SD contour setup, as usual
% Different contours for x2 >< 0, which is a development

% Number of quadrature points
N_C = 500;

% x2 > 0
Minfp = Up.f(deltap)./cp.f(deltap);
k0p = omega./cp.f(deltap);
kb1p = -Minfp/(1-Minfp^2)*k0p + sqrt(k0p^2 + k3^2*(1-Minfp^2))./(1-Minfp^2);
kb2p = -Minfp/(1-Minfp^2)*k0p - sqrt(k0p^2 + k3^2*(1-Minfp^2))./(1-Minfp^2);

% x2 < 0
Minfm = Um.f(deltam)./cm.f(deltam);
k0m = omega./cm.f(deltam);
kb1m = -Minfm/(1-Minfm^2)*k0m + sqrt(k0m^2 + k3^2*(1-Minfm^2))./(1-Minfm^2);
kb2m = -Minfm/(1-Minfm^2)*k0m - sqrt(k0m^2 + k3^2*(1-Minfm^2))./(1-Minfm^2);

% Various adjustment paramters
sep = 0.1;
angle = 3*pi/8;
curvature = 0.01;
apex = 0.1*omega;

% Contour: for x2 > 0
Cp1 = Finite_Symmetric_Hyperbola(N_C,kb2p,-apex,curvature,angle);
Cm1 = Finite_Symmetric_Hyperbola(N_C,kb1p,apex,curvature,angle);

% Contour: for x2 < 0
Cp2 = Finite_Symmetric_Hyperbola(N_C,kb2m,-apex,curvature,angle);
Cm2 = Finite_Symmetric_Hyperbola(N_C,kb1m,apex,curvature,angle);


%% Critical layer contour setup
% Whilst CL contour is the same independently of the sign of x2, it is
% driven by the "longest" critical-layer cut

N_CL = 1000; % High precision required
a = min(omega./Up.f(deltap),omega./Up.f(deltap));
if Up.f(0) == 0
    b = 100*a;
else
    b = omega./Up.f(0);
end
shift1 = 0.1*omega; % Displacement of contour from either end of critical layer
shift2 = 0.1*omega;

[C_CLp,C_CLm] = Finite_CL_Ellipse(N_CL,a,b,shift1,shift2);

disp('Contours set up');

%% Upstream behaviour
disp('Computing upstream behaviour');
% Agressive changes from the previous
[G_01_A] = Upstream_Distributed_A(y2,kappa1,shift_y2,omega,k3,Up,cp,deltap,10*N_bl,L1.l1,L1.l2,1);
[G_01] = Upstream_Distributed_B(x2(x2>=0),G_01_A,omega,k3,Up,cp,deltap);
[phi_0] = Upstream_Distributed_Integrated(x1,omega_c,G_01,dkappa1);

I_upstream = phi_0;

%% Upstream wall pressure

% Incident pressure at the wall
[I0] = Distributed_WallPressure(y2,kappa1,shift_y2,omega,k3,Up,cp,deltap,N_bl,L1.l1,L1.l2,1);

%% x1-dependence
disp('Computed upstream behaviour. Computing physical-space behaviour')
Xp1 = Fourier_Nov2018(x1(x1<0),Cp1.k);
Xm1 = Fourier_Nov2018(x1(x1>=0),Cm1.k);

Xp2 = Fourier_Nov2018(x1(x1<0),Cp2.k);
Xm2 = Fourier_Nov2018(x1(x1>=0),Cm2.k);

X_CLp = Fourier_Nov2018(x1(x1>=0),C_CLp.k);
X_CLm = Fourier_Nov2018(x1(x1>=0),C_CLm.k);

%% x2-dependence, including derivatives, a rate limiting step
[Yp1] = Y_Continuous_Scattering_Nov2018_Mod(x2(x2>=0),Cp1.k,k3,omega,Up,cp,deltap,N_bl);
[Ym1] = Y_Continuous_Scattering_Nov2018_Mod(x2(x2>=0),Cm1.k,k3,omega,Up,cp,deltap,N_bl);

[Yp2] = Y_Continuous_Scattering_Nov2018_Mod(sort(-x2(x2<0)),Cp2.k,k3,omega,Um,cm,deltam,N_bl);
[Ym2] = Y_Continuous_Scattering_Nov2018_Mod(sort(-x2(x2<0)),Cm2.k,k3,omega,Um,cm,deltam,N_bl);

[Y_CLp1] = Y_Continuous_Scattering_CL_Nov2018_Mod(x2(x2>=0),C_CLp.k,k3,omega,Up,cp,deltap,N_bl,1);
[Y_CLm1] = Y_Continuous_Scattering_CL_Nov2018_Mod(x2(x2>=0),C_CLm.k,k3,omega,Up,cp,deltap,N_bl,-1);

[Y_CLp2] = Y_Continuous_Scattering_CL_Nov2018_Mod(sort(-x2(x2<0)),C_CLp.k,k3,omega,Um,cm,deltam,N_bl,1);
[Y_CLm2] = Y_Continuous_Scattering_CL_Nov2018_Mod(sort(-x2(x2<0)),C_CLm.k,k3,omega,Um,cm,deltam,N_bl,-1);

%% Dispersion relationship setup
disp('Computing dispersion relationships along contours');
DZ1 = @(k1,HP) Dispersion_Generic_Nov2018(k1,k3,omega,Up,cp,L1.l1,L1.l2,deltap,N_bl,HP);
DZ2_temp = @(k1,HP) Dispersion_Generic_Nov2018(k1,k3,omega,Um,cm,L2.l1,L2.l2,deltam,N_bl,HP);
DZ2 = @(k1,HP) -DZ2_temp(k1,HP); % Not sure this is robust.

Dh1 = @(k1,HP) Dispersion_Generic_Nov2018(k1,k3,omega,Up,cp,V1.l1,V1.l2,deltap,N_bl,HP);
Dp1 = @(k1,HP) Dispersion_Generic_Nov2018(k1,k3,omega,Up,cp,P.l1,P.l2,deltap,N_bl,HP);
Dh2 = @(k1,HP) Dispersion_Generic_Nov2018(k1,k3,omega,Um,cm,V2.l1,V2.l2,deltam,N_bl,HP);
Dp2 = @(k1,HP) Dispersion_Generic_Nov2018(k1,k3,omega,Um,cm,P.l1,P.l2,deltam,N_bl,HP);

Dw = @(k1,HP) -Dp1(k1,HP).*Dh2(k1,HP) - Dp2(k1,HP).*Dh1(k1,HP);
% This is more intensive process than necessarily required.

%% Mode location determination
disp('Locating modes')

smol = 1e-12;
dDw = @(k1) (Dw(k1+smol,0) - Dw(k1-smol,0))./(2*smol);

% Assuming, for now, the only modes are in the gLHP and thus only the wake
% ones are important. We guess where they might be % omega/Up.f(deltap) + 
[zeros_wake_m_guess] = NewtonRaphson(@(k1) Dw(k1,0),dDw,[0.9+ 0.35i]);
% For standard asymmetric, omega = 1, 0.9+0.35i

% zeros_wake_m_guess = 10.7202395 + 0.2366665i;
% For standard asymmetric, omega = 5, 11 + 0.2i
zeros_wake_m = [zeros_wake_m_guess,conj(zeros_wake_m_guess)];

dDw_wake_m = dDw(zeros_wake_m);

disp('Located modes');
%% Dispersion relationships along specific contours

% Along SD contours
DZ1p = DZ1(Cp1.k,0);
DZ2p = DZ2(Cp2.k,0);

DZ1m = DZ1(Cm1.k,0);
DZ2m = DZ2(Cm2.k,0);

Dw1m = Dw(Cm1.k,0);
Dw2m = Dw(Cm2.k,0);

% Along CL contours
DZ1_CLp = DZ1(C_CLp.k,1);
DZ1_CLm = DZ1(C_CLm.k,-1);
DZ2_CLp = DZ2(C_CLp.k,1);
DZ2_CLm = DZ2(C_CLm.k,-1);
Dw_CLp = Dw(C_CLp.k,1);
Dw_CLm = Dw(C_CLm.k,-1);

disp('Computed dispersion relationships along integration contours')

%% Kernel factorisation
disp('Finding and factorising Wiener-Hopf kernel');
% Kernel is the same across all problems, predictably, which helps somewhat
z1 = Contour_Semicircle(deltap,1); % So branch cut in LHP, as with zp earlier
z2 = Contour_Semicircle(deltam,1);

clear prec
prec.N = 2000;
% Number of integration points for kernel factorisation

[contour.Cp,contour.dCp] = TanhContour(0,1,-1,1);
[contour.Cm,contour.dCm] = TanhContour(0,1,-1,1);
% Other contours are available, but the symmetric one does the job nicely.

% Branch cut for multiplicative factorisation, requiring logs. For the
% hard-soft factorisation this is normally trivial.
branch_cut = pi;

% "Exact" integrated kernel
Kernel = @(k1) Wake_Kernel(k1,k3,omega,Up,cp,deltap,Um,cm,deltam,Z,N_bl,z1,z2);

% Far-field kernel
Kernel_ff = FFWake_Kernel(k3,omega,Up.f(0),cp.f(0),Up.df(0),-Um.df(0),Z,0);

[K_factorised] = Fact_Mult_Cont_Scaled_Finite(Kernel,Kernel_ff,contour,prec,branch_cut);
Kernel_p_p1 = K_factorised.p(Cp1.k);
Kernel_m_m1 = K_factorised.m(Cm1.k);
Kernel_p_p2 = K_factorised.p(Cp2.k);
Kernel_m_m2 = K_factorised.m(Cm2.k);

Kernel_m_CLp = K_factorised.m(C_CLp.k);
Kernel_m_CLm = K_factorised.m(C_CLm.k);

% Some points require specific evaluation, which will depend on their
% respective generalised half-planes

Km_kappa1 = K_factorised.m(kappa1);

disp('Kernel found and factorised, setting up integrals')

%% Computing upstream pressure-term
Km_distributed = K_factorised.m(kappa1);
Int = -omega*cp.f(y2).^2.*I0.*omega_c./(Km_distributed.*Up.f(y2).^2);
F_bar = @(k1) Fbar(k1,kappa1,dkappa1,Int);

%% Cauchy-like term
Cauchy_Cp1 = F_bar(Cp1.k);
Cauchy_Cp2 = F_bar(Cp2.k);
Cauchy_Cm1 = F_bar(Cm1.k);
Cauchy_Cm2 = F_bar(Cm2.k);
Cauchy_CLp = F_bar(C_CLp.k);
Cauchy_CLm = F_bar(C_CLm.k);
disp('Computed forcing from upstream terms');

%% Leading coefficient
A = 1./(1i); % Much simpler than before

%% k1-dependence

% The "steepest descent" contours
K_Cp1 = A*Cauchy_Cp1./(DZ1p.*Kernel_p_p1);
K_Cp2 = A*Cauchy_Cp2./(DZ2p.*Kernel_p_p2);
K_Cm1 = A*Cauchy_Cm1.*DZ2m.*Kernel_m_m1./(Dw1m);
K_Cm2 = A*Cauchy_Cm2.*DZ1m.*Kernel_m_m2./(Dw2m);

% The "critical-layer" contours
K_CLp1 = A*Cauchy_CLp.*DZ2_CLp.*Kernel_m_CLp./(Dw_CLp);
K_CLm1 = A*Cauchy_CLm.*DZ2_CLm.*Kernel_m_CLm./(Dw_CLm);
K_CLp2 = A*Cauchy_CLp.*DZ1_CLp.*Kernel_m_CLp./(Dw_CLp);
K_CLm2 = A*Cauchy_CLm.*DZ1_CLm.*Kernel_m_CLm./(Dw_CLm);

disp('Constructed integrands');

%% Integration: acoustic

disp('Computing steepest descent integral');
I_generic = struct('phi',[],'dphi',[],'p',[],'v',[]);
physical_variables = fieldnames(I_generic);
physical_sign = [1,-1,1,-1];
% For x2 < 0, some of the signs need to be flipped. This does it.
% If any fields are deleted, delete their sign, too.

I_acoustic_Cp1 = I_generic;
I_acoustic_Cp2 = I_generic;
I_acoustic_Cm1 = I_generic;
I_acoustic_Cm2 = I_generic;
I_acoustic_1 = I_generic;
I_acoustic_2 = I_generic;
I_acoustic = I_generic;
for j = 1:numel(physical_variables)
    s = physical_variables{j};
    I_acoustic_Cp1.(s) = Integral_Nov2018(Cp1.wk,K_Cp1,Xp1,Yp1.(s));
    I_acoustic_Cp2.(s) = Integral_Nov2018(Cp2.wk,K_Cp2,Xp2,Yp2.(s));
    I_acoustic_Cm1.(s) = Integral_Nov2018(Cm1.wk,K_Cm1,Xm1,Ym1.(s));
    I_acoustic_Cm2.(s) = Integral_Nov2018(Cm2.wk,K_Cm2,Xm2,Ym2.(s));
    I_acoustic_1.(s) = [I_acoustic_Cp1.(s);I_acoustic_Cm1.(s)];
    I_acoustic_2.(s) = physical_sign(j)*[I_acoustic_Cp2.(s)(:,end:-1:1);I_acoustic_Cm2.(s)(:,end:-1:1)];
    I_acoustic.(s) = [I_acoustic_2.(s),I_acoustic_1.(s)];
end

%% Integration: critical-layer
disp('Computing critical-layer integrals');

I_CLp1 = I_generic;
I_CLp2 = I_generic;
I_CLm1 = I_generic;
I_CLm2 = I_generic;
I_CL1 = I_generic;
I_CL2 = I_generic;
I_CL = I_generic;

for j = 1:numel(physical_variables)
    s = physical_variables{j};
    I_CLp1.(s) = Integral_Nov2018(C_CLp.wk,K_CLp1,X_CLp,Y_CLp1.(s));
    I_CLp2.(s) = Integral_Nov2018(C_CLp.wk,K_CLp2,X_CLp,Y_CLp2.(s));
    I_CLm1.(s) = Integral_Nov2018(C_CLm.wk,K_CLm1,X_CLm,Y_CLm1.(s));
    I_CLm2.(s) = Integral_Nov2018(C_CLm.wk,K_CLm2,X_CLm,Y_CLm2.(s));
    I_CL1.(s) = I_CLp1.(s)+I_CLm1.(s);
    I_CL2.(s) = physical_sign(j)*I_CLp2.(s)(:,end:-1:1)+physical_sign(j)*I_CLm2.(s)(:,end:-1:1);
    I_CL.(s) = [I_CL2.(s),I_CL1.(s)];
end


%% Integration: modal contributions
disp('Determining modal contributions')


% Downstream wake modes
N_wake_m = numel(zeros_wake_m);

[Y1_wake_m] = Y_Continuous_Scattering_Nov2018_Mod(x2(x2>=0),zeros_wake_m,k3,omega,Up,cp,deltap,N_bl);
[Y2_wake_m] = Y_Continuous_Scattering_Nov2018_Mod(sort(-x2(x2<0)),zeros_wake_m,k3,omega,Um,cm,deltam,N_bl);

X_wake_m = Fourier_Nov2018(x1(x1>=0),zeros_wake_m);

K1_wake_m = -2*pi*1i*A*DZ2(zeros_wake_m,0).*K_factorised.m(zeros_wake_m).*F_bar(zeros_wake_m)./(dDw_wake_m);
K2_wake_m = -2*pi*1i*A*DZ1(zeros_wake_m,0).*K_factorised.m(zeros_wake_m).*F_bar(zeros_wake_m)./(dDw_wake_m);

I_wake_m1 = cell(N_wake_m,1);
I_wake_m2 = cell(N_wake_m,1);
I_wake_m = cell(N_wake_m,1);

I_wake_m_total = I_generic;
for j = 1:numel(physical_variables)
    s = physical_variables{j};
    I_wake_m_total.(s) = zeros(numel(x1(x1>=0)),numel(x2));
end

for k = 1:N_wake_m
    I_wake_m1{k} = I_generic;
    I_wake_m2{k} = I_generic;
    I_wake_m{k} = I_generic;
    for j = 1:numel(physical_variables)
        s = physical_variables{j};
        I_wake_m1{k}.(s) = K1_wake_m(k)*X_wake_m(:,k)*Y1_wake_m.(s)(:,k).';
        I_wake_m2{k}.(s) = physical_sign(j)*K2_wake_m(k)*X_wake_m(:,k)*Y2_wake_m.(s)(:,k).';
        I_wake_m{k}.(s) = [I_wake_m2{k}.(s)(:,end:-1:1),I_wake_m1{k}.(s)];
        I_wake_m_total.(s) = I_wake_m_total.(s) + I_wake_m{k}.(s);
    end
end


%% Combinations

Z11 = zeros(numel(x1(x1<0)),numel(x2(x2>=0)));
Z12 = zeros(numel(x1(x1<0)),numel(x2(x2<0)));
Z21 = zeros(numel(x1(x1>=0)),numel(x2(x2>=0)));
Z22 = zeros(numel(x1(x1>=0)),numel(x2(x2<0)));

% Firstly: extending various things to cover all x1, x2
I_modal = I_generic;
I_crit = I_generic;
I_upstream_extended = I_generic;

for j = 1:numel(physical_variables)
    s = physical_variables{j};
    I_modal.(s) = [ [Z12,Z11];I_wake_m_total.(s)] ;
    I_crit.(s) = [ [Z12,Z11];I_CL.(s)];
    I_upstream_extended.(s) = [ [Z12;Z22],I_upstream.(s)];
end

% Secondly: adding things together in a variety of ways
I_scattered = I_generic;
I_total = I_generic;

for j = 1:numel(physical_variables)
    s = physical_variables{j};
    I_scattered.(s) = I_modal.(s) + I_crit.(s) + I_acoustic.(s) ;
    I_total.(s) = I_scattered.(s) + I_upstream_extended.(s);
end