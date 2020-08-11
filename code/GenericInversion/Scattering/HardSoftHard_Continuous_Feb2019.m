%% Double scattering problem: hard - soft - hard, with approximations

%% Parameter setup
delta = 1;
%U = Struct_U_Linear(0,1,delta);
U = Struct_U_Parabolic(0.3,1,delta);
c = Struct_c_Const(5);
omega = 1;
k3 = 0;
N_bl = 1000;

% Vortex sheet location
y2 = 0.4;
kappa1 = omega./U.f(y2);

% Grid
x1 = -15:0.1:15;
x2 = 0:0.02:4;

% Junction locations
x1_u = -5;
x1_d = 0;

l_junction = x1_d - x1_u;

disp('Parameters set up');

%% Contour setup
N_C = 1000;

Minf = U.f(delta)./c.f(delta);
k0 = omega./c.f(delta);
kb1 = -Minf/(1-Minf^2)*k0 + sqrt(k0^2 + k3^2*(1-Minf^2))./(1-Minf^2);
kb2 = -Minf/(1-Minf^2)*k0 - sqrt(k0^2 + k3^2*(1-Minf^2))./(1-Minf^2);

sep = 0.1;

angle = pi/4;
curvature = 0.01;
apex = 0.1*omega;
Cp = Finite_Symmetric_Hyperbola(N_C,kb2,-apex,curvature,angle);
Cm = Finite_Symmetric_Hyperbola(N_C,kb1,apex,curvature,angle);

%% Critical layer contours

N_CL = 1000;
a = omega./U.f(delta);
if U.f(0)~= 0
    b = omega./U.f(0);
else
    b = 100*a; % A fudge, but useful.
end
shift1 = 0.1*omega;
shift2 = 0.1*omega;

[C_CLp,C_CLm] = Finite_CL_Ellipse(N_CL,a,b,shift1,shift2);
disp('Contours set up');

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
    error('Unstable pressure-release mode not found');
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
N_bly = 10*N_bl;

Y_upstream = Y_GreensHard(x2,y2,kappa1,k3,omega,U,c,delta,N_bly,-1);
X_upstream = exp(-1i*x1.'*kappa1);

I_upstream.phi = Q0*X_upstream.*Y_upstream.phi.';
I_upstream.dphi = Q0*X_upstream.*Y_upstream.dphi.';
I_upstream.p = Q0*X_upstream.*Y_upstream.p.';
I_upstream.v = Q0*X_upstream.*Y_upstream.v.';

% Some scalings for the scattered thing:
Wall_Upstream = Y_GreensHard(0,y2,kappa1,k3,omega,U,c,delta,N_bly,-1);
Pressure_Upstream = Q0*Wall_Upstream.p; % i.e. phid(y2)*ph(0): easy enough to genericise

%% X-dependence - differs for the two scattering problems

% Upstream scattering x1: translated so x1_u = 0 is upstream junction
x1_U = x1 - x1_u;
Xp_u = Fourier_Nov2018(x1_U(x1_U<0),Cp.k);
Xm_u = Fourier_Nov2018(x1_U(x1_U>=0),Cm.k);
X_CLp_u = Fourier_Nov2018(x1_U(x1_U>=0),C_CLp.k);
X_CLm_u = Fourier_Nov2018(x1_U(x1_U>=0),C_CLm.k);

% Downstream scattering: similarly
x1_D = x1 - x1_d;
Xp_d = Fourier_Nov2018(x1_D(x1_D<0),Cp.k);
Xm_d = Fourier_Nov2018(x1_D(x1_D>=0),Cm.k);
X_CLp_d = Fourier_Nov2018(x1_D(x1_D>=0),C_CLp.k);
X_CLm_d = Fourier_Nov2018(x1_D(x1_D>=0),C_CLm.k);


%% Y-dependence - generic for all scattering
% A rate-limiting routine, particularly if x2 is large. Fortunately,
% appropriate for all integrals.
disp('Computing wall-normal behaviour');
[Yp] = Y_Continuous_Scattering_Nov2018_Mod(x2,Cp.k,k3,omega,U,c,delta,N_bl);
[Ym] = Y_Continuous_Scattering_Nov2018_Mod(x2,Cm.k,k3,omega,U,c,delta,N_bl);
[Y_CLp] = Y_Continuous_Scattering_CL_Nov2018_Mod(x2,C_CLp.k,k3,omega,U,c,delta,N_bl,1);
[Y_CLm] = Y_Continuous_Scattering_CL_Nov2018_Mod(x2,C_CLm.k,k3,omega,U,c,delta,N_bl,-1);
disp('Computed wall-normal behaviour');

%% Dispersion relationships
% Contours: z0 along real axis, zp and zm semicircles
z0.f = @(s) delta*(1-s);
z0.df = @(s) -delta;
zp = Contour_Semicircle(delta,1);
zm = Contour_Semicircle(delta,-1);

% Dispersion relationships without deformation
decp = Comp_Dec_Aug18A_skip(Cp.k,k3,omega,U,c,delta,N_bl,z0);
decm = Comp_Dec_Aug18A_skip(Cm.k,k3,omega,U,c,delta,N_bl,z0);
Dh_p = decp.v;
Dh_m = decm.v;
Dp_p = decp.p;
Dp_m = decm.p;

% Dispersion relationships with deformation
dec_CLp = Comp_Dec_Aug18A_skip(C_CLp.k,k3,omega,U,c,delta,N_bl,zp);
dec_CLm = Comp_Dec_Aug18A_skip(C_CLm.k,k3,omega,U,c,delta,N_bl,zm);
Dh_CLp = dec_CLp.v;
Dp_CLp = dec_CLp.p;
Dh_CLm = dec_CLm.v;
Dp_CLm = dec_CLm.p;

% Specific points of interest (note already have dDp_kappap from above)
dec_kappap = Comp_Dec_Aug18A_skip(kappap,k3,omega,U,c,delta,N_bl,zp);
Dh_kappap = dec_kappap.v;

%% Kernel
% Kernel is the "same" for both scattering problems, which helps somewhat
z_kernel = Contour_Semicircle(delta,1); % So branch cut in LHP, as with zp earlier

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
Kernel = @(k1) Kernel_HardSoft(k1,k3,omega,U,c,delta,N_bl,z_kernel);

% Far-field kernel, with allowance for non-slip.
if U.f(0) ~= 0
    Kernel_ff = FFKernel_HardSoft_Slip(omega,U.f(0),c.f(0),k3,0);
else
    Kernel_ff = FFKernel_HardSoft_NoSlip(omega,U.df(0),c.f(0),k3,0);
end

[K_factorised] = Fact_Mult_Cont_Scaled_Finite(Kernel,Kernel_ff,contour,prec,branch_cut);
Kernel_p_p = K_factorised.p(Cp.k);
Kernel_m_m = K_factorised.m(Cm.k);

Kernel_m_CLp = K_factorised.m(C_CLp.k);
Kernel_m_CLm = K_factorised.m(C_CLm.k);

% Some points require specific evaluation. All are in the generalised LHP,
% so all use Km.
Km_kappa1 = K_factorised.m(kappa1);
Kp_kappa1 = Kernel(kappa1)./Km_kappa1;
Km_kappap = K_factorised.m(kappap);

Km_zeds = K_factorised.m(zeds);
disp('Kernel found and factorised, setting up integrals')

%% Cauchy-like terms
Cauchy_kappa1_p = 1./(Cp.k - kappa1);
Cauchy_kappa1_m = 1./(Cm.k - kappa1);
Cauchy_kappa1_CLp = 1./(C_CLp.k - kappa1);
Cauchy_kappa1_CLm = 1./(C_CLm.k - kappa1);

Cauchy_kappap_p = 1./(Cp.k - kappap);
Cauchy_kappap_m = 1./(Cm.k - kappap);
Cauchy_kappap_CLp = 1./(C_CLp.k - kappap);
Cauchy_kappap_CLm = 1./(C_CLm.k - kappap);

Cauchy_kappap_kappa1 = 1./(kappap - kappa1);


%% Scaling (inc l-scaling)
% Note: the 2pi is in the Fourier term.
A_s1 = -Pressure_Upstream./(1i*Km_kappa1); 
A_s2_1 = ( Pressure_Upstream.*Dh_kappap.*Km_kappap.^2.*exp(-1i*kappap*l_junction) ) ...
    ./ (1i*Km_kappa1.*(kappap - kappa1).*dDp_kappap);
A_s2_2 = ( Pressure_Upstream*exp(-1i*kappa1*l_junction) )./(1i*Kp_kappa1);

%% K-combination
% Upstream junction
Kp_s1 = A_s1*Cauchy_kappa1_p./(Kernel_p_p.*Dh_p);
Km_s1 = A_s1*Cauchy_kappa1_m.*Kernel_m_m./Dp_m;
K_CLp_s1 = A_s1*Cauchy_kappa1_CLp.*Kernel_m_CLp./Dp_CLp;
K_CLm_s1 = A_s1*Cauchy_kappa1_CLm.*Kernel_m_CLm./Dp_CLm;

% Downstream junction: modal unscattering
Kp_s2_1 = A_s2_1*Cauchy_kappap_p.*Kernel_p_p./(Dp_p);
Km_s2_1 = A_s2_1*Cauchy_kappap_m./(Kernel_m_m.*Dh_m);
K_CLp_s2_1 = A_s2_1*Cauchy_kappap_CLp./(Kernel_m_CLp.*Dh_CLp);
K_CLm_s2_1 = A_s2_1*Cauchy_kappap_CLm./(Kernel_m_CLm.*Dh_CLm);

% Downstream junction: forcing term
Kp_s2_2 = A_s2_2*Cauchy_kappa1_p.*Kernel_p_p./(Dp_p);
Km_s2_2 = A_s2_2*Cauchy_kappa1_m./(Kernel_m_m.*Dh_m);
K_CLp_s2_2 = A_s2_2*Cauchy_kappa1_CLp./(Kernel_m_CLp.*Dh_CLp);
K_CLm_s2_2 = A_s2_2*Cauchy_kappa1_CLm./(Kernel_m_CLm.*Dh_CLm);


%% Integration: acoustic

disp('Computing steepest descent integral');
I_generic = struct('phi',[],'dphi',[],'p',[],'v',[]);
physical_variables = fieldnames(I_generic);

I_s1_p = I_generic; % Scattering due to first junction, upstream of first junction
I_s1_m = I_generic;

I_s2_1_p = I_generic ; % Upstream of second junction
I_s2_1_m = I_generic ;
I_s2_2_p = I_generic ;
I_s2_2_m = I_generic ;

% Summations of the above
I_s1 = I_generic;
I_s2_1 = I_generic;
I_s2_2 = I_generic;
I_acoustic  = I_generic;

% This loops through all physical variables, which at least cleans up the
% code, even if speed isn't a particular characteristic of this code.
for j = 1:numel(physical_variables)
    s = physical_variables{j};
    I_s1_p.(s) = Integral_Nov2018(Cp.wk,Kp_s1,Xp_u,Yp.(s));
    I_s1_m.(s) = Integral_Nov2018(Cm.wk,Km_s1,Xm_u,Ym.(s));
    I_s2_1_p.(s) = Integral_Nov2018(Cp.wk,Kp_s2_1,Xp_d,Yp.(s));
    I_s2_1_m.(s) = Integral_Nov2018(Cm.wk,Km_s2_1,Xm_d,Ym.(s));
    I_s2_2_p.(s) = Integral_Nov2018(Cp.wk,Kp_s2_2,Xp_d,Yp.(s));
    I_s2_2_m.(s) = Integral_Nov2018(Cm.wk,Km_s2_2,Xm_d,Ym.(s));
    I_s1.(s) = [I_s1_p.(s);I_s1_m.(s)];
    I_s2_1.(s) = [I_s2_1_p.(s);I_s2_1_m.(s)];
    I_s2_2.(s) = [I_s2_2_p.(s);I_s2_2_m.(s)];
    I_acoustic.(s) = I_s1.(s) + I_s2_1.(s) + I_s2_1.(s);
end

%% Integration: critical-layer
disp('Computing critical-layer integral');

% Again, a large amount of different functions floating around
I_CL_s1_p = I_generic;
I_CL_s1_m = I_generic;
I_CL_s1 = I_generic;

I_CL_s2_1_p = I_generic;
I_CL_s2_1_m = I_generic;
I_CL_s2_1 = I_generic;

I_CL_s2_2_p = I_generic;
I_CL_s2_2_m = I_generic;
I_CL_s2_2 = I_generic;

I_CL = I_generic;

for j = 1:numel(physical_variables)
    s = physical_variables{j};
    I_CL_s1_p.(s) = Integral_Nov2018(C_CLp.wk,K_CLp_s1,X_CLp_u,Y_CLp.(s));
    I_CL_s1_m.(s) = Integral_Nov2018(C_CLm.wk,K_CLm_s1,X_CLm_u,Y_CLm.(s));
    I_CL_s2_1_p.(s) = Integral_Nov2018(C_CLp.wk,K_CLp_s2_1,X_CLp_d,Y_CLp.(s));
    I_CL_s2_1_m.(s) = Integral_Nov2018(C_CLm.wk,K_CLm_s2_1,X_CLm_d,Y_CLm.(s));
    I_CL_s2_2_p.(s) = Integral_Nov2018(C_CLp.wk,K_CLp_s2_2,X_CLp_d,Y_CLp.(s));
    I_CL_s2_2_m.(s) = Integral_Nov2018(C_CLm.wk,K_CLm_s2_2,X_CLm_d,Y_CLm.(s));
    I_CL_s1.(s) = [zeros(size(I_s1_p.p));I_CL_s1_p.(s) + I_CL_s1_m.(s)];
    I_CL_s2_1.(s) = [zeros(size(I_s2_1_p.p));I_CL_s2_1_p.(s) + I_CL_s2_1_m.(s)];
    I_CL_s2_2.(s) = [zeros(size(I_s2_2_p.p));I_CL_s2_2_p.(s) + I_CL_s2_2_m.(s)];
    I_CL.(s) = I_CL_s1.(s) + I_CL_s2_1.(s) + I_CL_s2_1.(s);
end


%% Modal contributions

disp('Computing modal contributions');
% Two modal contributions, only in the "middle" region, from zeros of Dp
x1_mid = x1( (x1 >= x1_u) & (x1 < x1_d) ) - x1_u;

X_mid = Fourier_Nov2018(x1_mid,zeds);
Y_mid = Y_Continuous_Scattering_Nov2018_Mod(x2,zeds,k3,omega,U,c,delta,N_bl);

Res_weight = zeros(size(zeds));
for j = 1:numel(zeds)
    Res_weight(j) = -2*pi*1i*A_s1*Km_zeds(j)./((zeds(j) - kappa1).*dDp_zeds(j));
end

% Assuming two zeros, can probably be generalised
I_mid_1 = I_generic;
I_mid_2 = I_generic;
I_mid = I_generic;
I_res = I_generic;
for j = 1:numel(physical_variables)
    s = physical_variables{j};
    I_mid_1.(s) = (Res_weight(1)*Y_mid.(s)(:,1)*X_mid(:,1).').';
    I_mid_2.(s) = (Res_weight(2)*Y_mid.(s)(:,2)*X_mid(:,2).').';
    I_mid.(s) = I_mid_1.(s) + I_mid_2.(s);
    I_res.(s) = [zeros(size(I_s1_p.p));I_mid.(s);zeros(size(I_s2_1_m.p))];
end


%% Combinations (other combinations are available)
disp('Adding it all together');
% Full scattered solution
I_scattered = I_generic;

% Full upstream+scattered solution
I_full = I_generic;

% Only the first scattered solution (with and without incident, a consistency check)
I_scattered1 = I_generic;
I_full1 = I_generic;

for j = 1:numel(physical_variables)
    s = physical_variables{j};
    I_scattered.(s) = I_res.(s) + I_CL.(s) + I_acoustic.(s);
    I_full.(s) = I_upstream.(s) + I_scattered.(s);
    I_scattered1.(s) = I_res.(s) + I_CL_s1.(s) + I_s1.(s);
    I_full1.(s) = I_upstream.(s) + I_scattered1.(s);
end

disp('All done - well done');