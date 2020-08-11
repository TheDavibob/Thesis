%% Hopefully pings out the various things we want for a (slipping or not) hard/soft transition

%% Parameter setup
delta = 1;
%U = Struct_U_Linear(0,1,delta);
U = Struct_U_Parabolic(0,1,delta);
c = Struct_c_Const(5);
omega = 1;
k3 = 0;
N_bl = 1000;

Zu = 1-1i; % Can be infininty
Zd = 2*(1-1i); % As above
[Lu,Ld] = ImpedanceUpDown(Zu,Zd,omega,U.f(0),U.df(0));
% Generates boundary condition operators, as Lu.l1, Lu.l2 etc.

y2 = 0.5;
kappa1 = omega./U.f(y2);

x1 = -5:0.02:5;
x2 = 0:0.02:5;
% x2 = [0,0];

%% Contour setup
N_C = 1000;

Minf = U.f(delta)./c.f(delta);
k0 = omega./c.f(delta);
kb1 = -Minf/(1-Minf^2)*k0 + sqrt(k0^2 + k3^2*(1-Minf^2))./(1-Minf^2);
kb2 = -Minf/(1-Minf^2)*k0 - sqrt(k0^2 + k3^2*(1-Minf^2))./(1-Minf^2);

sep = 0.1;

%Cp = Finite_Parabola_3Points(N_C,kb2*(1+sep),kb2*(1+1i*sep),kb2*(1-sep));
%Cm = Finite_Parabola_3Points(N_C,kb1*(1-sep),kb1*(1+1i*sep),kb1*(1+sep));

angle = pi/4;
curvature = 0.01;
apex = 0.1*omega;
Cp = Finite_Symmetric_Hyperbola(N_C,kb2,-apex,curvature,angle);
Cm = Finite_Symmetric_Hyperbola(N_C,kb1,apex,curvature,angle);

%Cm = Finite_Parabola_3Points(N_C,0,0.75*omega/U.f(0) + 0.2i,1.5*omega/U.f(0));

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

%% Upstream terms (which, in turn, will fix the scalings
Q0 = 1e-10; % Depends on y2 in a variety of reasonably subtle ways, but scales out

% Useful previous work: easy to replace with Generic, too. Note:
% deformation of critical layer into lower HP to "preserve causality" in
% some bollox sense.
N_bly = 10*N_bl;

Y_upstream = Y_GreensGeneric(x2,y2,kappa1,k3,omega,U,c,Lu.l1,Lu.l2,delta,N_bly,-1);
X_upstream = exp(-1i*x1.'*kappa1);

I_upstream.phi = Q0*X_upstream.*Y_upstream.phi.';
I_upstream.dphi = Q0*X_upstream.*Y_upstream.dphi.';
I_upstream.p = Q0*X_upstream.*Y_upstream.p.';
I_upstream.v = Q0*X_upstream.*Y_upstream.v.';

% Some scalings for the scattered thing:
Wall_Upstream = Y_GreensGeneric(0,y2,kappa1,k3,omega,U,c,Lu.l1,Lu.l2,delta,N_bly,-1);
Incident_Upstream = Q0*( Ld.l1(kappa1).*Wall_Upstream.dphi + Ld.l2(kappa1).*Wall_Upstream.phi );
% i.e. phid(y2)*Ld(phi)(0)


%% X-depdendence
Xp = Fourier_Nov2018(x1(x1<0),Cp.k);
Xm = Fourier_Nov2018(x1(x1>=0),Cm.k);
%X = [Xp;Xm];

X_CLp = Fourier_Nov2018(x1(x1>=0),C_CLp.k);
X_CLm = Fourier_Nov2018(x1(x1>=0),C_CLm.k);

%% Y-dependence
% Decidedly not quick
[Yp] = Y_Continuous_Scattering_Nov2018_Mod(x2,Cp.k,k3,omega,U,c,delta,N_bl);
[Ym] = Y_Continuous_Scattering_Nov2018_Mod(x2,Cm.k,k3,omega,U,c,delta,N_bl);
%Y = [Yp;Ym];

[Y_CLp] = Y_Continuous_Scattering_CL_Nov2018_Mod(x2,C_CLp.k,k3,omega,U,c,delta,N_bl,1);
[Y_CLm] = Y_Continuous_Scattering_CL_Nov2018_Mod(x2,C_CLm.k,k3,omega,U,c,delta,N_bl,-1);

%% Dispersion relationship
% Going to recompute, as it makes it easier to flip any boundary conditions
clear z
z.f = @(s) delta*(1-s);
z.df = @(s) -delta;
decp = Comp_Dec_Aug18A_skip(Cp.k,k3,omega,U,c,delta,N_bl,z);
Du_p = Lu.l1(Cp.k).*decp.dphi + Lu.l2(Cp.k).*decp.phi;
decm = Comp_Dec_Aug18A_skip(Cm.k,k3,omega,U,c,delta,N_bl,z);
% Dh_m = decm.v;
Dd_m = Ld.l1(Cm.k).*decm.dphi + Ld.l2(Cm.k).*decm.phi;

zp = Contour_Semicircle(delta,1);
zm = Contour_Semicircle(delta,-1);
dec_CLp = Comp_Dec_Aug18A_skip(C_CLp.k,k3,omega,U,c,delta,N_bl,zp);
% Dh_CLp = dec_CLp.v;
Dd_CLp = Ld.l1(C_CLp.k).*dec_CLp.dphi + Ld.l2(C_CLp.k).*dec_CLp.phi;
dec_CLm = Comp_Dec_Aug18A_skip(C_CLm.k,k3,omega,U,c,delta,N_bl,zm);
% Dh_CLm = dec_CLm.v;
Dd_CLm = Ld.l1(C_CLm.k).*dec_CLm.dphi + Ld.l2(C_CLm.k).*dec_CLm.phi;


%% Forcing pole
pole_p = 1./(Cp.k - kappa1);
pole_m = 1./(Cm.k - kappa1);

pole_CLp = 1./(C_CLp.k - kappa1);
pole_CLm = 1./(C_CLm.k - kappa1);

%% Kernel
%z_kernel.f = @(s) delta*(1-s);
%z_kernel.df = @(s) -delta*s ;
z_kernel = Contour_Semicircle(delta,1); % So branch cut in LHP

%z_down = Contour_Semicircle(delta,-1); % For finding kernel at kappa1.
clear prec
prec.N = 2000;
%[contour.Cp,contour.dCp] = TanhContour(-sep,3*sep,0.5*(kb1 + kb2),kb1);
%[contour.Cm,contour.dCm] = TanhContour(sep,3*sep,kb2,0.5*(kb1 + kb2));
[contour.Cp,contour.dCp] = TanhContour(0,1,-1,1);
[contour.Cm,contour.dCm] = TanhContour(0,1,-1,1);

% branch_cut = pi/2;
branch_cut = pi;

% Potentially going to need "curved" log cut, eg.
% branch_cut = @(r) pi/2 * ( r - 1)./(r + 1) - 2*pi;

Kernel = @(k1) Kernel_Generic(k1,k3,omega,U,c,delta,Lu,Ld,N_bl,z_kernel);
% Kernel_down = @(k1) Kernel_HardSoft(k1,k3,omega,U,c,delta,N_bl,z_down);

if U.f(0) ~= 0
    Kernel_ff = FFKernel_Impedance(Zu,Zd,omega,U.f(0),c.f(0),k3,0);
else
    Kernel_ff = FFKernel_Impedance(Zu,Zd,omega,U.f(0),c.f(0),k3,0); % This should be fine.
end

[K_factorised] = Fact_Mult_Cont_Scaled_Finite(Kernel,Kernel_ff,contour,prec,branch_cut);
Kernel_p_p = K_factorised.p(Cp.k);
%Kernel_p_m = K_factorised.f(Cm.k)./K_factorised.m(Cm.k);
%Kernel_p_m =  K_factorised.p(Cm.k);
Kernel_m_m = K_factorised.m(Cm.k); % Going to use this one

%Kernel_p_CLp = K_factorised.f(C_CLp.k)./K_factorised.m(C_CLp.k);
%Kernel_p_CLm = K_factorised.f(C_CLm.k)./K_factorised.m(C_CLm.k);
%Kernel_p_CLp = K_factorised.p(C_CLp.k);
%Kernel_p_CLm = K_factorised.p(C_CLm.k);
Kernel_m_CLp = K_factorised.m(C_CLp.k);
Kernel_m_CLm = K_factorised.m(C_CLm.k);

Km_kappa1 = K_factorised.m(kappa1);


%% Scaling bits
A = -1i*Incident_Upstream./Km_kappa1;

%% Combination of k dependent bits
Kp = A*pole_p./(Kernel_p_p.*Du_p);
%Km = A*pole_m./(Kernel_p_m.*Dh_m);
Km = A*pole_m.*Kernel_m_m./Dd_m; % For analyticity reasons, this is the best one.

% K_CLp = A*pole_CLp./(Kernel_p_CLp.*Dh_CLp);
% K_CLm = A*pole_CLm./(Kernel_p_CLm.*Dh_CLm);
K_CLp = A*pole_CLp.*Kernel_m_CLp./Dd_CLp;
K_CLm = A*pole_CLm.*Kernel_m_CLm./Dd_CLm;

%% Integration
clear Ip Im I
Ip.phi = Integral_Nov2018(Cp.wk,Kp,Xp,Yp.phi);
Im.phi = Integral_Nov2018(Cm.wk,Km,Xm,Ym.phi);
Ip.dphi = Integral_Nov2018(Cp.wk,Kp,Xp,Yp.dphi);
Im.dphi = Integral_Nov2018(Cm.wk,Km,Xm,Ym.dphi);
Ip.p = Integral_Nov2018(Cp.wk,Kp,Xp,Yp.p);
Im.p = Integral_Nov2018(Cm.wk,Km,Xm,Ym.p);
Ip.v = Integral_Nov2018(Cp.wk,Kp,Xp,Yp.v);
Im.v = Integral_Nov2018(Cm.wk,Km,Xm,Ym.v);

I.phi = [Ip.phi;Im.phi];
I.dphi = [Ip.dphi;Im.dphi];
I.p = [Ip.p;Im.p];
I.v = [Ip.v;Im.v];

%% Critical layer integration

I_CLp.phi = Integral_Nov2018(C_CLp.wk,K_CLp,X_CLp,Y_CLp.phi);
I_CLm.phi = Integral_Nov2018(C_CLm.wk,K_CLm,X_CLm,Y_CLm.phi);
I_CL.phi = I_CLp.phi + I_CLm.phi;

I_CLp.dphi = Integral_Nov2018(C_CLp.wk,K_CLp,X_CLp,Y_CLp.dphi);
I_CLm.dphi = Integral_Nov2018(C_CLm.wk,K_CLm,X_CLm,Y_CLm.dphi);
I_CL.dphi = I_CLp.dphi + I_CLm.dphi;

I_CLp.p = Integral_Nov2018(C_CLp.wk,K_CLp,X_CLp,Y_CLp.p);
I_CLm.p = Integral_Nov2018(C_CLm.wk,K_CLm,X_CLm,Y_CLm.p);
I_CL.p = I_CLp.p + I_CLm.p;

I_CLp.v = Integral_Nov2018(C_CLp.wk,K_CLp,X_CLp,Y_CLp.v);
I_CLm.v = Integral_Nov2018(C_CLm.wk,K_CLm,X_CLm,Y_CLm.v);
I_CL.v = I_CLp.v + I_CLm.v;


%% Adding it all together, again
I_tot.phi = [Ip.phi; Im.phi + I_CL.phi];
I_tot.dphi = [Ip.dphi; Im.dphi + I_CL.dphi];
I_tot.p = [Ip.p; Im.p + I_CL.p];
I_tot.v = [Ip.v; Im.v + I_CL.v];

%% Upstream + Downstream
I_full.phi = I_tot.phi + I_upstream.phi;
I_full.dphi = I_tot.dphi + I_upstream.dphi;
I_full.p = I_tot.p + I_upstream.p;
I_full.v = I_tot.v + I_upstream.v;
