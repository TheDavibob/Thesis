%% Hopefully pings out the various things we want for a (slipping or not) hard/soft transition

%% Parameter setup
delta = 1;
%U = Struct_U_Linear(0,1,delta);
U = Struct_U_Parabolic(0.3,1,delta);
c = Struct_c_Const(5);
omega = 20;
k3 = 0;
N_bl = 1000;

y2 = 0.1;
kappa1 = omega./U.f(y2);

x1 = -5:0.05:5;
x2 = 0:0.05:5;
% x2 = [0.2,0.3];

%% Contour setup
N_C = 1000;

Minf = U.f(delta)./c.f(delta);
k0 = omega./c.f(delta);
kb1 = -Minf/(1-Minf^2)*k0 + sqrt(k0^2 + k3^2*(1-Minf^2))./(1-Minf^2);
kb2 = -Minf/(1-Minf^2)*k0 - sqrt(k0^2 + k3^2*(1-Minf^2))./(1-Minf^2);

sep = 0.1;

%Cp = Finite_Parabola_3Points(N_C,kb2*(1+sep),kb2*(1+1i*sep),kb2*(1-sep));
%Cm = Finite_Parabola_3Points(N_C,kb1*(1-sep),kb1*(1+1i*sep),kb1*(1+sep));

angle = 3*pi/8;
curvature = 0.01;
apex = 0.1*omega;
Cp = Finite_Symmetric_Hyperbola(N_C,kb2,-apex,curvature,angle);
Cm = Finite_Symmetric_Hyperbola(N_C,kb1,apex,curvature,angle);

%Cm = Finite_Parabola_3Points(N_C,0,0.75*omega/U.f(0) + 0.2i,1.5*omega/U.f(0));

%% Critical layer integral
a = omega./U.f(delta);
if U.f(0)~= 0
    b = omega./U.f(0);
else
    b = 100*a; % Currently a fudge: this crosses the critical layer.
end
shift1 = 0.5;
shift2 = 0.5;

[C_CLp,C_CLm] = Finite_CL_Ellipse(N_C,a,b,shift1,shift2);

% a fudge, a test, please delet soon
% dispersion_zeros = [1.2202 + 0.5378i,1.2202 - 0.5378i]; % Only parabolic, boring frequency, etc. etc.
% C_pole1.k = dispersion_zeros(1) + 0.01*exp(-1i*(0:0.01:2*pi));
% C_pole1.wk = -0.01*0.01*1i*exp(-1i*(0:0.01:2*pi));
% 
% C_pole2.k = dispersion_zeros(2) + 0.01*exp(-1i*(0:0.01:2*pi));
% C_pole2.wk = -0.01*0.01*1i*exp(-1i*(0:0.01:2*pi));
% 
% C_CLp.k = [C_CLp.k,C_pole1.k];
% C_CLp.wk = [C_CLp.wk,C_pole1.wk];
% 
% C_CLm.k = [C_CLm.k,C_pole2.k];
% C_CLm.wk = [C_CLm.wk,C_pole2.wk];


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
Dh_p = decp.v;
decm = Comp_Dec_Aug18A_skip(Cm.k,k3,omega,U,c,delta,N_bl,z);
% Dh_m = decm.v;
Dp_m = decm.p;

zp = Contour_Semicircle(delta,1);
zm = Contour_Semicircle(delta,-1);
dec_CLp = Comp_Dec_Aug18A_skip(C_CLp.k,k3,omega,U,c,delta,N_bl,zp);
% Dh_CLp = dec_CLp.v;
Dp_CLp = dec_CLp.p;
dec_CLm = Comp_Dec_Aug18A_skip(C_CLm.k,k3,omega,U,c,delta,N_bl,zm);
% Dh_CLm = dec_CLm.v;
Dp_CLm = dec_CLm.p;


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

branch_cut = pi/2;
% branch_cut = pi;

% Potentially going to need "curved" log cut, eg.
% branch_cut = @(r) pi/2 * ( r - 1)./(r + 1) - 2*pi;

Kernel = @(k1) Kernel_HardSoft(k1,k3,omega,U,c,delta,N_bl,z_kernel);
%Kernel_down = @(k1) Kernel_HardSoft(k1,k3,omega,U,c,delta,N_bl,z_down);

if U.f(0) ~= 0
    Kernel_ff = FFKernel_HardSoft_Slip(omega,U.f(0),c.f(0),k3,0);
else
    Kernel_ff = FFKernel_HardSoft_NoSlip(omega,U.df(0),c.f(0),k3,0);
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
A = -1i*Pressure_Upstream./Km_kappa1;

%% Combination of k dependent bits
Kp = A*pole_p./(Kernel_p_p.*Dh_p);
%Km = A*pole_m./(Kernel_p_m.*Dh_m);
Km = A*pole_m.*Kernel_m_m./Dp_m; % For analyticity reasons, this is the best one.

% K_CLp = A*pole_CLp./(Kernel_p_CLp.*Dh_CLp);
% K_CLm = A*pole_CLm./(Kernel_p_CLm.*Dh_CLm);
K_CLp = A*pole_CLp.*Kernel_m_CLp./Dp_CLp;
K_CLm = A*pole_CLm.*Kernel_m_CLm./Dp_CLm;

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

%% Residues at poles, in LHk1P only
% smol = 1e-11;
% if exist('dispersion_zeros','var') == 1
%     dec_dispersion_zeros_plus = Comp_Dec_Aug18A_skip(dispersion_zeros + smol,k3,omega,U,c,delta,N_bl,z);
%     dec_dispersion_zeros_minus = Comp_Dec_Aug18A_skip(dispersion_zeros - smol,k3,omega,U,c,delta,N_bl,z);
%     dDp_dispersion_zeros = (dec_dispersion_zeros_plus.p - dec_dispersion_zeros_minus.p)./(2*smol);
%     
%     poles_dispersion_zeros = 1./(dispersion_zeros - kappa1);
%     Km_dispersion_zeros = K_factorised.m(dispersion_zeros);
%     
%     K_dispersion_zeros = -2*pi*1i*A.*Km_dispersion_zeros.*poles_dispersion_zeros./dDp_dispersion_zeros ;
%     % Note traversal is clockwise
%     
%     [Y_dispersion_zeros] = Y_Continuous_Scattering_Nov2018_Mod(x2,dispersion_zeros,k3,omega,U,c,delta,N_bl);
%     
%     X_dispersion_zeros = Fourier_Nov2018(x1(x1>=0),dispersion_zeros);
%     
%     clear I_res
%     for j = 1:numel(dispersion_zeros)
%         I_res.phi(:,:,j) = X_dispersion_zeros(:,j)*Y_dispersion_zeros.phi(:,j).';
%         I_res.dphi(:,:,j) = X_dispersion_zeros(:,j)*Y_dispersion_zeros.dphi(:,j).';
%         I_res.p(:,:,j) = X_dispersion_zeros(:,j)*Y_dispersion_zeros.p(:,j).';
%         I_res.v(:,:,j) = X_dispersion_zeros(:,j)*Y_dispersion_zeros.v(:,j).';
%         
%         I_CL.phi = I_CL.phi + I_res.phi(:,:,j);
%         I_CL.dphi = I_CL.dphi + I_res.dphi(:,:,j);
%         I_CL.p = I_CL.p + I_res.p(:,:,j);
%         I_CL.v = I_CL.v + I_res.v(:,:,j);
%     end
% end

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
