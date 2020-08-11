%% Computes the disturbances due to a point source above a discontinuous wake, by utilising the old routines

%% Parameters
deltap = 1;
U0 = 0.3;
Uinfp = 1;
% U = Struct_U_Linear(0.3,1,delta);
Up = Struct_U_Parabolic(U0,Uinfp,deltap);
cp = Struct_c_Const(5);

deltam = 1;
Uinfm = 1;
Um = Struct_U_Parabolic(U0,Uinfm,deltam);
cm = Struct_c_Const(5);

omega = 1;
k3 = 0;
N_bl = 1000;

q0 = 1;

y2 = 0.5;
kappa1 = omega./Up.f(y2);

x1 = -2:0.05:10;
x2 = -2:0.05:2;
%x2 = [0.4,0.5,0.6];
% x2 = [0.5,0.5];

%% Boundary conditions for "Upper" solution

C0 = @(k1) 1i*(omega - Up.f(0)*k1);
dC0p = @(k1) -1i*Up.df(0)*k1;
dC0m = @(k1) -1i*Um.df(0)*k1;


z_below.f = @(s) deltam*(1-s);
z_below.df = @(s) - deltam*ones(size(s));
dec_below_phi = @(k1) WakeDec_phi(k1,k3,omega,Um,cm,deltam,N_bl,z_below) ; % These two addn functions might be unnecessary.
dec_below_dphi = @(k1) WakeDec_dphi(k1,k3,omega,Um,cm,deltam,N_bl,z_below) ;

l1 = @(k1) dec_below_phi(k1);
l2 = @(k1) dec_below_dphi(k1) + 3*(dC0m(k1) + dC0p(k1))./C0(k1).*dec_below_phi(k1);

%% Scaling constant
A = 1; % to be fixed

%% Contour setup
N_C = 200;

Minf = Up.f(deltap)./cp.f(deltap);
k0 = omega./cp.f(deltap);
kb1 = -Minf/(1-Minf^2)*k0 + sqrt(k0^2 + k3^2*(1-Minf^2))./(1-Minf^2);
kb2 = -Minf/(1-Minf^2)*k0 - sqrt(k0^2 + k3^2*(1-Minf^2))./(1-Minf^2);

sep = 0.1;

%Cp = Finite_Parabola_3Points(N_C,kb2*(1+sep),kb2*(1+1i*sep),kb2*(1-sep));
%Cm = Finite_Parabola_3Points(N_C,kb1*(1-sep),kb1*(1+1i*sep),kb1*(1+sep));

angle = pi/4;
curvature = 0.01;
apex = 0.1;
Cp = Finite_Symmetric_Hyperbola(N_C,kb2,-apex,curvature,angle);
Cm = Finite_Symmetric_Hyperbola(N_C,kb1,apex,curvature,angle);

%Cm = Finite_Parabola_3Points(N_C,0,0.75*omega/U.f(0) + 0.2i,1.5*omega/U.f(0));

%% Critical layer contours

N_CL = 1000; % More precision needed, :(
a = omega./Up.f(deltap);
b = omega./Up.f(0);
%shift1 = 0.4;
%shift2 = 0.4;
shift1 = 0.1*omega;
shift2 = 0.1*omega;

[C_CLp,C_CLm] = Finite_CL_Ellipse(N_CL,a,b,shift1,shift2);

%% Extra, if zeros known
D_l = @(k1) Dispersion_Generic_Nov2018(k1,k3,omega,Up,cp,l1,l2,deltap,N_bl,0);
smol = 1e-12;
dD_l = @(k1) (D_l(k1 + smol) - D_l(k1 - smol)) ./(2*smol);

zeds(1) = NewtonRaphson(D_l,dD_l,1.5+0.4i);
zeds(2) = NewtonRaphson(D_l,dD_l,1.5-0.4i);

% zeds = [];
dispersion_zeros = zeds; %
zeros_UHP = zeds(imag(zeds)>0);
zeros_LHP = zeds(imag(zeds)<=0);

N_zeros_UHP = numel(zeros_UHP);
N_zeros_LHP = numel(zeros_LHP);

step = 0.01;
radius = 0.01;

for j = 1:N_zeros_UHP
    C_pole.k = zeros_UHP(j) + radius*exp(-1i*(0:step:2*pi));
    C_pole.wk = -step*radius*1i*exp(-1i*(0:step:2*pi));

    C_CLp.k = [C_CLp.k,C_pole.k];
    C_CLp.wk = [C_CLp.wk,C_pole.wk];
end

for j = 1:N_zeros_LHP
    C_pole.k = zeros_LHP(j) + radius*exp(-1i*(0:step:2*pi));
    C_pole.wk = -step*radius*1i*exp(-1i*(0:step:2*pi));

    C_CLm.k = [C_CLm.k,C_pole.k];
    C_CLm.wk = [C_CLm.wk,C_pole.wk];
end


%% X-depdendence
Xp = Fourier_Nov2018(x1(x1<0),Cp.k);
Xm = Fourier_Nov2018(x1(x1>=0),Cm.k);
%X = [Xp;Xm];

X_CLp = Fourier_Nov2018(x1(x1>=0),C_CLp.k);
X_CLm = Fourier_Nov2018(x1(x1>=0),C_CLm.k);

%% Y-dependence
% Decidedly not quick
[Yp] = Y_GreensWake(x2,y2,Cp.k,k3,omega,Up,cp,Um,cm,l1,l2,deltap,deltam,N_bl,0);
[Ym] = Y_GreensWake(x2,y2,Cm.k,k3,omega,Up,cp,Um,cm,l1,l2,deltap,deltam,N_bl,0);

[Y_CLp] = Y_GreensWake(x2,y2,C_CLp.k,k3,omega,Up,cp,Um,cm,l1,l2,deltap,deltam,N_bl,1);
[Y_CLm] = Y_GreensWake(x2,y2,C_CLm.k,k3,omega,Up,cp,Um,cm,l1,l2,deltap,deltam,N_bl,-1);

%% k-dependence
Kp =  K_Generic(y2,Cp.k,k3,omega,Up,cp,l1,l2,deltap,N_bl,0);
Km =  K_Generic(y2,Cm.k,k3,omega,Up,cp,l1,l2,deltap,N_bl,0);

K_CLp = K_Generic(y2,C_CLp.k,k3,omega,Up,cp,l1,l2,deltap,N_bl,1);
K_CLm = K_Generic(y2,C_CLm.k,k3,omega,Up,cp,l1,l2,deltap,N_bl,-1);

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
