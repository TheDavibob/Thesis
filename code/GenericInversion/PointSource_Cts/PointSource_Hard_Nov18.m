%% Computes the disturbances due to a point source above a hard wall

%% Parameters
delta = 1;
% U = Struct_U_Linear(0.2,1,delta);
% U = Struct_U_Parabolic(0.2,1,delta);
U = Struct_U_Linear(0.2,1,delta);
c = Struct_c_Const(5);
omega = 0.1;
k3 = 0;
N_bl = 1000;

q0 = 1;

y2 = 0.2;
kappa1 = omega./U.f(y2);

x1 = -2:0.01:5;
x2 = 0:0.01:2;
% x2 = [0.4,0.5,0.6];

%% Scaling constant
A = 1; % to be fixed

%% Contour setup
N_C = 100;

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

%% Critical layer contours

N_CL = 1000; % More precision needed, :(
a = omega./U.f(delta);
b = omega./U.f(0);
shift1 = 0.4*omega;
shift2 = 0.4*omega;

[C_CLp,C_CLm] = Finite_CL_Ellipse(N_CL,a,b,shift1,shift2);

%% X-depdendence
Xp = Fourier_Nov2018(x1(x1<0),Cp.k);
Xm = Fourier_Nov2018(x1(x1>=0),Cm.k);
%X = [Xp;Xm];

X_CLp = Fourier_Nov2018(x1(x1>=0),C_CLp.k);
X_CLm = Fourier_Nov2018(x1(x1>=0),C_CLm.k);

%% Y-dependence
% Decidedly not quick
[Yp] = Y_GreensHard(x2,y2,Cp.k,k3,omega,U,c,delta,N_bl,0);
[Ym] = Y_GreensHard(x2,y2,Cm.k,k3,omega,U,c,delta,N_bl,0);

[Y_CLp] = Y_GreensHard(x2,y2,C_CLp.k,k3,omega,U,c,delta,N_bl,1);
[Y_CLm] = Y_GreensHard(x2,y2,C_CLm.k,k3,omega,U,c,delta,N_bl,-1);

%% k-dependence
Kp =  K_Hard(y2,Cp.k,k3,omega,U,c,delta,N_bl,0);
Km =  K_Hard(y2,Cm.k,k3,omega,U,c,delta,N_bl,0);

K_CLp = K_Hard(y2,C_CLp.k,k3,omega,U,c,delta,N_bl,1);
K_CLm = K_Hard(y2,C_CLm.k,k3,omega,U,c,delta,N_bl,-1);

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
