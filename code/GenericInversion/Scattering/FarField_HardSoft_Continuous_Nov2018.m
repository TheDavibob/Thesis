%% The far-field noise for the hard-soft scattering problem, slip or not

%% Parameter setup
delta = 1;
%U = Struct_U_Linear(0,1,delta);
U = Struct_U_Parabolic(0,1,delta);
c = Struct_c_Const(5);
omega = 1;
k3 = 0;
N_bl = 2000;

y2 = 0.5;
kappa1 = omega./U.f(y2);

%% Far-field bit
N_theta = 1000;
theta = linspace(0,pi,N_theta);

x1 = cos(theta);
x2 = sin(theta);

[k1] = Steepest_Descent_Contour(x1,x2,omega,U.f(delta),c.f(delta),k3,1);
Theta = diag(k1.Theta).';
Saddle = k1.K0.*cos(Theta) - (k1.M).*k1.k0/(1-k1.M^2);
% The above is the steepest descent point as a function of Mach-modified
% angle, Theta.



%% Two alternative factorisations
z_kernel = Contour_Semicircle(delta,1); % So branch cut in LHP
prec.N = 5000;

shift = 0.5*(0.5.*(Saddle(1) - Saddle(end)));

[contour.Cp,contour.dCp] = TanhContour(0,3,-1+shift,1+shift);
[contour.Cm,contour.dCm] = TanhContour(0,3,-1-shift,1-shift);

Kernel = @(k1) Kernel_HardSoft(k1,k3,omega,U,c,delta,N_bl,z_kernel);
if U.f(0) ~= 0
    Kernel_ff = FFKernel_HardSoft_Slip(omega,U.f(0),c.f(0),k3,0);
else
    Kernel_ff = FFKernel_HardSoft_NoSlip(omega,U.df(0),c.f(0),k3,0);
end
branch_cut = pi;
[K_factorised] = Fact_Mult_Cont_Scaled_Finite(Kernel,Kernel_ff,contour,prec,branch_cut);
% has both + and - parts which work properly on a strip.

% If upstream scaling wanted
%[phid_y2] = Y_Continuous_Scattering_CL_Nov2018_Mod(y2,kappa1,k3,omega,U,c,delta,N_bl,1);
%Dh = phid_y2.v;
%q0 = 1;
%Q0 = q0.*c.f(y2).^2./(U.f(y2).*Dh);

Q0 = 1e-10;
Wall_Upstream = Y_GreensHard(0,y2,kappa1,k3,omega,U,c,delta,N_bl,-1);
Pressure_Upstream = Q0*Wall_Upstream.p; % i.e. phid(y2)*ph(0): easy enough to genericise
Km_kappa1 = K_factorised.m(kappa1);
A = -1i*Pressure_Upstream./Km_kappa1;


%% Functional dependence
pole = 1./(Saddle - kappa1);

z.f = @(s) delta*(1-s);
z.df = @(s) -delta*ones(size(s));
dec = Comp_Dec_Aug18A_skip(Saddle,k3,omega,U,c,delta,N_bl,z);
Dh = dec.v;
Dp = dec.p;

Kp = A*pole./(K_factorised.p(Saddle).*Dh);
Km = A*pole.*K_factorised.m(Saddle)./Dp;
% These two should precisely agree for saddle points between the two
% contours. Km is valid downstream and Kp is valid upstream.

% Only because it contains the derivatives
[Y] = Y_Continuous_Scattering_Nov2018_Mod(delta,Saddle,k3,omega,U,c,delta,10);
FFp.p = Kp.*sin(Theta).*Y.p;
FFp.v = Kp.*sin(Theta).*Y.v;
FFp.phi = Kp.*sin(Theta).*Y.phi;
FFm.p = Km.*sin(Theta).*Y.p;
FFm.v = Km.*sin(Theta).*Y.v;
FFm.phi = Km.*sin(Theta).*Y.phi;

FF = FFp;
FF.p(theta < pi/2) = FFm.p(theta<pi/2);
FF.v(theta < pi/2) = FFm.v(theta<pi/2);
FF.phi(theta < pi/2) = FFm.phi(theta<pi/2);
polarplot(theta,abs(FF.p.'));
