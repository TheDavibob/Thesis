%% Computes the far-field noise due to a point source above a wake, possibly asymmetric

%% Parameters
deltap = 1;
c0 = 5;
U0 = 1;
Uinfp = 1;
% U = Struct_U_Linear(0.3,1,delta);
Up = Struct_U_Parabolic(U0,Uinfp,deltap);
cp = Struct_c_Const(c0);

deltam = 1;
Uinfm = 1;
Um = Struct_U_Parabolic(U0,Uinfm,deltam);
cm = Struct_c_Const(c0);

omega = 10;
k3 = 0;
N_bl = 1000;

q0 = 1;

y2 = 0.5;
kappa1 = omega./Up.f(y2);


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

%% Far-field bit
N_theta_p = 200;
N_theta_m = 200;
theta_p = linspace(0,pi,N_theta_p);
theta_m = linspace(pi,2*pi,N_theta_m);

x1_p = cos(theta_p);
x2_p = sin(theta_p);

x1_m = cos(theta_m);
x2_m = sin(theta_m);

[k1_p] = Steepest_Descent_Contour(x1_p,x2_p,omega,Up.f(deltap),cp.f(deltap),k3,1);
[k1_m] = Steepest_Descent_Contour(x1_m,-x2_m,omega,Um.f(deltam),cm.f(deltam),k3,1);

Theta_p = diag(k1_p.Theta).';
Saddle_p = k1_p.K0.*cos(Theta_p) - (k1_p.M).*k1_p.k0/(1-k1_p.M^2);

Theta_m = diag(k1_m.Theta).';
Saddle_m = k1_m.K0.*cos(Theta_m) - (k1_m.M).*k1_m.k0/(1-k1_m.M^2);
% The above is the steepest descent point as a function of Mach-modified
% angle, Theta.


%% Functional dependence
% "x2"-dependence (i.e. phi_ell(y2), with any derivatives included here too
Y_p = Y_BoundaryGeneric(y2,Saddle_p,k3,omega,Up,cp,l1,l2,deltap,N_bl,0);
K_p = K_Generic(y2,Saddle_p,k3,omega,Up,cp,l1,l2,deltap,N_bl,0);

Y_m = Y_BoundaryGeneric(y2,Saddle_m,k3,omega,Up,cp,l1,l2,deltap,N_bl,0);
K_m = K_Generic(y2,Saddle_m,k3,omega,Up,cp,l1,l2,deltap,N_bl,0);

FF.p = [K_p.*sin(Theta_p).*Y_p.p,K_m.*sin(Theta_m).*Y_m.p];
FF.v = [K_p.*sin(Theta_p).*Y_p.v,K_m.*sin(Theta_m).*Y_m.v];
FF.phi = [K_p.*sin(Theta_p).*Y_p.phi,K_m.*sin(Theta_m).*Y_m.phi];

theta = [theta_p,theta_m];

polarplot(theta,abs(FF.p.'));


%% Plotting bits:
% Deltam = 0.5:0.5:2;
% F = zeros(numel(Deltam),numel(theta));
% G = zeros(numel(Deltam),numel(theta));
% for j = 1:numel(Deltam)
% deltam = Deltam(j);
% FarField_WakePointSource_Feb19;
% F(j,:) = FF.p;
% G(j,:) = FF.v;
% end

% clf
% N = numel(G(:,1));
% clear h;
% h(1) = polar(theta,abs(F(1,:))./max(abs(F(1,:))));
% h(1).LineWidth = 2;
% for j = 2:N
% hold on
% h(j) = polar(theta,abs(F(j,:))./max(abs(F(j,:))));
% hold off
% h(j).LineWidth = 2;
% end
% set(gca,'YLim',[-1.7250,1.7250])