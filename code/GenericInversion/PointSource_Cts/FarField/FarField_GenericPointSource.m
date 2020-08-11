% Parameters as passed from the first part of PointSource_Generic

%% Parameters
delta = 1;
% U = Struct_U_Linear(0.3,1,delta);
U = Struct_U_Parabolic(0.3,1,delta);
c = Struct_c_Const(5);
omega = 20;
k3 = 0;
N_bl = 1000;

q0 = 1;

y2 = 0.01;
kappa1 = omega./U.f(y2);

%% Boundary conditions

C0 = @(k1) 1i*(omega - U.f(0)*k1);
dC0 = @(k1) -1i*U.df(0)*k1;

% Pressure release
l1 = @(k1) zeros(size(k1));
l2 = @(k1) -C0(k1).^3;

% Hard-wall
% l1 = @(k1) -C0(k1).^2;
% l2 = @(k1) -3*C0(k1).*dC0(k1);

% Impedance, iwZv' - Cp' = 0
%Z = exp(-1i*pi/4);
%if isinf(Z) ~= 1
%     l1 = @(k1) 1i*omega*Z*(-C0(k1).^2) - 0;
%     l2 = @(k1) 1i*omega*Z*(-3*C0(k1).*dC0(k1)) - C0(k1).*(-C0(k1).^3);
% else
%     l1 = @(k1) 1i*omega*(-C0(k1).^2);
%     l2 = @(k1) 1i*omega*(-3*C0(k1).*dC0(k1));
% end
%% Scaling constant
A = 1; % to be fixed

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


%% Functional dependence
% "x2"-dependence (i.e. phi_ell(y2), with any derivatives included here too
Y = Y_BoundaryGeneric(y2,Saddle,k3,omega,U,c,l1,l2,delta,N_bl,0);
K = K_Generic(y2,Saddle,k3,omega,U,c,l1,l2,delta,N_bl,0);

FF.p = K.*sin(Theta).*Y.p;
FF.v = K.*sin(Theta).*Y.v;
FF.phi = K.*sin(Theta).*Y.phi;

polarplot(theta,abs(FF.p.'));


%% Example code: Multiple running
% comment out c
% theta = linspace(0,pi,1000);
% c0 = [2,5,10,20,50,100];
% c = Struct_c_Const(c0(1));
% F = zeros(numel(c0),numel(theta));
% G = zeros(numel(c0),numel(theta));
% for j = 1:numel(c0)
%     c = Struct_c_Const(c0(j));
%     FarField_GenericPointSource;
%     F(j,:) = FF.p;
%     G(j,:) = FF.v;
% end

% or: (comment out omega)

% Omega = 10.^(-2:1:3);
% theta = linspace(0,pi,1000);
% F = zeros(numel(Omega),numel(theta));
% G = zeros(numel(Omega),numel(theta));
% for j = 1:numel(Omega)
% omega = Omega(j);
% FarField_GenericPointSource;
% F(j,:) = FF.p;
% G(j,:) = FF.v;
% end
% polarplot(theta,abs(G)./abs(G(:,500)))

% Y2 = 0.1:0.1:0.9;
% F = zeros(numel(Y2),numel(theta));
% G = zeros(numel(Y2),numel(theta));
% for j = 1:numel(Y2)
% y2 = Y2(j);
% FarField_GenericPointSource;
% F(j,:) = FF.p;
% G(j,:) = FF.v;
% end

% Some plotting bollox
% clf
% N = numel(G(:,1));
% h(N) = polar(theta,abs(F(N,:))./max(abs(F(N,:))));
% h(N).LineWidth = 2;
% for j = N-1:-1:1
% hold on
% h(j) = polar(theta,abs(F(j,:))./max(abs(F(j,:))));
% hold off
% h(j).LineWidth = 2;
% end
% set(gca,'YLim',[0,1.7250])
