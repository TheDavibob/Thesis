function [I_FF] = Oct_18_FF_SD(theta,y2,k3,omega,U,c,delta)
%OCT_18_FF_SD (Hard-wall) Far-field sound due to source above hard wall
%  	The usual, first order contribution only.

% theta is the real coordinate
x1 = cos(theta);
x2 = sin(theta);

N_SD = 1; % Dummy, don't care.
[k1] = Steepest_Descent_Contour(x1,x2,omega,U.f(delta),c.f(delta),k3,N_SD);
Theta = k1.Theta; % Mach modified theta.
Theta = diag(Theta).';
K0 = k1.K0; % Modified k0 = omega/c0

M = k1.M;
k0 = k1.k0;
K1 = K0*cos(Theta) - M*k0/(1-M^2);

% Constructing the integrand
N_BL = 1000; % consitency
% Vortex sheet term
    zh.f = @(t) t*y2;
    zh.df = @(t) y2;
    [hard_y2] = Comp_Hard_Aug18A_skip(K1,k3,omega,U,c,N_BL,zh);
    disp('Vortex sheet contribution computed')
% Dispersion relationship and Wronskian
    z.f = @(t) delta*(1-t);
    z.df = @(t) -delta*t;
    [dec_0] = Comp_Dec_Aug18A_skip(K1,k3,omega,U,c,delta,N_BL,z);
    Disp_hard = dec_0.v;
    Wronskian = Disp_hard./(c.f(0).^2.*(1i*(omega - U.f(0)*K1)).^5);
    disp('Hard-wall contribution computed')
% Cauchy term
    kappa1 = omega./U.f(y2);
    Cauchy = 1./(2*pi*1i*(K1- kappa1));
% Convected ratio
    C_y2 = omega - U.f(y2)*K1;
    C_0 = omega - U.f(0)*K1;
    Conv_rat = C_y2.^4./C_0.^4;
    
% And, together (without sqrt(pi/(RK0) )
I_FF.phi = Cauchy.*hard_y2.phi./Wronskian.*Conv_rat;
I_FF.p = Cauchy.*hard_y2.p./Wronskian.*Conv_rat;
I_FF.v = Cauchy.*hard_y2.v./Wronskian.*Conv_rat;
I_FF.dphi = Cauchy.*hard_y2.dphi./Wronskian.*Conv_rat;

I_FF.Theta = Theta;
I_FF.Disp_hard = Disp_hard;
I_FF.Cauchy = Cauchy;
I_FF.hard_y2 = hard_y2;
I_FF.Wronskian = Wronskian;
    disp('Solutions combined')

end

