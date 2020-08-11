function [G_01_A] = Upstream_Distributed_A(y2,kappa1,shift_y2,omega,k3,U,c,delta,N_bl,l1,l2,HP)
% Careful computation of the Greens function for an upstream disturbance.
% Since going to want to integrate over some weighting, an absolute pain.

% Various warnings, for example the singularity at y2 causing issues.

% Note: this does not yet depend on x2. x2 is taken care of in _B, which
% interpolates on what we have here.

if HP > 0
    I = 1;
elseif HP < 0
    I = -1;
else
    I = 0;
end

%% Dispersion relationship
if I == 0
        z_0.f = @(s) delta*(1-s); % so from delta to 0
        z_0.df = @(s) ( - delta)*ones(size(s));
else
        z_0.f = @(s) ((delta)/2)*(1+ exp(I*1i*pi*s)) ;
        z_0.df = @(s) ((delta)/2)*I*1i*pi*exp(I*1i*pi*s);
end
[dec_0] = Comp_Dec_Aug18A_skip(kappa1,k3,omega,U,c,delta,N_bl,z_0) ;
DZ = (l1(kappa1).*dec_0.dphi + l2(kappa1).*dec_0.phi);

%% Convective quantities, small in shift_y2
C_y2_p = 1i*(omega - U.f(y2+shift_y2).*kappa1);
C_y2_m = 1i*(omega - U.f(y2-shift_y2).*kappa1);
C_0 = 1i*(omega - U.f(0).*kappa1);
Q0_p = -C_y2_p.^3./(c.f(0).^2.*C_0.^4.*DZ);
Q0_m = -C_y2_m.^3./(c.f(0).^2.*C_0.^4.*DZ);
Q0 = 1./(c.f(0).^2.*C_0.^4.*DZ); % Removal of the pressure term should clean this up a bit.

dec_y2_phi = zeros(size(y2));
generic_y2_phi = zeros(size(y2));
dec_y2_p = zeros(size(y2));
generic_y2_p = zeros(size(y2));

h = waitbar(0,'Computing functions at y2');
for j = 1:numel(y2)
%     if I == 0
        z_outer_lim.f = @(s) delta*(1-s) + (y2(j)+shift_y2)*s; % so from delta to y
        z_outer_lim.df = @(s) ((y2(j)+shift_y2) - delta)*ones(size(s));
        z_inner_lim.f = @(s) (y2(j)-shift_y2)*s;
        z_inner_lim.df = @(s) (y2(j)-shift_y2)*ones(size(s));
%     else
%         z_outer_lim.f = @(s) y2(j) + ((delta-(y2(j)+shift_y2))/2)*(1+ exp(I*1i*pi*s)) ;
%         z_outer_lim.df = @(s) ((delta-(y2(j)+shift_y2))/2)*I*1i*pi*exp(I*1i*pi*s);
%         z_inner_lim.f = @(s) (y2(j)-shift_y2)*(1-exp(-I*1i*pi*s))/2;
%         z_inner_lim.df = @(s) 1i*pi*I*(y2(j)-shift_y2)*exp(-I*1i*pi*s)/2;
%     end
    [dec_y2_temp] = Comp_Dec_Aug18A_skip(kappa1(j),k3,omega,U,c,delta,N_bl,z_outer_lim) ;
    [generic_y2_temp] = Comp_Generic_Aug18A_skip(kappa1(j),k3,omega,U,c,l1,l2,N_bl,z_inner_lim) ;
    dec_y2_phi(j) = dec_y2_temp.phi;
    dec_y2_p(j) = dec_y2_temp.p;
    generic_y2_phi(j) = generic_y2_temp.phi;
    generic_y2_p(j) = generic_y2_temp.p;
    waitbar(j./numel(y2),h);
end
close(h)

dec_y2.phi = dec_y2_phi;
generic_y2.phi = generic_y2_phi;

dec_y2.p = dec_y2_p;
generic_y2.p = generic_y2_p;

%% Above source
z_outer.f = @(s) delta*(1-s);
z_outer.df = @(s) (- delta)*ones(size(s));
[dec] = Comp_Dec_Aug18A(kappa1,k3,omega,U,c,delta,N_bl,z_outer);
% Interpolate from dec.z and dec.variable for y2<x2<delta

%% Below source
z_inner.f = @(s) delta*s;
z_inner.df = @(s) delta*ones(size(s));
[generic] = Comp_Generic_Aug18A(kappa1,k3,omega,U,c,l1,l2,N_bl,z_inner);
% Interpolate from generic.z and generic.variable for 0<x2<y2

%% Outputs (for insertion into _B)
G_01_A.dec = dec;
G_01_A.generic = generic;
G_01_A.DZ = DZ;
G_01_A.dec_y2 = dec_y2;
G_01_A.generic_y2 = generic_y2;
G_01_A.Q0_p = Q0_p;
G_01_A.Q0_m = Q0_m;
G_01_A.Q0 = Q0;

G_01_A.y2 = y2;
G_01_A.kappa1 = kappa1;

G_01_A.omega = omega;
G_01_A.U = U;
G_01_A.c = c;


end