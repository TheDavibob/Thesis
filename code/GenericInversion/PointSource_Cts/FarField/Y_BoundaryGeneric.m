function [Y] = Y_BoundaryGeneric(y2,k1,k3,omega,U,c,l1,l2,delta,N_bl,HP)
%Y_BoundaryGeneric Computes phi_ell, evaluated at y2 for generic bc.

if y2 > delta
    error('This code is calibrated only for sources within the boundary layer');
end

Y.phi = zeros(1,numel(k1));
Y.dphi = zeros(1,numel(k1));
Y.p = zeros(1,numel(k1));
Y.v = zeros(1,numel(k1));

%N_y = numel(x2);
%N_ybl = numel(x2(x2<delta));

% Going to brute-loop the decaying solution within the boundary layer
%N_bl = 1000;
%h1 = waitbar(0,'Computing x2 dependence');

if HP == 0
    zly.f = @(s) y2*s;
    zly.df = @(s) y2*ones(size(s));
else
    if HP > 0
        I = 1;
    else
        I = -1;
    end
    zly.f = @(s) (y2/2)*(1-exp(-I*1i*pi*s));
    zly.df = @(s) (y2/2)*I*1i*pi*exp(-I*1i*pi*s);
end

phil_y2 = Comp_Generic_Aug18A_skip(k1,k3,omega,U,c,l1,l2,N_bl,zly);
Ell = phil_y2.phi;
%dec_y2 = Comp_Dec_Aug18A_skip(k1,k3,omega,U,c,delta,N_bl,zdy) ;
%d = dec_y2.phi;

% Outwith the boundary-layer: "known" results
%disp('Computing auxiliary functions outside the boundary-layer');
gammainf = Gamma_FF(k3,omega,c.f(delta),U.f(delta));
G_inf = gammainf(k1);
C_inf = 1i*(omega - U.f(delta)*k1);
%dC_inf = -1i*U.df(delta)*k1;

Y.phi = Ell;
Y.dphi = -Ell.*G_inf;
Y.p = -C_inf.^3.*Y.phi;
Y.v = -C_inf.^2.*Y.dphi;


end

