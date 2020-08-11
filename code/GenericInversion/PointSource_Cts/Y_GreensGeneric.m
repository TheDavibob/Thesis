function [Y] = Y_GreensGeneric(x2,y2,k1,k3,omega,U,c,l1,l2,delta,N_bl,HP)
%Y_GREENSL Computes Y for the generic linear BC case

% Y is simply phi_><(x2)phi_<>(y2), with the careful jump stuff absorbed
% into the k bit.
% Done with a brute-loop in x2, as with the scattering case. This makes it
% slower than it could be, but needs only be computed "once".

Y.phi = zeros(numel(x2),numel(k1));
Y.dphi = zeros(numel(x2),numel(k1));
Y.p = zeros(numel(x2),numel(k1));
Y.v = zeros(numel(x2),numel(k1));

N_y = numel(x2);
N_ybl = numel(x2(x2<delta));

% Going to brute-loop the decaying solution within the boundary layer
%N_bl = 1000;
h1 = waitbar(0,'Computing x2 dependence');

if HP == 0
    zdy.f = @(s) delta*(1-s) + y2*s;
    zdy.df = @(s) (y2 - delta)*ones(size(s));
    zly.f = @(s) y2*s;
    zly.df = @(s) y2*ones(size(s));
else
    if HP > 0
        I = 1;
    else
        I = -1;
    end
    zdy.f = @(s) y2 + ((delta-y2)/2)*(1+ exp(I*1i*pi*s)) ;
    zdy.df = @(s) ((delta-y2)/2)*I*1i*pi*exp(I*1i*pi*s);
    zly.f = @(s) (y2/2)*(1-exp(-I*1i*pi*s));
    zly.df = @(s) (y2/2)*I*1i*pi*exp(-I*1i*pi*s);
end

phil_y2 = Comp_Generic_Aug18A_skip(k1,k3,omega,U,c,l1,l2,N_bl,zly);
Ell = phil_y2.phi;
dec_y2 = Comp_Dec_Aug18A_skip(k1,k3,omega,U,c,delta,N_bl,zdy) ;
d = dec_y2.phi;

for j = 1:N_ybl
    X2 = x2(j);
    if X2 > y2
        if HP == 0
            z.f = @(s) delta*(1-s) + X2*s; % so from delta to x2
            z.df = @(s) (X2 - delta)*ones(size(s));
        else
            z.f = @(s) X2 + ((delta-X2)/2)*(1+ exp(I*1i*pi*s)) ;
            z.df = @(s) ((delta-X2)/2)*I*1i*pi*exp(I*1i*pi*s);
        end
        [dec] = Comp_Dec_Aug18A_skip(k1,k3,omega,U,c,delta,N_bl,z) ;
        Y.phi(j,:) = Ell.*dec.phi;
        Y.dphi(j,:) = Ell.*dec.dphi;
        Y.p(j,:) = Ell.*dec.p;
        Y.v(j,:) = Ell.*dec.v;
    else
        if HP == 0
            z.f = @(s) X2*s;
            z.df = @(s) X2*ones(size(s));
        else
            z.f = @(s) (X2/2)*(1-exp(-I*1i*pi*s));
            z.df = @(s) (X2/2)*I*1i*pi*exp(-I*1i*pi*s);
        end
        [phil] = Comp_Generic_Aug18A_skip(k1,k3,omega,U,c,l1,l2,N_bl,z);
        Y.phi(j,:) = d.*phil.phi;
        Y.dphi(j,:) = d.*phil.dphi;
        Y.p(j,:) = d.*phil.p;
        Y.v(j,:) = d.*phil.v;
    end
    waitbar(j./N_ybl,h1,'Computing x2 dependence');
end
close(h1);

Y.phi(abs(Y.phi) < 5*eps) = 0;
Y.dphi(abs(Y.dphi) < 5*eps) = 0;
Y.p(abs(Y.p) < 5*eps) = 0;
Y.v(abs(Y.v) < 5*eps) = 0;

% Outwith the boundary-layer: "known" results
%disp('Computing auxiliary functions outside the boundary-layer');
gammainf = Gamma_FF(k3,omega,c.f(delta),U.f(delta));
G_inf = gammainf(k1);
C_inf = 1i*(omega - U.f(delta)*k1);
%dC_inf = -1i*U.df(delta)*k1;

Y.phi(x2 >= delta,:) = Ell.*exp(-x2(x2>=delta).'*G_inf);
Y.dphi(x2 >= delta,:) = -Ell.*G_inf.*exp(-x2(x2>=delta).'*G_inf);
Y.p(x2 >= delta,:) = -C_inf.^3.*Y.phi(x2 >= delta,:);
Y.v(x2 >= delta,:) = -C_inf.^2.*Y.dphi(x2 >= delta,:);


end

