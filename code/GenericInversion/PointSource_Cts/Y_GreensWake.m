function [Y] = Y_GreensWake(x2,y2,k1,k3,omega,Up,cp,Um,cm,l1,l2,deltap,deltam,N_bl,HP)
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
N_ybl = numel(x2((x2<deltap) & (-x2<deltam)));

% Going to brute-loop the decaying solution within the boundary layer
%N_bl = 1000;
h1 = waitbar(0,'Computing x2 dependence');

if HP == 0
    zdy.f = @(s) deltap*(1-s) + y2*s;
    zdy.df = @(s) (y2 - deltap)*ones(size(s));
    zly.f = @(s) y2*s;
    zly.df = @(s) y2*ones(size(s));
else
    if HP > 0
        I = 1;
    else
        I = -1;
    end
    zdy.f = @(s) y2 + ((deltap-y2)/2)*(1+ exp(I*1i*pi*s)) ;
    zdy.df = @(s) ((deltap-y2)/2)*I*1i*pi*exp(I*1i*pi*s);
    zly.f = @(s) (y2/2)*(1-exp(-I*1i*pi*s));
    zly.df = @(s) (y2/2)*I*1i*pi*exp(-I*1i*pi*s);
end

phil_y2 = Comp_Generic_Aug18A_skip(k1,k3,omega,Up,cp,l1,l2,N_bl,zly);
Ell = phil_y2.phi;
dec_y2 = Comp_Dec_Aug18A_skip(k1,k3,omega,Up,cp,deltap,N_bl,zdy) ;
d = dec_y2.phi;

jmin = numel(x2(find(x2<-deltam)));
for j = jmin + 1 + (1:N_ybl)
    X2 = x2(j);
    if X2 > y2
        if HP == 0
            z.f = @(s) deltap*(1-s) + X2*s; % so from delta to x2
            z.df = @(s) (X2 - deltap)*ones(size(s));
        else
            z.f = @(s) X2 + ((deltap-X2)/2)*(1+ exp(I*1i*pi*s)) ;
            z.df = @(s) ((deltap-X2)/2)*I*1i*pi*exp(I*1i*pi*s);
        end
        [dec] = Comp_Dec_Aug18A_skip(k1,k3,omega,Up,cp,deltap,N_bl,z) ;
        Y.phi(j,:) = Ell.*dec.phi;
        Y.dphi(j,:) = Ell.*dec.dphi;
        Y.p(j,:) = Ell.*dec.p;
        Y.v(j,:) = Ell.*dec.v;
    elseif X2 >=0
        if HP == 0
            z.f = @(s) X2*s;
            z.df = @(s) X2*ones(size(s));
        else
            z.f = @(s) (X2/2)*(1-exp(-I*1i*pi*s));
            z.df = @(s) (X2/2)*I*1i*pi*exp(-I*1i*pi*s);
        end
        [phil] = Comp_Generic_Aug18A_skip(k1,k3,omega,Up,cp,l1,l2,N_bl,z);
        Y.phi(j,:) = d.*phil.phi;
        Y.dphi(j,:) = d.*phil.dphi;
        Y.p(j,:) = d.*phil.p;
        Y.v(j,:) = d.*phil.v;
    else
        X2 = -X2;
        if HP == 0
            z.f = @(s) deltam*(1-s) + X2*s; % so from delta to x2
            z.df = @(s) (X2 - deltam)*ones(size(s));
        else
            z.f = @(s) X2 + ((deltam-X2)/2)*(1+ exp(I*1i*pi*s)) ;
            z.df = @(s) ((deltam-X2)/2)*I*1i*pi*exp(I*1i*pi*s);
        end
        [dec] = Comp_Dec_Aug18A_skip(k1,k3,omega,Um,cm,deltam,N_bl,z) ;
        Y.phi(j,:) = d.*dec.phi;
        Y.dphi(j,:) = -d.*dec.dphi;
        Y.p(j,:) = d.*dec.p;
        Y.v(j,:) = -d.*dec.v;
    end
    waitbar((j-(1+jmin))./N_ybl,h1,'Computing x2 dependence');
end
close(h1);

%Y.phi(abs(Y.phi) < 5*eps) = 0;
%Y.dphi(abs(Y.dphi) < 5*eps) = 0;
%Y.p(abs(Y.p) < 5*eps) = 0;
%Y.v(abs(Y.v) < 5*eps) = 0;

% Outwith the boundary-layer: "known" results
%disp('Computing auxiliary functions outside the boundary-layer');
gammainf = Gamma_FF(k3,omega,cp.f(deltap),Up.f(deltap));
G_inf = gammainf(k1);
C_inf = 1i*(omega - Up.f(deltap)*k1);
%dC_inf = -1i*U.df(delta)*k1;

Y.phi(x2 >= deltap,:) = Ell.*exp(-x2(x2>=deltap).'*G_inf);
Y.dphi(x2 >= deltap,:) = -Ell.*G_inf.*exp(-x2(x2>=deltap).'*G_inf);
Y.p(x2 >= deltap,:) = -C_inf.^3.*Y.phi(x2 >= deltap,:);
Y.v(x2 >= deltap,:) = -C_inf.^2.*Y.dphi(x2 >= deltap,:);

gammainfm = Gamma_FF(k3,omega,cm.f(deltam),Um.f(deltam));
G_inf_m = gammainfm(k1);
C_inf_m = 1i*(omega - Um.f(deltam)*k1);
%dC_inf = -1i*U.df(delta)*k1;

Y.phi(-x2 >= deltam,:) = d.*exp(x2(-x2>=deltap).'*G_inf_m);
Y.dphi(-x2 >= deltam,:) = d.*G_inf_m.*exp(x2(-x2>=deltam).'*G_inf_m);
Y.p(-x2 >= deltam,:) = -C_inf_m.^3.*Y.phi(-x2 >= deltam,:);
Y.v(-x2 >= deltam,:) = -C_inf_m.^2.*Y.dphi(-x2 >= deltam,:);


end

