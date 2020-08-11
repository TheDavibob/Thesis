function [Y] = Y_Continuous_Scattering_Nov2018(x2,k1,k3,omega,U,c,delta,N_bl)
%The x2 dependent part, i.e. just phid. This is easy!

Y.phi = zeros(numel(x2),numel(k1));
Y.dphi = zeros(numel(x2),numel(k1));
Y.p = zeros(numel(x2),numel(k1));
Y.v = zeros(numel(x2),numel(k1));

N_y = numel(x2);
N_ybl = numel(x2(x2<delta));

% Going to brute-loop the decaying solution within the boundary layer
%N_bl = 1000;
h1 = waitbar(0,'Computing x2 dependence');
for j = 1:N_ybl
    y = x2(j);
    z.f = @(s) delta*(1-s) + y*s; % so from delta to y
    z.df = @(s) (y - delta)*ones(size(s));
    [dec] = Comp_Dec_Aug18A_skip(k1,k3,omega,U,c,delta,N_bl,z) ;
    Y.phi(j,:) = dec.phi;
    Y.dphi(j,:) = dec.dphi;
    Y.p(j,:) = dec.p;
    Y.v(j,:) = dec.v;
    waitbar(j./N_ybl,h1,'Computing x2 dependence');
end
close(h1);

% Outwith the boundary-layer: "known" results
%disp('Computing auxiliary functions outside the boundary-layer');
gammainf = Gamma_FF(k3,omega,c.f(delta),U.f(delta));
G_inf = gammainf(k1);
C_inf = 1i*(omega - U.f(delta)*k1);
%dC_inf = -1i*U.df(delta)*k1;

Y.phi(x2 >= delta,:) = exp(-x2(x2>=delta).'*G_inf);
Y.dphi(x2 >= delta,:) = -G_inf.*exp(-x2(x2>=delta).'*G_inf);
Y.p(x2 >= delta,:) = -C_inf.^3.*Y.phi(x2 >= delta,:);
Y.v(x2 >= delta,:) = -C_inf.^2.*Y.dphi(x2 >= delta,:);


end

