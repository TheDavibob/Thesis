function [phi_0] = Upstream_Distributed_Integrated(x1,omega_c,G_01,dkappa1)
% Integrates G_01 over some vorticity distribution omega_c.

% x2 grid contained within G_01, as is y2, kappa1 grid (though dkappa1 has
% not yet been prescribed).

% G_01 is size x2 x y2

omega_c = omega_c(:).'; % size 1 x y2
dkappa1 = dkappa1(:).'; % size 1 x y2

y2 = G_01.y2;
kappa1 = G_01.kappa1;
x2 = G_01.x2;

x1 = x1(:).'; % size 1 x x1

omega_weighted = G_01.omega*omega_c.*G_01.c.f(y2).^2.*dkappa1./G_01.U.f(y2).^2;

fourier = exp(-1i*kappa1.'*x1); % size y2 x x1

I.phi = zeros(numel(x1),numel(x2),numel(y2));
I.dphi = zeros(numel(x1),numel(x2),numel(y2));
I.p = zeros(numel(x1),numel(x2),numel(y2));
I.v = zeros(numel(x1),numel(x2),numel(y2));


% The integration shall be looped. Being lazy with variables.
h = waitbar(0,'Integrating vorticity against Green''s function');
for j = 1:numel(y2)
    I.phi(:,:,j) = omega_weighted(j)*fourier(j,:).'*G_01.phi(:,j).' ;
    I.dphi(:,:,j) = omega_weighted(j)*fourier(j,:).'*G_01.dphi(:,j).' ;
    I.p(:,:,j) = omega_weighted(j)*fourier(j,:).'*G_01.p(:,j).' ;
    I.v(:,:,j) = omega_weighted(j)*fourier(j,:).'*G_01.v(:,j).' ;
    waitbar(j./numel(y2),h);
end
close(h);
phi_0.phi = sum(I.phi,3);
phi_0.dphi = sum(I.dphi,3);
phi_0.p = sum(I.p,3);
phi_0.v = sum(I.v,3);

end