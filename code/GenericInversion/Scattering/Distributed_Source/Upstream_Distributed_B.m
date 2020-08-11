function [G_01] = Upstream_Distributed_B(x2,G_01_A,omega,k3,U,c,delta)
% Computes the Greens function for a range of y2.

y2 = G_01_A.y2;

G_01.phi = zeros(numel(x2),numel(y2));
G_01.dphi = zeros(numel(x2),numel(y2));
G_01.p = zeros(numel(x2),numel(y2));
G_01.v = zeros(numel(x2),numel(y2));

Gamma = Gamma_FF(k3,omega,c.f(delta),U.f(delta));
gamma = Gamma(G_01_A.kappa1);

h = waitbar(0,'Determining upstream Green''s function');
for j = 1:numel(y2)
    x2_a = x2(x2<y2(j));
    x2_y2 = x2(x2 == y2(j)); % a "tricky" case
    x2_b = x2((x2>y2(j)) & (x2<delta));
    x2_c = x2(x2>=delta);
    
    % Interpolate in each region
    phi_a = interp1(G_01_A.generic.z,G_01_A.generic.phi(j,:),x2_a,'pchip');
    phi_b = interp1(G_01_A.dec.z,G_01_A.dec.phi(j,:),x2_b,'pchip');
    dphi_a = interp1(G_01_A.generic.z,G_01_A.generic.dphi(j,:),x2_a,'pchip');
    dphi_b = interp1(G_01_A.dec.z,G_01_A.dec.dphi(j,:),x2_b,'pchip');
    
    p_a = interp1(G_01_A.generic.z,G_01_A.generic.p(j,:),x2_a,'pchip');
    p_b = interp1(G_01_A.dec.z,G_01_A.dec.p(j,:),x2_b,'pchip');
    
    v_a = interp1(G_01_A.generic.z,G_01_A.generic.v(j,:),x2_a,'pchip');
    v_b = interp1(G_01_A.dec.z,G_01_A.dec.v(j,:),x2_b,'pchip');
    
    G_01_a_phi = G_01_A.Q0(j).*G_01_A.dec_y2.p(j).*phi_a;
    G_01_a_dphi = G_01_A.Q0(j).*G_01_A.dec_y2.p(j).*dphi_a;
    
    G_01_a_p = G_01_A.Q0(j).*G_01_A.dec_y2.p(j).*p_a;
    G_01_a_v = G_01_A.Q0(j).*G_01_A.dec_y2.p(j).*v_a;
    
    G_01_b_phi = G_01_A.Q0(j).*G_01_A.generic_y2.p(j).*phi_b;
    G_01_b_dphi = G_01_A.Q0(j).*G_01_A.generic_y2.p(j).*dphi_b;
    
    G_01_b_p = G_01_A.Q0(j).*G_01_A.generic_y2.p(j).*p_b;
    G_01_b_v = G_01_A.Q0(j).*G_01_A.generic_y2.p(j).*v_b;
    
    if numel(x2_y2 == 1)
        G_01_y2_phi = G_01_A.Q0(j).*G_01_A.generic_y2.p(j).*G_01_A.dec_y2.phi(j);
        G_01_y2_dphi = 0; % Ill-defined, a pain.
        G_01_y2_p = G_01_A.Q0(j).*G_01_A.generic_y2.p(j).*G_01_A.dec_y2.p(j);
        G_01_y2_v = G_01_A.Q0(j).*G_01_A.generic_y2.p(j).*G_01_A.dec_y2.v(j);
    else
        G_01_y2_phi = [];
        G_01_y2_dphi = [];
        G_01_y2_p = [];
        G_01_y2_v = [];
        
    end
        
    G_01_c_phi = G_01_A.Q0_m(j).*G_01_A.generic_y2.phi(j).*exp(-gamma(j)*x2_c);
    G_01_c_dphi = -gamma(j)*G_01_A.Q0_m(j).*G_01_A.generic_y2.phi(j).*exp(-gamma(j)*x2_c);
    
    G_01_c_p = 1i*(omega - U.f(delta).*G_01_A.kappa1(j)).^3.*G_01_c_phi;
    G_01_c_v = (omega - U.f(delta).*G_01_A.kappa1(j)).^2.*G_01_c_dphi;
        
    G_01.phi(:,j) = [G_01_a_phi,G_01_y2_phi,G_01_b_phi,G_01_c_phi];
    G_01.dphi(:,j) = [G_01_a_dphi,G_01_y2_dphi,G_01_b_dphi,G_01_c_dphi];
    G_01.p(:,j) = [G_01_a_p,G_01_y2_p,G_01_b_p,G_01_c_p];
    G_01.v(:,j) = [G_01_a_v,G_01_y2_v,G_01_b_v,G_01_c_v];
    waitbar(j./numel(y2),h);
end
close(h);

C_inner = 1i*(omega - U.f(x2(x2<delta)).'*G_01_A.kappa1);
C_outer = 1i*(omega - (ones(size(x2(x2>=delta)))*U.f(delta)).'*G_01_A.kappa1);
C = [C_inner;C_outer];

dC_inner = -1i*U.df(x2(x2<delta)).'*G_01_A.kappa1;
dC_outer = zeros(numel(x2(x2>=delta)),numel(G_01_A.kappa1));
dC = [dC_inner;dC_outer];

% G_01.p = -C.^3.*G_01.phi;
% G_01.v = -3*C.*dC.*G_01.phi - C.^2.*G_01.dphi;

% An important change: since continuity was forced via pressure, which is
% O(1), we're going to flip the definition to define phi. This ensures the
% singularity is well understood

G_01.phi = -G_01.p./(C.^3);
G_01.dphi = -(G_01.v + 3*C.*dC.*G_01.phi)./(C.^2);

G_01.x2 = x2;
G_01.y2 = G_01_A.y2;
G_01.kappa1 = G_01_A.kappa1;

G_01.omega = G_01_A.omega;
G_01.c = G_01_A.c;
G_01.U = G_01_A.U;

end