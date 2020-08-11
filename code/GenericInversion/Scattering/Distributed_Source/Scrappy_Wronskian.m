y2 = 0.01:0.01:0.99;
W = zeros(size(y2));
A = zeros(size(y2));
B = zeros(size(y2));
C = zeros(size(y2));
D = zeros(size(y2));

for j = 1:numel(y2)
    shift_y2 = 0.001;
    if I == 0
        z_outer.f = @(s) delta*(1-s) + (y2(j)+shift_y2)*s; % so from delta to y
        z_outer.df = @(s) ((y2(j)+shift_y2) - delta)*ones(size(s));
        z_inner.f = @(s) (y2(j)-shift_y2)*s;
        z_inner.df = @(s) (y2(j)-shift_y2)*ones(size(s));
    else
        z_outer.f = @(s) y2(j) + ((delta-(y2(j)+shift_y2))/2)*(1+ exp(I*1i*pi*s)) ;
        z_outer.df = @(s) ((delta-(y2(j)+shift_y2))/2)*I*1i*pi*exp(I*1i*pi*s);
        z_inner.f = @(s) ((y2(j)-shift_y2)/2)*(1-exp(-I*1i*pi*s));
        z_inner.df = @(s) ((y2(j)-shift_y2)/2)*I*1i*pi*exp(-I*1i*pi*s);
    end
    kappa1 = omega./U.f(y2(j))+shift*1i;
    phi_dec = Comp_Dec_Aug18A_skip(kappa1,k3,omega,U,c,delta,N_bl,z_outer) ;
    phi_inner = Comp_Generic_Aug18A_skip(kappa1,k3,omega,U,c,l1,l2,N_bl,z_inner);
    
    A(j) = phi_dec.phi;
    B(j) = phi_dec.dphi;
    C(j) = phi_inner.phi;
    D(j) = phi_inner.dphi;
    W(j) = phi_dec.dphi.*phi_inner.phi - phi_dec.phi.*phi_inner.dphi;
end