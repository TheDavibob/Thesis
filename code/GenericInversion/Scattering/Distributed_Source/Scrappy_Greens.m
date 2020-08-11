U0 = 0.3;
U = Struct_U_Parabolic(U0,1,1);
c = Struct_c_Const(5);
delta = 1;
N_bl = 2000;
omega = 1;
k3 = 0;

I = 0;

C0 = @(k1) 1i*(omega - U.f(0)*k1);
dC0 = @(k1) -1i*U.df(0)*k1;
l1 = @(k1) (-C0(k1).^2);
l2 = @(k1) -3*C0(k1).*dC0(k1);

y2 = 0.5;
kappa1 = omega./U.f(y2);

x2 = 0.0:0.01:1;
W = zeros(size(x2));
A = zeros(size(x2));
B = zeros(size(x2));
C = zeros(size(x2));
D = zeros(size(x2));

for j = 1:numel(x2)
    shift_y2 = 0.00;
    if I == 0
        z_outer.f = @(s) delta*(1-s) + (x2(j)+shift_y2)*s; % so from delta to y
        z_outer.df = @(s) ((x2(j)+shift_y2) - delta)*ones(size(s));
        z_inner.f = @(s) (x2(j)-shift_y2)*s;
        z_inner.df = @(s) (x2(j)-shift_y2)*ones(size(s));
    else
        z_outer.f = @(s) x2(j) + ((delta-(x2(j)+shift_y2))/2)*(1+ exp(I*1i*pi*s)) ;
        z_outer.df = @(s) ((delta-(x2(j)+shift_y2))/2)*I*1i*pi*exp(I*1i*pi*s);
        z_inner.f = @(s) ((x2(j)-shift_y2)/2)*(1-exp(-I*1i*pi*s));
        z_inner.df = @(s) ((x2(j)-shift_y2)/2)*I*1i*pi*exp(-I*1i*pi*s);
    end
    phi_dec = Comp_Dec_Aug18A_skip(kappa1,k3,omega,U,c,delta,N_bl,z_outer) ;
    phi_inner = Comp_Generic_Aug18A_skip(kappa1,k3,omega,U,c,l1,l2,N_bl,z_inner);
    
    A(j) = phi_dec.phi;
    B(j) = phi_dec.dphi;
    C(j) = phi_inner.phi;
    D(j) = phi_inner.dphi;
    W(j) = phi_dec.dphi.*phi_inner.phi - phi_dec.phi.*phi_inner.dphi;
end