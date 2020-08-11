function I0 = Scrappy_I0(k1,y2,I)

% Given y2, finds I0(k1) where k1 isn't yet kappa1.

U0 = 0.3;
U = Struct_U_Linear(U0,1,1);
c = Struct_c_Const(5);
delta = 1;
N_bl = 6000;
omega = 1;
k3 = 0;

% I = -1;

C0 = @(k1) 1i*(omega - U.f(0)*k1);
dC0 = @(k1) -1i*U.df(0)*k1;
l1 = @(k1) (-C0(k1).^2);
l2 = @(k1) -3*C0(k1).*dC0(k1);

% y2 = 0.5;

%% Three bits

shift_y2 = 0;
if I == 0
        z_outer.f = @(s) delta*(1-s) + (y2+shift_y2)*s; % so from delta to y
        z_outer.df = @(s) ((y2+shift_y2) - delta)*ones(size(s));
else
        z_outer.f = @(s) y2 + ((delta-(y2+shift_y2))/2)*(1+ exp(I*1i*pi*s)) ;
        z_outer.df = @(s) ((delta-(y2+shift_y2))/2)*I*1i*pi*exp(I*1i*pi*s);
end
[dec_y2] = Comp_Dec_Aug18A_skip(k1,k3,omega,U,c,delta,N_bl,z_outer) ;

if I == 0
        z_0.f = @(s) delta*(1-s); % so from delta to y
        z_0.df = @(s) ( - delta)*ones(size(s));
else
        z_0.f = @(s) ((delta)/2)*(1+ exp(I*1i*pi*s)) ;
        z_0.df = @(s) ((delta)/2)*I*1i*pi*exp(I*1i*pi*s);
end
[dec_0] = Comp_Dec_Aug18A_skip(k1,k3,omega,U,c,delta,N_bl,z_0) ;

C_0 = 1i*(omega - U.f(0)*k1);
C_y2 = 1i*(omega - U.f(y2+shift_y2)*k1);

P_0 = -C_0.^3.*l1(k1);

Rat = C_y2.^3./C_0.^4;
I0 = Rat.*P_0.*dec_y2.phi./(l1(k1).*dec_0.dphi + l2(k1).*dec_0.phi);

end