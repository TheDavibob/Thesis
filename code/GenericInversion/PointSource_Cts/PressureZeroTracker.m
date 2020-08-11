%% A simple zero finding script, as a parameter changes.

c0 = 6:1:100;

delta = 1;
U = Struct_U_Linear(0.2,1,1);

omega = 1;
k3 = 0;
N_bl = 1000;

q0 = 1;

y2 = 0.5;
kappa1 = omega./U.f(y2);

smol = 1e-14;

% BCs: here don't depend on c, fortunately
C0 = @(k1) 1i*(omega - U.f(0)*k1);
dC0 = @(k1) -1i*U.df(0)*k1;

% Pressure release
l1 = @(k1) zeros(size(k1));
l2 = @(k1) -C0(k1).^3;

z = zeros(size(c0));

h = waitbar(0,'Here we go');

c = Struct_c_Const(c0(1));
D_l = @(k1) Dispersion_Generic_Nov2018(k1,k3,omega,U,c,l1,l2,delta,N_bl,1);
dD_l = @(k1) (D_l(k1 + smol) - D_l(k1 - smol)) ./(2*smol);

z(1) = NewtonRaphson(D_l,dD_l,1.4+0.3i);

for j = 2:numel(c0)
    c = Struct_c_Const(c0(j));
    D_l = @(k1) Dispersion_Generic_Nov2018(k1,k3,omega,U,c,l1,l2,delta,N_bl,1);
    dD_l = @(k1) (D_l(k1 + smol) - D_l(k1 - smol)) ./(2*smol);
    z(j) = NewtonRaphson(D_l,dD_l,z(j-1));
    waitbar(j./numel(c0),h,'Here we go');
end

close(h);