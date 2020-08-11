function [f] = ODEFun_PB_Aug18A(x2,Phi_jm1,k1,omega,c,U,k3)
%ODEFUN_PB_AUG18 The ode function for the PB equation
%   Call as @(x2,Phi_jm1,k1) ODEFun...(...,params)
%   U and c are structures containing f and its derivatives.
%   k3, omega are fixed, singleton parameters as, implicitly, is x2



%% 
c0 = c.f(x2);
dc0 = c.df(x2);
u=U.f(x2);
du=U.df(x2);

C = 1i*(omega - u*k1);
dC = -1i*du*k1;
ddCc2 = -1i*U.d2f(x2).*c0.^2*k1 -2*dC.*c0.*dc0;

p = 2*dc0./c0 + 4*dC./C;
q = 3*ddCc2./(c0.^2.*C) - (k1.^2 + k3.^2) - C.^2/c0.^2;

f=zeros(2,numel(k1)); % with Phi having two components
f(1,:) = Phi_jm1(2,:) ;
f(2,:) = -q.*Phi_jm1(1,:)-p.*Phi_jm1(2,:);


end

