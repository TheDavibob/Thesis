function [dec] = Dec_Aug18A_GridX_skip(x2_min,k1,k3,omega,U,c,delta,N)
K1 = k1(:).';
X2_min = x2_min(:).';

% Decaying solution (within boundary layer)
phidinf = @(k1) IC_Decaying_PB_Aug18A(k1,k3,omega,U,c,delta);
fd = @(x2,phi,k1) ODEFun_PB_Aug18A_GridX(delta-x2,phi,k1,omega,c,U,k3) ;

[phi] = Fixed_Step_IVP_GridX_skip(@(x2,phi,k1) -fd(x2,phi,k1),phidinf,delta - x2_min,N,K1);

C_end = 1i*(omega - K1.*U.f(X2_min));
dC_end = -1i*K1.*U.df(X2_min);

dec.phi = reshape(phi(1,:),size(k1));
dec.dphi = reshape(phi(2,:),size(k1));

dec.p = reshape(-C_end.^3.*phi(1,:),size(k1));
dec.v = reshape(-C_end.^2.*phi(2,:) - 3*dC_end.*C_end.*phi(1,:),size(k1));

end