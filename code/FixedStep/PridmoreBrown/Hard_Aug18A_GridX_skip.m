function [hard] = Hard_Aug18A_GridX_skip(x2_max,k1,k3,omega,U,c,~,N)

K1 = k1(:).';
X2_max = x2_max(:).';

k1 = k1(:).';


% Hard-wall solution (within boundary layer)
phih0 = @(k1) IC_Hard_PB_Aug18A(k1,omega,U,c);
fh = @(x2,phi,k1) ODEFun_PB_Aug18A_GridX(x2,phi,k1,omega,c,U,k3) ;

[phi] = Fixed_Step_IVP_GridX_skip(fh,phih0,x2_max,N,K1);
C_end = 1i*(omega - K1.*U.f(X2_max));
dC_end = -1i*K1.*U.df(X2_max);

hard.phi = reshape(phi(1,:),size(k1));
hard.dphi = reshape(phi(2,:),size(k1));

hard.p = reshape(-C_end.^3.*phi(1,:),size(k1));
hard.v = reshape(-C_end.^2.*phi(2,:) - 3*dC_end.*C_end.*phi(1,:),size(k1));
end