function [hard] = Hard_Aug18A(k1,k3,omega,U,c,delta,N)
% Check this solution is actually correct...


k1 = k1(:).';


% Hard-wall solution (within boundary layer)
phih0 = @(k1) IC_Hard_PB_Aug18A(k1,omega,U,c);
fh = @(x2,phi,k1) ODEFun_PB_Aug18A(x2,phi,k1,omega,c,U,k3) ;

[phi2,x2] = Fixed_Step_IVP(fh,phih0,delta,N,k1);
hard.x = x2;
hard.phi = permute(phi2(:,1,:),[1,3,2]) ;
hard.dphi = permute(phi2(:,2,:),[1 3 2]) ;

hard.phi = transpose(hard.phi);
hard.dphi = transpose(hard.dphi);

C = 1i*(omega - k1(:)*U.f(hard.x));
dC = -1i*k1(:)*U.df(hard.x);

hard.p = -C.^3.*hard.phi;
hard.v = - 3*dC.*C.*hard.phi - C.^2.*hard.dphi;

end