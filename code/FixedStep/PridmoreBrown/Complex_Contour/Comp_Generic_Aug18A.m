function [phil] = Comp_Generic_Aug18A(k1,k3,omega,U,c,l1,l2,N,z)

% integrates from s = 0 to s = 1 along contour z.

% l1, l2 anonymous functions of k1 (but not of omega, since assumed fixed).
% They should automatically be evaluated at x2 = 0.

k1 = k1(:).';

% Decaying solution (within boundary layer)
phil0 = @(k1) IC_Generic_PB_Aug18A(k1,l1,l2);
fl = @(s,phi,k1) Comp_ODEFun_PB_Aug18A(s,phi,k1,omega,c,U,k3,z) ;

[phi1,x1] = Fixed_Step_IVP(fl,phil0,1,N,k1);

phil.s = x1;
phil.z = z.f(phil.s);

phil.phi = permute(phi1(:,1,:),[1,3,2]) ;
phil.dphi = permute(phi1(:,2,:),[1 3 2]) ;

phil.phi = transpose(phil.phi);
phil.dphi = transpose(phil.dphi);

C = 1i*(omega - k1(:)*U.f(phil.z));
dC = -1i*k1(:)*U.df(phil.z);

phil.p = -C.^3.*phil.phi;
phil.v = - 3*dC.*C.*phil.phi - C.^2.*phil.dphi;

end