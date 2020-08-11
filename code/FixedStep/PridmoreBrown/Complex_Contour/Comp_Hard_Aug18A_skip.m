function [hard] = Comp_Hard_Aug18A_skip(k1,k3,omega,U,c,N,z)

% skip: evaluates only at s = 1

%N = 1000;

k1 = k1(:).';


% Decaying solution (within boundary layer)
phih0 = @(k1) IC_Hard_PB_Aug18A(k1,omega,U,c);
fh = @(s,phi,k1) Comp_ODEFun_PB_Aug18A(s,phi,k1,omega,c,U,k3,z) ;

[phi1] = Fixed_Step_IVP_skip(fh,phih0,1,N,k1);

%dec.s = x1;
%dec.z = z.f(dec.s);
hard.phi = phi1(1,:) ;
hard.dphi = phi1(2,:) ;

%dec.phi = transpose(dec.phi);
%dec.dphi = transpose(dec.dphi);

C = 1i*(omega - k1*U.f(z.f(1)));
dC = -1i*k1*U.df(z.f(1));

hard.p = -C.^3.*hard.phi;
hard.v = - 3*dC.*C.*hard.phi - C.^2.*hard.dphi;

end