function [dec] = Comp_Dec_Aug18A_skip(k1,k3,omega,U,c,delta,N,z)

% skip: evaluates only at s = 1

%N = 1000;

k1 = k1(:).';


% Decaying solution (within boundary layer)
phidinf = @(k1) Comp_IC_Decaying_PB_Aug18A(k1,k3,omega,U,c,delta);
fd = @(s,phi,k1) Comp_ODEFun_PB_Aug18A(s,phi,k1,omega,c,U,k3,z) ;

[phi1] = Fixed_Step_IVP_skip(fd,phidinf,1,N,k1);

%dec.s = x1;
%dec.z = z.f(dec.s);
dec.phi = phi1(1,:) ;
dec.dphi = phi1(2,:) ;

%dec.phi = transpose(dec.phi);
%dec.dphi = transpose(dec.dphi);

C = 1i*(omega - k1*U.f(z.f(1)));
dC = -1i*k1*U.df(z.f(1));

dec.p = -C.^3.*dec.phi;
dec.v = - 3*dC.*C.*dec.phi - C.^2.*dec.dphi;

end