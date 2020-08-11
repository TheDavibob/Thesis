function [dec] = Dec_Aug18A_GridX(x2_min,k1,k3,omega,U,c,delta,N)

%N = 1000;

k1 = k1(:).';


% Decaying solution (within boundary layer)
phidinf = @(k1) IC_Decaying_PB_Aug18A(k1,k3,omega,U,c,delta);
fd = @(x2,phi,k1) ODEFun_PB_Aug18A_GridX(delta-x2,phi,k1,omega,c,U,k3) ;

[phi1,x1,phi1_end] = Fixed_Step_IVP_GridX(@(x2,phi,k1) -fd(x2,phi,k1),phidinf,delta - x2_min,N,k1);

dec.x = delta - x1;
dec.phi = permute(phi1(end:-1:1,1,:),[1,3,2]) ;
dec.dphi = permute(phi1(end:-1:1,2,:),[1 3 2]) ;

dec.phi = transpose(dec.phi);
dec.dphi = transpose(dec.dphi);

C = 1i*(omega - k1(:).*U.f(dec.x));
dC = -1i*k1(:).*U.df(dec.x);

dec.p = -C.^3.*dec.phi;
dec.v = - 3*dC.*C.*dec.phi - C.^2.*dec.dphi;

dec.phi_end = reshape(phi1_end(1,:),size(k1));
dec.dphi_end = reshape(phi1_end(2,:),size(k1));

end