function K = Kernel_Generic(k1,k3,omega,U,c,delta,Lu,Ld,N,z)

% Uses Comp_Dec_Aug18A_skip to compute the kernel for transition from Lu
% phi = 0 to Ld phi = 0.

% Lj structures, with components l1 and l2, so that L phi = l1 phi' + l2 phi = 0

% z a suitably choen contour: an UHP one works swimmingly, which should
% finish at 0.

[dec] = Comp_Dec_Aug18A_skip(k1,k3,omega,U,c,delta,N,z);

% Dispersion relations
Du = Lu.l1(k1).*dec.dphi + Lu.l2(k1).*dec.phi;
Dd = Ld.l1(k1).*dec.dphi + Ld.l2(k1).*dec.phi;

K = Dd./Du;

end

