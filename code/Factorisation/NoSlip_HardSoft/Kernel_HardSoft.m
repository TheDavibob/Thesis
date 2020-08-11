function K = Kernel_HardSoft(k1,k3,omega,U,c,delta,N,z)

% Uses Comp_Dec_Aug18A_skip to compute the kernel for (slipping or not)
% hard to soft transition. Can easily be anonymised.

% z a suitably choen contour: an UHP one works swimmingly.

[dec] = Comp_Dec_Aug18A_skip(k1,k3,omega,U,c,delta,N,z);
K = dec.p ./dec.v ;

end
