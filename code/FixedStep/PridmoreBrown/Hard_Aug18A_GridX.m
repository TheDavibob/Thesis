function [hard] = Hard_Aug18A_GridX(x2_max,k1,k3,omega,U,c,N)
% Check this solution is actually correct...


k1 = k1(:).';


% Hard-wall solution (within boundary layer)
phih0 = @(k1) IC_Hard_PB_Aug18A(k1,omega,U,c);
fh = @(x2,phi,k1) ODEFun_PB_Aug18A_GridX(x2,phi,k1,omega,c,U,k3) ;

[phi2,x2,phi2_end] = Fixed_Step_IVP_GridX(fh,phih0,x2_max,N,k1);
hard.x = x2;
hard.phi = permute(phi2(:,1,:),[1,3,2]) ;
hard.dphi = permute(phi2(:,2,:),[1 3 2]) ;

hard.phi = transpose(hard.phi);
hard.dphi = transpose(hard.dphi);

C = 1i*(omega - k1(:).*U.f(hard.x));
dC = -1i*k1(:).*U.df(hard.x);

hard.p = -C.^3.*hard.phi;
hard.v = - 3*dC.*C.*hard.phi - C.^2.*hard.dphi;

hard.phi_end = reshape(phi2_end(1,:),size(k1));
hard.dphi_end = reshape(phi2_end(2,:),size(k1));

C_end = 1i*(omega - k1(:).'.*U.f(x2_max));
dC_end = -1i*k1(:).'.*U.df(x2_max);

hard.p_end = reshape(-C_end.^3.*phi2_end(1,:),size(k1));
hard.v_end = reshape(-C_end.^2.*phi2_end(2,:) - 3*dC_end.*C_end.*phi2_end(1,:),size(k1));

end