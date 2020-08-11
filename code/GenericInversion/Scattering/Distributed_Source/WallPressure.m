function [I0] = WallPressure(y2,omega,k3,U,c,delta,l1,l2,N_bl,HP)
% Given background profile and boundary condition, computes the pressure on
% the wall for a range of y2, possibly with some deformation required

I0 = zeros(size(y2));
kappa1 = omega./U.f(y2); % a vector of kappa1s

Q0 = (c.f(y2).^2.*(omega - U.f(y2).*kappa1).^4)./(c.f(0).^2.*(omega - U.f(0).*kappa1).^4);
if HP == 0
    [phid] = Y_Continuous_Scattering_Nov2018_Mod([0,0],kappa1,k3,omega,U,c,delta,N_bl);
elseif HP > 0
    I = 1;
    [phid] = Y_Continuous_Scattering_CL_Nov2018_Mod([0,0],kappa1,k3,omega,U,c,delta,N_bl,I);
else
    I = -1;
    [phid] = Y_Continuous_Scattering_CL_Nov2018_Mod([0,0],kappa1,k3,omega,U,c,delta,N_bl,I);
end
D_l = l1(kappa1).*phid.dphi(1,:) + l2(kappa1).*phid.phi(1,:);

for j = 1:numel(I0)
    Y = Y_GreensGeneric(0,y2(j),kappa1(j),k3,omega,U,c,l1,l2,delta,N_bl,HP);
    I0(j) = Q0(j).*Y.p./D_l(j);
end

end
