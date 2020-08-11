function K = Wake_Kernel(k1,k3,omega,U1,c1,delta1,U2,c2,delta2,Z,N,z1,z2)

% Upstream we have a lined plate with impedance Z, downstream have a wake
% with continuity of pressure and of velocity. Z can be infinite, or 0
% (though the latter is a clearly distinct case, so be careful).

% Flow above the plate: U1, c1, delta1 and below U2, c2, delta2, though U =
% U2(-x2) etc required.

% N number of boundary-layer integration steps, z the integration contour
% (which might be complex)

[dec1] = Comp_Dec_Aug18A_skip(k1,k3,omega,U1,c1,delta1,N,z1);
[dec2] = Comp_Dec_Aug18A_skip(k1,k3,omega,U2,c2,delta2,N,z2);

% Dispersion relations
C0 = 1i*(omega - U1.f(0)*k1);

Dw = - dec1.p.*dec2.v - dec2.p.*dec1.v;
if Z == inf
    DZ1 = dec1.v;
    DZ2 = - dec2.v;
    K = Dw./(DZ1.*DZ2);
elseif Z == 0 % a fundamentally different case, so be careful
    DZ1 = dec1.p;
    DZ2 = dec2.p;
    K = -Dw./(DZ1.*DZ2);
else
    DZ1 = 1i*omega*Z*dec1.v - C0.*dec1.p;
    DZ2 = - 1i*omega*Z*dec2.v - C0.*dec2.p;
    K = 1i*omega*Z*Dw./(DZ1.*DZ2);
end

end

