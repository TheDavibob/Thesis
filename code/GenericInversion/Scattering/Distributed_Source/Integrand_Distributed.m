function Int = Integrand_Distributed(omega,U,c,I0,omega_distributed,Km_distributed)
% Note implicit assumption of constant speed of sound

Int = -omega.*c.f(0)^2.*I0.*omega_distributed./(Km_distributed.*U.f(y2).^2);

end