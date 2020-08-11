function [K] = K_Hard(y2,k1,k3,omega,U,c,delta,N_bl,HP)
%K_HARD The position invariant contributions for a point source in ctsly
%sheared flow.

% HP determines which half-plane the integration contour passes through.

% Three components:
% a) Cauchy

kappa1 = omega./U.f(y2);
K_Cauchy = 1/(2*pi*1i) .* 1./(k1(:).' - kappa1);

% b) Dispersion relationship
if HP == 0
    z.f = @(s) delta*(1-s);
    z.df = @(s) -delta*ones(size(s));
else
    if HP > 0
        I = 1;
    else
        I = -1;
    end
    z.f = @(s) (delta/2)*(1+ exp(I*1i*pi*s)) ;
    z.df = @(s) (delta/2)*I*1i*pi*exp(I*1i*pi*s);
end
[dec] = Comp_Dec_Aug18A_skip(k1,k3,omega,U,c,delta,N_bl,z) ;

D_h = dec.v(:).';

% c) Convected thingies: ignoring what's happing with the speed of sound
% for now.
C = @(x2) 1i*(omega - U.f(x2)*k1(:).'); % for scalar x2
Wscale_y2 = c.f(y2).^2.*C(y2).^4;
Wscale_0 = c.f(0).^2.*C(0).^4;
Phi_Scale = -C(0).^5; % oopsie: this is due to how phih was defined.
K_Wronskian = Wscale_y2.*Phi_Scale./(Wscale_0.*D_h);

K = K_Wronskian.*K_Cauchy;

end

