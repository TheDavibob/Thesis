function [I0] = Distributed_WallPressure(y2,kappa1,shift_y2,omega,k3,U,c,delta,N_bl,l1,l2,HP)
% Computes the wall pressure due to a flow above a wall with linear
% boundary conditions l1 phi' + l2 phi = 0, at convected wavenumber.

% Note, y2 and kappa1 should be linked via kappa1 = omega/U(y2). Since
% kappa1 is generally what we're integrating over, this will need to be
% inverted.

% THIS CODE IS SLOW: we really only want to run it once. Once it's done,
% it's done, however.

if HP > 0
    I = 1;
elseif HP < 0
    I = -1;
else
    I = 0;
end

I0 = zeros(size(y2));

h = waitbar(0,'Setting up distributive integral: wall-pressure');
for j = 1:numel(y2)
    if I == 0
        z_outer.f = @(s) delta*(1-s) + (y2(j)+shift_y2)*s; % so from delta to y
        z_outer.df = @(s) ((y2(j)+shift_y2) - delta)*ones(size(s));
    else
        z_outer.f = @(s) (y2(j)+shift_y2) + ((delta-(y2(j)+shift_y2))/2)*(1+ exp(I*1i*pi*s)) ;
        z_outer.df = @(s) ((delta-(y2(j)+shift_y2))/2)*I*1i*pi*exp(I*1i*pi*s);
    end
    [dec_y2] = Comp_Dec_Aug18A_skip(kappa1(j),k3,omega,U,c,delta,N_bl,z_outer) ;

    if I == 0
            z_0.f = @(s) delta*(1-s); % so from delta to y
            z_0.df = @(s) ( - delta)*ones(size(s));
    else
            z_0.f = @(s) ((delta)/2)*(1+ exp(I*1i*pi*s)) ;
            z_0.df = @(s) ((delta)/2)*I*1i*pi*exp(I*1i*pi*s);
    end
    [dec_0] = Comp_Dec_Aug18A_skip(kappa1(j),k3,omega,U,c,delta,N_bl,z_0) ;
  
    
    C_0 = 1i*(omega - U.f(0)*kappa1(j));
    C_y2 = 1i*(omega - U.f(y2(j)+shift_y2)*kappa1(j));

    P_0 = l1(kappa1(j)); % Just canceling out the C_0^3

%    Rat = C_y2.^3./C_0.^4;
    Rat = 1./(C_0.*c.f(0).^2); % Absorbing y2 into p, rather than phi.
    I0(j) = Rat.*P_0.*dec_y2.p./(l1(kappa1(j)).*dec_0.dphi + l2(kappa1(j)).*dec_0.phi);
    
    waitbar(j./numel(y2),h);
end
close(h);

