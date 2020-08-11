function D_l = Dispersion_Generic_Nov2018(k1,k3,omega,U,c,l1,l2,delta,N_bl,HP)

% This is a throwaway, for zero finding

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

D_l = l1(k1).*dec.dphi(:).' + l2(k1).*dec.phi(:).';

end

