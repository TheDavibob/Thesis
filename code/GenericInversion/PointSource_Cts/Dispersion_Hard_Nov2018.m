function Dh = Dispersion_Hard_Nov2018(k1,k3,omega,U,c,delta,N_bl,HP)

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

Dh = dec.v(:).';
Dh = reshape(Dh, size(k1));

end

