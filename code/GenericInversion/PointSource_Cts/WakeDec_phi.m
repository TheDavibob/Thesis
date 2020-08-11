function dec_below = WakeDec_phi(k1,k3,omega,Um,cm,deltam,N_bl,z_below)

dec = Comp_Dec_Aug18A_skip(k1,k3,omega,Um,cm,deltam,N_bl,z_below) ;
dec_below = dec.phi;

end