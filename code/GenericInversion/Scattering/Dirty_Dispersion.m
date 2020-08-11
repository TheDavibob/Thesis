function [Dp] = Dirty_Dispersion(k1,k3,omega,U,c,delta,N_bl)
%DIRTY_DISPERSION Summary of this function goes here
%   Detailed explanation goes here

z.f = @(s) delta*(1-s);
z.df = @(s) -delta;
dec = Comp_Dec_Aug18A_skip(k1,k3,omega,U,c,delta,N_bl,z);
Dp = dec.p;

end

