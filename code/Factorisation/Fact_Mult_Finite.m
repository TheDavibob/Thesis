function [Gp,Gm] = Fact_Mult_Finite(G,contour,prec,branch_cut)

% Factorises G (with G -> 1 at infinity) multiplicatively.

% Assumes log(G) is reasonably behaved... (i.e. no looping around the
% origin, and even stronger not crossing the branch cut specified

% See Fact_Additive for contour, prec setup.

if isa(branch_cut,'function_handle') % This functionality should be added to the others. Added Nov2018
    F = @(z) Log_CurvedBranch(G(z),branch_cut);
else
    F = @(z) LogA(G(z),branch_cut);
end

[Fp,Fm] = Fact_Additive_Finite(F,contour,prec) ;

Gp = @(k) exp(Fp(k));
Gm = @(k) exp(Fm(k));

end