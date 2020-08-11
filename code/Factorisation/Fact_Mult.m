function [Gp,Gm] = Fact_Mult(G,contour,prec,branch_cut)

% Factorises G (with G -> 1 at infinity) multiplicatively.

% Assumes log(G) is reasonably behaved... (i.e. no looping around the
% origin, and even stronger not crossing the branch cut specified

% See Fact_Additive for contour, prec setup.

F = @(z) LogA(G(z),branch_cut);

[Fp,Fm] = Fact_Additive(F,contour,prec) ;

Gp = @(k) exp(Fp(k));
Gm = @(k) exp(Fm(k));

end