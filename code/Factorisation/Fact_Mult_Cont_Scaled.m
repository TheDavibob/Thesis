function [g] = Fact_Mult_Cont_Scaled(G,Gff,contour,prec,branch_cut)
% Multiplicatively factorised G, with far-field behaviour of G known (but
% not unity)

% Gff a structure with Gff.p and Gff.m  factorisations of Gff.f, known
% analytically.

% Output g contains g.f, g.p and g.m in a nice confined form, analytically
% continued.

H = @(k) G(k)./Gff.f(k); % Should scale to unity at infinity

[Hplus,Hminus] = Fact_Mult_Cont(H,contour,prec,branch_cut);

g.f = @(k) G(k);
g.p = @(k) Hplus(k).*Gff.p(k);
g.m = @(k) Hminus(k).*Gff.m(k);

end

