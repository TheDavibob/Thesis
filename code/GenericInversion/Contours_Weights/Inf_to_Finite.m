function [C] = Inf_to_Finite(N)
%INF_TO_FINITE Gauss-Legendre weights for an infinite integral, mapped to
%the interval (-1,1) through a transformation controlled by alpha.

[s,w_s] = lgwt(N,-1,1);
s = s(end:-1:1).';
w_s = w_s.';

C.s = s;
C.ws = w_s;
C.k = s./sqrt(1-s.^2);
C.dk = w_s./(sqrt(1-s.^2).^3);

end

