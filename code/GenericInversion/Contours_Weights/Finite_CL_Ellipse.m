function [C_CLp,C_CLm] = Finite_CL_Ellipse(N,a,b,shift1,shift2)
%FINITE_CL_ELLIPSE Finds integration contour around a branch cut

% Branch cut from a to b, shifted at either end and upwwards/downwards

A = 1/2*(a + b);
B = 1/2*(b-a) + shift1;
C = 1i*shift2;

[s,ws] = lgwt(N,-1,1);
s = s(end:-1:1).';
ws = ws.';

C_CLp.k = A + B*sin(pi*s/2) + C*cos(pi*s/2);
C_CLm.k = A - B*sin(pi*s/2) - C*cos(pi*s/2);

C_CLp.wk = pi/2*ws.*(B*cos(pi*s/2) - C*sin(pi*s/2));
C_CLm.wk = pi/2*ws.*(- B*cos(pi*s/2) + C*sin(pi*s/2));

end

