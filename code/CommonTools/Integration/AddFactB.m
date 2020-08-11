function [Fp] = AddFactB(F,C,dC,N1,N2,N3,t1,t2)
%AddFact - Finds the part of (suitably decaying) F that is analytic in either
% half-plane, half-planes defined by the contour C

% Slightly modified original, this created 2017-07-27
% Major sign correction 2017-09-07

G=@(t) F(C(t)).*dC(t); % Note "arbitrary" choice of where dC lies
f=@(t,x) CauchyKernel(C(t),x); % this has the wrong sign

Fp2=FullIntervalIntA(f,G,N1,N2,N3,t1,t2);
Fp=@(k) -Fp2(k);
% THIS HAS BEEN UPDATED TO ACCOUNT FOR THE ERROR IN THE CAUCHY KERNEL
end
