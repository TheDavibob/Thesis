function [C,dC] = TanhContourB(a,bm,bp,tm,tp)
%TanhContourB - A tanh shaped (in the complex plane) contour, split at 0.

% This allows for different behaviour in the right and left half planes, at
% the loss of smoothness at 0.

% both tp and tm should be positive, they are just scaling factors.

fp=@(t) a+(bp-a)*tanh(t/tp); % -> bp as t -> infty
fm=@(t) a-(bm-a)*tanh(t/tm); % -> bm as t -> -infty

dfp=@(t) (bp-a)*sech(t/tp).^2/tp;
dfm=@(t) (bm+a)*sech(t/tm).^2/tm;

C=@(t) t+1i*( heaviside(t).*fp(t)+heaviside(-t).*fm(t) );
dC=@(t) 1+1i*( heaviside(t).*dfp(t)+heaviside(-t).*dfm(t) );
end

