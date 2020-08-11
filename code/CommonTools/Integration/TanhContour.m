function [C,dC] = TanhContour(a,b,tm,tp)
%TanhContour - A tanh shaped (in the complex plane) contour

% The contour centres between tp and tm, tending towards a pm b at pm
% infinity ( so ideally want b>a to avoid things)

C=@(t) t + 1i*( a + b*tanh( 2/(tp-tm) * (t - (tp+tm)/2) ) );
dC=@(t) 1+1i*b*2/(tp-tm)*sech( 2/(tp-tm) * (t - (tp+tm)/2) ).^2;

end

