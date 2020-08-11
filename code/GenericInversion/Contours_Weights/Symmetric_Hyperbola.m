function z = Symmetric_Hyperbola(z0r,apex,curvature,angle)
% A complex hyperbola, written as a function of real coordinate

% z0r is the centre on the real axis.
% Apex is the minimal or maximal value of z0i. For this lazy code, it is
% assumed positive apex leads to a curve tending to the lower half plane,
% and vice versa
% angle the argument of z for large argument
% C is some measure of the curvature of the curve. Sharp curve: C small,
% smooth curve, C large

C2 = curvature.^2;
A2 = (tan(angle)).^2;

if apex >= 0
    zi = @(t) (apex + sqrt(C2)) - sqrt(C2 + A2*(t - z0r).^2);
    dzi = @(t) - A2*(t - z0r)./sqrt(C2 + A2*(t - z0r).^2);
else
    zi = @(t) (apex - sqrt(C2)) + sqrt(C2 + A2*(t - z0r).^2);
    dzi = @(t) A2*(t - z0r)./sqrt(C2 + A2*(t - z0r).^2);
end

z.f = @(t) (t+z0r) + 1i*zi(t+z0r);
z.df = @(t) 1 + 1i*dzi(t + z0r);

end