function N = NumberZeros_Box(fun,x,y,Nh,Nv,precision)

% Computes the number of zeros of a (assumed complex analytic) function
% within a complex box

% Corners given by x,y as per ComplexIntegral_Box. "precision" gives an
% expected error on the integrals

% fun is now a structure array with f and df components. df may well have
% previously been computed numerically:
% fun.f = @(z) F(z);
% fun.df = @(z) dF/dz(z);

% Section one: Brutish analyticity check, using Cauchy's theorem
I_analytic = ComplexIntegral_Box(fun.f,x,y,Nh,Nv);
if abs(I_analytic) > precision
    error('function is not analytic in the box')
end

% Zero counter: Uses the argument principle
fun.dff = @(z) fun.df(z)./(2*pi*1i*fun.f(z)); % f'/f, scaled out the 2pi i we expect.

N2 = ComplexIntegral_Box(fun.dff,x,y,Nh,Nv);
if abs(N2 - round(N2,0)) > precision
    error('function seems to have a fractional number of zeros: either accuracy is insufficient or branch cuts have been crossed')
end
N = round(N2,0);

end