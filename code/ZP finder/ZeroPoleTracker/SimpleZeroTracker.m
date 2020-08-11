function [z] = SimpleZeroTracker(f,p,z0)
%SimpleZeroTracker Tracks a zero of f(z,p) as p is varied

% Finds a zero from guess z0, for first element of p. For successive
% elements, finds zero using previous position as guess.

z=zeros(size(p));
F=@(z) f(z,p(1));
dF=@(z) PointwiseDerivative(F,z);

h=waitbar(0,'Hold on');
z(1)=NewtonRaphson(F,dF,z0);


for j=2:length(p)
    F=@(z) f(z,p(j));
    dF=@(z) PointwiseDerivative(F,z);
    z(j)=NewtonRaphson(F,dF,z(j-1));
    if isnan(z(j))==1
        z(j+1:end)=NaN;
        break
    end
    plot(real(z(1:j)),imag(z(1:j)));
    drawnow limitrate
    waitbar(j/length(p));
end

plot(real(z),imag(z));
close(h) 
end

