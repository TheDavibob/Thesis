function [p]=Poles(f,p0)

% Finds numel(p0) poles of f, using a NewtonRaphson method, with p0 being
% the original guess, probably preferably real. The guess needs to be good
% or else will not work.


F=@(z) 1./f(z);
dF=@(z) PointwiseDerivative(F,z);

p=zeros(size(p0));
p(1) = NewtonRaphson(F,dF,p0(1));

for j=2:numel(p0)
    F=@(z) F(z)./(z-p(j-1)); % removing the previous pole
    dF=@(z) PointwiseDerivative(F,z);
    p(j)  = NewtonRaphson(F,dF,p0(j));
end

end