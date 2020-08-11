function [I] = PointwiseTrapezium(f,x)
%PointwiseTrapezium  Integrates a pointwise defined function (eg. via an
%integral) along an set of points, given by vector x.

%   Straightforward if inefficient use of the trapezium rule.

F=zeros(size(x));
for j=1:numel(F)
    F(j)=f(x(j));
end

dI=0.5*(F(2:end)+F(1:end-1)).*(x(2:end)-x(1:end-1));
I=sum(dI(:));

end

