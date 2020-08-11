function [Z]=NewtonRaphson(f,df,z0)

% Implements a Newton-Raphson scheme to find zeros.

% Accuracy currently fixed. f,df pointwise anon functions.

epsilon=1e-10;
Z=z0;

while abs(f(Z))>epsilon
    if abs(f(Z))>abs(f(z0))
        Z=NaN; % BREAK CLAUSE NEEDS CHECKING
        break
    end
    Z=Z-f(Z)./df(Z);
end

end