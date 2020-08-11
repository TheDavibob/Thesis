function fun = Numerical_Derivatives(f,acc)
% Creates a structure of numerical derivatives of f

% fun.f is the original function, with fun.df the first derivative and
% fun.d2f the second. Could feasibly be made more generic.

% acc is the accuracy of the derivative

fun.f = @(z) f(z);
fun.df = @(z) (f(z+acc) - f(z-acc))./(2*acc) ;
fun.d2f = @(z) ( f(z+acc) - 2*f(z) + f(z-acc) )./(acc.^2) ;

end