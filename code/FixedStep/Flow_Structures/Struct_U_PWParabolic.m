function [U] = Struct_U_PWParabolic(U_points,delta_points)
% Constructs a profile consisting of two parabolas, matched at some
% intermediate point, with no jumps in shear

U0 = U_points(1);
U1 = U_points(2);
Uinf = U_points(3);

delta0 = delta_points(1);
delta1 = delta_points(2);
deltainf = delta_points(3);

A = (U1 - Uinf)/((delta1 - deltainf).^2);
sigma = -2*A*(deltainf - delta1);
B = ( (U0 - U1) - sigma*(delta0 - delta1) ) ./((delta0 - delta1).^2);

Up.f = @(x2) U1 + sigma*(x2- delta1) + A*(x2-delta1).^2;
Up.df = @(x2) sigma + 2*A*(x2-delta1);
Up.d2f = @(x2) 2*A*ones(size(x2));

Um.f = @(x2) U1 + sigma*(x2- delta1) + B*(x2-delta1).^2;
Um.df = @(x2) sigma + 2*B*(x2-delta1);
Um.d2f = @(x2) 2*B*ones(size(x2));

x2p = @(x2) (x2 > delta1);
x2m = @(x2) (x2 <= delta1);

U.f = @(x2) Up.f(x2).*x2p(x2) + Um.f(x2).*x2m(x2);
U.df = @(x2) Up.df(x2).*x2p(x2) + Um.df(x2).*x2m(x2);
U.d2f = @(x2) Up.d2f(x2).*x2p(x2) + Um.d2f(x2).*x2m(x2);

