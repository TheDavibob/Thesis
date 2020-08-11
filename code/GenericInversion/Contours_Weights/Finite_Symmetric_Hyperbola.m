function C = Finite_Symmetric_Hyperbola(N,z0r,apex,curvature,angle)
% Finds a hyperbola (see Symmetric_Hyperbola for deets), maps it to a finite
% interval and finds weights, etc.

z = Symmetric_Hyperbola(z0r,apex,curvature,angle);
C_int = Inf_to_Finite(N);

C.k = z.f(C_int.k); % Points in the complex plane of the hyperbola.
C.wk = z.df(C_int.k).*C_int.dk; % The weights.

end