function C = Finite_Parabola_3Points(N,zm,z0,zp)
% Finds a parabola passing through three points, maps it to a finite
% interval and finds weights, etc.

z = Parabola_3Points(zm,z0,zp);
C_int = Inf_to_Finite(N);

C.k = z.f(C_int.k); % Points in the complex plane of the parabola.
C.wk = z.df(C_int.k).*C_int.dk; % The weights.

end