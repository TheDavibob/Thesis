function z = Parabola_3Points(zm,z0,zp)
% A (complex parabola), with z( pm 1 ) = zpm, and z(0) = z0, as a function
% of t runnning along the real axis.

A = 1/2*(zm + zp) - z0;
B = 1/2*(zp - zm);
C = z0;

z.f = @(t) A*t.^2 + B*t + C;
z.df = @(t) 2*A*t + B;
z.d2f = @(t) 2*A*ones(size(t));

end