function [U] = Struct_U_Parabolic_B(u,delta)
%STRUCT_U_PARABOLIC_B A parabola given by three points, with U(deltaj) = uj

C = ( (u(3) - u(1))/(delta(3) - delta(1)) - (u(2) - u(1))/(delta(2) - delta(1)))./(delta(3) - delta(2));
B = (u(2) - u(1))/(delta(2) - delta(1)) - C*(delta(2) - delta(1)) ;


U.f = @(x) u(1) + B*(x - delta(1)) + C*(x - delta(1)).^2;
U.df = @(x) B + 2*C*(x-delta(1));
U.d2f = @(x) 2*C*ones(size(x));

if sign(B) ~= sign(B + 2*C*(delta(3) - delta(1)))
    disp('Parabola is not montonic across the boundary-layer');
end

end

