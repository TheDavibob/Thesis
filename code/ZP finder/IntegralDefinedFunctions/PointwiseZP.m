function [ZP] = PointwiseZP(f,x,y)
% PointwiseZP - uses PWTrapBox to compute zeros-poles within a box (bounded
% by x,y as PWTrapBox

df=@(z) PointwiseDerivative(f,z);
fun=@(z) df(z)/(2*pi*1i*f(z));

ZP=PointwiseTrapeziumBox(fun,x,y);

end