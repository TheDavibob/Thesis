function I = ComplexIntegral_Box(fun,x,y,Nh,Nv)
% ComplexIntegral_Box - Integrates around a box, orientated with the real
% and imaginary axes.

% Uses standard Gaussian quadrature methods to compute int fun dz around
% the box with opposite corners [x1,y1], [y2,x2]; Uses four line integrals,
% the horizontal with accuracy Nh, the vertical with accuracy Nv.

if numel(x) ~= 2
    error('Box must be given by two x and y coordinates')
elseif numel(y) ~= 2
    error('Box must be given by two x and y coordinates')
end

X = sort(x);
Y = sort(y); 

[X2,Y2] = ndgrid(X(:),Y(:));
Z = X2 + 1i*Y2; 
% So Z(i,j) = x(i) + i*y(j);

% Integrate each side of the box in turn, progressing anticlockwise
I1 = ComplexIntegral_Line(fun,Z(1,1),Z(2,1),Nh) ;
I2 = ComplexIntegral_Line(fun,Z(2,1),Z(2,2),Nv) ;
I3 = ComplexIntegral_Line(fun,Z(2,2),Z(1,2),Nh) ;
I4 = ComplexIntegral_Line(fun,Z(1,2),Z(1,1),Nv) ;

I = I1+I2+I3+I4;
% A little clunky, but this isn't going to need generalisation

end