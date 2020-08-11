function [I] = PointwiseTrapeziumBox(f,x,y)
%PointwiseTrapeziumBox Integrates around the box z=x+iy

% x,y increasing vectors, need not have constant spacing.

% Direct extension of PointwiseTrapezium. Integrates anticlockwise.

I1=PointwiseTrapezium(f,x+1i*y(1));
I2=PointwiseTrapezium(f,x(end)+1i*y);
I3=-PointwiseTrapezium(f,x+1i*y(end)); % note direction of integration
I4=-PointwiseTrapezium(f,x(1)+1i*y);
I=I1+I2+I3+I4;

end

