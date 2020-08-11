function [I1,I2]=BranchCutIntegral(f,B,delta)

% Integrates on either side of the branch cut, C1 on the left/upper side
% from bottom left to top right, C2 the other way around. If there is no
% branch cut, ie B describes a simple line of points, we would expect C1+C2
% to vanish.

norm=-1i*(B(2:end)-B(1:end-1));
norm=[norm, norm(end)];
norm=norm./abs(norm);

C2=B+delta*norm;
C1=B-delta*norm;

% Prefers C1 to run from top to bottom and C2 from bottom to top. If not
% possible, C1 runs from left to right.
if imag(C1(1))<imag(C1(end))
    C1=C1(end:-1:1);
elseif imag(C1(1))>imag(C1(end))
    C2=C2(end:-1:1);
elseif real(C1(1))>real(C1(end))
    C1=C1(end:-1:1);
else
    C2=C2(end:-1:1);
end
    
F1=f(C1); F2=f(C2);

i1=0.5*(C1(2:end)-C1(1:end-1)).*(F1(2:end)+F1(1:end-1));
I1=sum(i1);

i2=0.5*(C2(2:end)-C2(1:end-1)).*(F2(2:end)+F2(1:end-1));
I2=sum(i2);

end