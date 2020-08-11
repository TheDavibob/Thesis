function [z]=SimpleZeroFinder(fun,dfun,x,y)

% finds up to 2 zeros of fun within the box of interest. Only works if
% there are no poles (though, feasibly, could extend to include pole
% finding too.)

Int=BoxIntegral(fun,x,y);
if abs(Int) > 0.001
    error('the function contains singularities within the box')
end

Z=ZerosPoles(fun,dfun,x,y);

if Z>2.001
    error('more than two zeros found within the box - try moving or shrinking the box')
end

if Z>1.001
    z=[0 0];
    d1=SumZerosj(fun,dfun,x,y,1);
    d2=SumZerosj(fun,dfun,x,y,2);
    z(1)=(d1+sqrt(2*d2-d1^2))/2;
    z(2)=(d1-sqrt(2*d2-d1^2))/2;
else
    d1=SumZerosj(fun,dfun,x,y,1);
    z=d1;
end

end