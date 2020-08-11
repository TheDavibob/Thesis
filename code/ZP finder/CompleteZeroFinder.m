function [z]=CompleteZeroFinder(fun,dfun,x,y,Plot)

% finds all zeros of fun within the box of interest. Only works if
% there are no poles (though, feasibly, could extend to include pole
% finding too.) Uses a bisectiony sort of method.

% Need x=[x1 x2] etc

Int=BoxIntegral(fun,x,y);
if abs(Int) > 0.1
    error('the function contains singularities within the box')
end

    function zb=ZeroFinder2(xb,yb) %Finds up to two zeros in a given box
        zb=[0 0];
        d1=SumZerosj(fun,dfun,xb,yb,1);
        d2=SumZerosj(fun,dfun,xb,yb,2);
        zb(1)=(d1+sqrt(2*d2-d1^2))/2;
        zb(2)=(d1-sqrt(2*d2-d1^2))/2;
    end


x(1,:)=x; y(1,:)=y;
Z(1)=ZerosPoles(fun,dfun,x(1,:),y(1,:));
X=['Number of zeros in region ',num2str(Z(1))];
disp(X);
NoZeros=round(Z(1));

while max(Z)>2.001
    [~,I]=max(Z);
    x(I+4:end+3,:)=x(I+1:end,:); y(I+4:end+3,:)=y(I+1:end,:);
    Z(I+4:end+3)=Z(I+1:end);
    x1=x(I,1); y1=y(I,1); x2=x(I,2); y2=y(I,2);
    xm=(x1+x2)/2; ym=(y1+y2)/2;
    x(I,:)=[x1,xm]; y(I,:)=[y1,ym];
    x(I+1,:)=x(I,:); y(I+1,:)=[ym,y2];
    x(I+2,:)=[xm,x2]; y(I+2,:)=y(I,:);
    x(I+3,:)=x(I+2,:); y(I+3,:)=y(I+1,:);
    for j=0:3;
        Z(I+j)=ZerosPoles(fun,dfun,x(I+j,:),y(I+j,:));
        if min(abs(x(I+j,2)-x(I+j,1)),abs(y(I+j,2)-y(I+j,1)))<0.00001
            Z(I+j)=-1;
        end
    end
end
% The last bit is a break clause in case of large multiplicity
% Run into problems if integrating across the zero

z=zeros(1,NoZeros);
j=1;
for k=1:length(Z)
    if Z(k)==-1
        Mult=ZerosPoles(fun,dfun,x(k,:),y(k,:));
        Mult=round(Mult);
        z(j:j+Mult-1)=(x(k,1)+1i*y(k,1))*ones(1,Mult);
        j=j+Mult;
    elseif Z(k)<0.001
    elseif Z(k)<1.001
        z(j)=SumZerosj(fun,dfun,x(k,:),y(k,:),1);
        j=j+1;
    else
        zb=ZeroFinder2(x(k,:),y(k,:));
        z(j)=zb(1); z(j+1)=zb(2);
        j=j+2;
    end
end

if Plot==1
    scatter(real(z),imag(z),'o');
end

end