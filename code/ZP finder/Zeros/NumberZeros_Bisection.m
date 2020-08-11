function [N,Z,sz,X,Y] = NumberZeros_Bisection(fun,x,y,Nh,Nv,recurse,precision)
% Computes the rough location of the zeros in a box, via a bisection method

% Finds the zeros count in a large box, splits it into quarters and
% recurseively finds the zeros in each following stage.

% sz is the sum of zeros in each of the boxes -- if N is one for that box,
% sz is the location of the single zero

% recurse an integer, indicating how many times the box is divided into
% two.

% fun a structure-type function, as per NumberZeros_Box.

fine = 2^recurse;

if numel(x) ~= 2
    error('Box must be given by two x and y coordinates')
elseif numel(y) ~= 2
    error('Box must be given by two x and y coordinates')
end

x = sort(x);
y = sort(y);

dx = (x(2)-x(1))/fine;
dy = (y(2)-y(1))/fine;

X = x(1):dx:x(2);
Y = y(1):dy:y(2);
[Xg,Yg] = ndgrid(X,Y);
Z2 = Xg +1i*Yg;
Z = 0.5 * (Z2(1:end-1,1:end-1) + Z2(2:end,2:end));

N2 = NumberZeros_Box(fun,x,y,Nh,Nv,precision) ; % Total number of zeros
N = N2*ones(fine,fine); % A big full grid of N. Going to recursively reduce.

for j = 1:recurse
    index = 1:fine/(2^j):fine+1;
    for k = 1:numel(index)-1
        for l = 1:numel(index)-1
            if N(index(k),index(l)) ~= 0
                N2 = NumberZeros_Box(fun,[X(index(k)),X(index(k+1))],[Y(index(l)),Y(index(l+1))],Nh,Nv,precision) ;
                N(index(k):1:index(k+1)-1,index(l):1:index(l+1)-1) = N2*ones(size(fine/(2^j),fine/(2^j)));
            end
        end
    end
    disp(['Completed ',num2str(j),'th bisection'])
end

if max(N(:)) > 1
    disp('Some boxes have multiple zeros')
    sz = Z;
else
    
    % Z(N==j) gives the centre of each box with j zeros
    
    fun.zdff = @(z) z.*fun.df(z)./(2*pi*1i*fun.f(z));
    
    sz = zeros(size(N));
    extra_prec = 1; % Extra precision in the integral to get them bang on. Currently removed as it's stupid.
    
    for k=1:fine
        for l=1:fine
            if N(k,l) ~= 0
                sz(k,l) = ComplexIntegral_Box(fun.zdff, [X(k),X(k+1)], [Y(l),Y(l+1)], extra_prec*Nh, extra_prec*Nv);
            end
        end
    end

end
% If N is only ever 1, sz(N == 1) is a more precise guess of the zeros.
% There might be problems with integration at this scale?

end