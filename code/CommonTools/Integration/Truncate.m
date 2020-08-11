function Lt=Truncate(L,x,N)

% truncates a function L(x) at +/- N, in case it goes NaNny. This means we
% don't have to worry about integration points outside this region too
% much.

Lt=zeros(size(x));
Lt(abs(x)<N)=L(x(abs(x)<N));

end