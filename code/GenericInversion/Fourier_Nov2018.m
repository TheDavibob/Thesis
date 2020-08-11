function [X] = Fourier_Nov2018(x,k)
% Fourier inversion factor, as an x by k array
% Convention, waves like -1i*k1.

% Conversion to row vectors
x = x(:).';
k = k(:).';

X = exp(-1i*x.'*k)/(2*pi);

end

