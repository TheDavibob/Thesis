function [I] = Integral_Nov2018(wk,K,X,Y)
%Computes an integral of the form int_C (K(k) X(x,k) Y(y,k) dk

% wk1 contains the relevant weights from the contour of integration: note
% we don't actually need the points, since they've already been used.

% K is of size of k, X of size x by k, and so on. All need to be
% precomputed, which should be fine.

% Integrates as a loop over k, which isn't too challenging.

N_k = numel(K);
N_x = floor(numel(X)./N_k);
N_y = floor(numel(Y)./N_k);

if numel(wk) ~= N_k
    error('Number of weights differs from number of points integrated over');
end

% Computes the integral as a loop over k. For this reason, N_k should be
% reasonably small (maybe 100), with points clustered near interesting
% locations.

I = zeros(N_x,N_y);
h = waitbar(0,'Performing integration loop');
for j = 1:N_k
    I = I + wk(j)*K(j)*X(:,j)*Y(:,j).';
    waitbar(j./N_k,h,'Performing integration loop');
end
close(h);

