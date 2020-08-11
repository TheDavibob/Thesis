function [phi]=Fixed_Step_IVP_GridX_skip(f,phi0,x_max,x_N,k)
% Integrates dphi/dx = f(x,phi;k1) using a RK method with fixed step size

% Instead of integrating to a single fixed x_max, integrates to a maximum
% which is the same size as k. The SKIP version doesn't save intermediate
% values, and thus reduces things being carried around unnecessarily.

% The equation is setup as d/dx phi = f(x,phi;k) where phi and f are
% phi_N-vectors and k is a parameter of arbitrary size k_N.

% Now need input X to be an k_N x x_N sized matrix, which shouldn't be too
% challenging, with output the same size.

% In practice, size(phi0) = phi_N x k_N ;
% size(phi) = x_N x phi_N x k_N ;
% size(f) = phi_N x k_N.

% Integrates from 0, with phi(0)=phi0(k), a function of k, to x_max, with
% x_N steps.

% Uses a fourth-order RK routine.

K = k(:);
X_max = x_max(:);
k_N = numel(K);

if numel(X_max) ~= numel(K)
    error('Discrepency between number of parameter values and end point values')
end

X_max_stretch = repmat(X_max,1,x_N);
K_stretch = repmat(K,1,x_N);
Scale_stretch = repmat(linspace(0,1,x_N),k_N,1);
X = X_max_stretch.*Scale_stretch;

%x = linspace(0,x_max,x_N) ; % x points evenly distributed
% h = x(2)-x(1) ; % step size

H = X(:,2) - X(:,1) ; % step size

%k_N=numel(k);

IC=phi0(K); % Needs N rows and K columns, where N is the order of the ode
phi_N = size(IC,1) ;

phi=zeros(phi_N,k_N); % sets up an array for phi
phi(:,:) = IC;

% Going to need a little loop
for j=2:x_N
    c1 = H.'.*f(X(:,j-1) , phi , K ) ;
    c2 = H.'.*f(X(:,j-1) + H/2, phi+c1/2 , K ) ;
    c3 = H.'.*f(X(:,j-1) + H/2, phi+c2/2 , K ) ;
    c4 = H.'.*f(X(:,j-1) + H , phi + c3 , K) ;
    % Note, these constants are phi_N x k_N arrays
    phi = phi + (c1+2*c2+2*c3+c4)/6 ;
end

end