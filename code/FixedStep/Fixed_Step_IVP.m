function [phi,x]=Fixed_Step_IVP(f,phi0,x_max,x_N,k)
% Integrates dphi/dx = f(x,phi;k1) using a RK method with fixed step size

% The equation is setup as d/dx phi = f(x,phi;k) where phi and f are
% phi_N-vectors and k is a parameter of arbitrary size k_N.

% In practice, size(phi0) = phi_N x k_N ;
% size(phi) = x_N x phi_N x k_N ;
% size(f) = phi_N x k_N.

% Integrates from 0, with phi(0)=phi0(k), a function of k, to x_max, with
% x_N steps.

% Uses a fourth-order RK routine.

x = linspace(0,x_max,x_N) ; % x points evenly distributed
h = x(2)-x(1) ; % step size

k_N=numel(k);

IC=phi0(k); % Needs N rows and K columns
phi_N = size(IC,1) ;

phi=zeros(x_N,phi_N,k_N); % sets up an array for phi

phi(1,:,:) = IC;

% Going to need a little loop
for j=2:x_N
    P = permute(phi(j-1,:,:),[2,3,1]) ; % Stupid Matlab is stupid. Not using squeeze in case either phi_N or k_N are one.
    c1 = h*f(x(j-1) , P , k ) ;
    c2 = h*f(x(j-1) + h/2, P+c1/2 , k ) ;
    c3 = h*f(x(j-1) + h/2, P+c2/2 , k ) ;
    c4 = h*f(x(j-1) + h , P + c3 , k) ;
    % Note, these constants are phi_N x k_N arrays
    phi(j,:,:) = P + (c1+2*c2+2*c3+c4)/6 ;
end


end