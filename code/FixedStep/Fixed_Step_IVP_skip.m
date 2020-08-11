function [phi]=Fixed_Step_IVP_skip(f,phi0,x_max,x_N,k)
% Integrates dphi/dx = f(x,phi;k1) using a RK method with fixed step size

% Skip: Only retains the last point of integration, i.e. the value at x_max

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

%k_N=numel(k);

IC=phi0(k); % Needs N rows and K columns
%phi_N = size(IC,1) ;

phi = IC;

% Going to need a little loop
for j=2:x_N
    c1 = h*f(x(j-1) , phi , k ) ;
    c2 = h*f(x(j-1) + h/2, phi+c1/2 , k ) ;
    c3 = h*f(x(j-1) + h/2, phi+c2/2 , k ) ;
    c4 = h*f(x(j-1) + h , phi + c3 , k) ;
    % Note, these constants are phi_N x k_N arrays
    phi = phi + (c1+2*c2+2*c3+c4)/6 ;
end


end