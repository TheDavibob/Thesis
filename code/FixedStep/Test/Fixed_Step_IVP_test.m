function [phi,x] = Fixed_Step_IVP_test(x_max,x_N,k)

% Uses a damped oscillator to test Fixed_Step_IVP

% Equation: y'' + 2k y' + y = 0 ; y(0)= 1, y'(0)=0 ;

% With phi(1) = phi and phi(2) = y', can construct f.

    function f = OdeFun(x,phi_jm1,k)
        f=zeros(2,numel(k)); % with phi having two components
        f(1,:) = phi_jm1(2,:) ;
        f(2,:) = -phi_jm1(1,:) - 2*k.*phi_jm1(2,:) ;
    end
% no clever tricks required

% Initial conditions - should be straightforward to make functions of k
    function phi0 = IC(k)
        phi0 = zeros(2,numel(k));
        phi0(1,:) = ones(1,numel(k)) ;
        phi0(2,:) = zeros(1,numel(k)) ;
%        phi0(1,:) = zeros(1,numel(k)) ;
%        phi0(2,:) = ones(1,numel(k)) ;
    end

[phi,x]=Fixed_Step_IVP(@(x,phi,k) OdeFun(x,phi,k) ,@(k) IC(k),x_max,x_N,k);

end