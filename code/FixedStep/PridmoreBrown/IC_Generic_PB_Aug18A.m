function [Phi0]=IC_Generic_PB_Aug18A(k1,l1,l2)
% Initial conditions for the generic PB problem in the simplest
% formulation

% on x2=0, l1 phi' + l2 phi = 0;
% l1, l2 anonymous functions of k1
L1 = l1(k1);
L2 = l2(k1);

Phi0 = zeros(2,numel(k1));
Phi0(1,:) = L1 ;
Phi0(2,:) = -L2;
% Normalised so that W(0) = D_l(0);


% TEST

end
