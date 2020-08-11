function [F_bar] = Fbar(k1,kappa1,dkappa1,Int)
%F_bar: A simple vector construction for F_bar, given Int prescribed for
%kappa1
% Most of the work is in computing Int, which should be done separately

% Makes everything the same size
Int = Int(:).';
k1 = k1(:).';
kappa1 = kappa1(:).';
dkappa1 = dkappa1(:).';

% Constructs an array of Cauchy like terms
[K1,Kappa1] = meshgrid(k1,kappa1);
Cauchy = 1./(K1 - Kappa1);
% Rows k1, columns kappa1;

% Integration
F_bar = (dkappa1.*Int)*Cauchy;
end

