function [k1] = Steepest_Descent_Contour(x1,x2,omega,Uinf,cinf,k3,N)
% Finds the steepest descent contour, as a function of location x1 and x2,
% with N quadrature points and corresponding weights.

% Not yet tailored for negative x2, which flips various things.

% Allows both vector and meshgrid input of coordinates.
if isvector(x1) == 1
    N1 = numel(x1);
    N2 = numel(x2);
    [X1,X2] = meshgrid(x1,x2) ;
else
    N12 = size(x1);
    N1 = N12(2); % Assuming meshgrid so coordinates flipped.
    N2 = N12(1);
    X1 = x1;
    X2 = x2;
end

% Flow derived parameters
M = Uinf./cinf;
k0 = omega./cinf;
beta = sqrt(1-M^2); % assumed positive and real

% The angle corresponding to X1 and X2, with Doppler shift
Theta = atan2(beta*X2,X1); % size N2 x N1
% between -+pi

K0 = k0/(1-M^2) *sqrt(1 - beta^2*k3^2/k0^2); % Modified acoustic wavenumber.
% Watch out: cut-off for large k3.

[t,w] = lgwt(N,-1,1); % Gauss-Legendre weights for quadrature
t = t(end:-1:1);


THETA = repmat(Theta,1,1,N); % size N2 x N1 x N
T = permute(repmat(t,1,N2,N1),[2,3,1]); % size N2 x N1 x N

% Steepest descent contour in transformed space.
K1_SD = K0*(cos(Theta) + T./sqrt(1-T.^2).*sin(THETA) - 1i*(T.^2./sqrt(1-T.^2) .* cos(THETA) - T.*sin(THETA) ));
dK1 = K0*(1./(sqrt(1-T.^2)).^3.*sin(THETA) + 1i*(T.*(T.^2 - 2)./(sqrt(1-T.^2)).^3 .* cos(THETA) + sin(THETA) ));


% Steepest descent contour in "real" k1 space
k1_SD = K1_SD - M*k0/(1-M^2);
dk1 = dK1;

W = permute(repmat(w,1,N2,N1),[2,3,1]); % size N2 x N1 x N, epanded weights
WdT = W.*dk1;

% Contour and space parameters
k1.SD = k1_SD;
k1.t = t;
k1.dt = dk1;
k1.wdt = WdT;
k1.w = w;
k1.X1 = X1;
k1.X2 = X2;
k1.Theta = Theta;

% Physical parameters
k1.k0 = k0;
k1.K0 = K0;
k1.beta = beta;
k1.M = M;

end
