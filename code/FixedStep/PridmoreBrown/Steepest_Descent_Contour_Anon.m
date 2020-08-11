function [f,F,K1_SD,k1_SD] = Steepest_Descent_Contour_Anon(omega,Uinf,cinf,k3)

% Determines the steepest descent contour for the free-stream, in terms of
% physical location and free stream variables, parameterised as per the
% thesis.

% Not finished, essentially a check the mathematics is correct. Don't use
% for complicated things.

% x1 and x2 are a meshgrid input

M = Uinf./cinf;
k0 = omega./cinf;

beta = sqrt(1-M^2); % assumed positive and real

% Theta = atan2(beta*x2,x1);
% K1 = @(k1) k1 + M*k0/(1-M^2);
K0 = k0/(1-M^2) *sqrt(1 - beta^2*k3^2/k0^2);

if imag(K0) ~= 0
    disp('Warning: Acoustic wavenumber is complex')
end

[gamma,~,~]=Gamma_FF(k3,omega,cinf,Uinf);

f = @(k1,Theta) 1i*k1*cos(Theta) + gamma(k1)*sin(Theta)./beta;

F = @(K1,Theta) 1i*K1*cos(Theta) + SqrtA(K1 - K0,3*pi/2).*SqrtA(K1 + K0,pi/2).*sin(Theta);
K1_SD = @(t,Theta) K0*(cos(Theta) + t./sqrt(1-t.^2).*sin(Theta) - 1i*(t.^2./sqrt(1-t.^2) * cos(Theta) - t.*sin(Theta) ));

k1_SD = @(t,Theta) K1_SD(t,Theta) - M*k0/(1-M^2);

end