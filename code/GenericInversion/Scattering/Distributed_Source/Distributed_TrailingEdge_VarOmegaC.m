%% Supposing only omega_c is changed, this recomputes everything that needs recomputing

% omega_c = ones(size(y2));
% omega_c(40) = 1;
% omega_c(60) = -1;

omega_c = y2.*(deltap - y2).*ones(size(y2));

%% Upstream
[phi_0] = Upstream_Distributed_Integrated(x1,omega_c,G_01,dkappa1);
I_upstream = phi_0;

%% F_bar
Int = -omega*cp.f(y2).^2.*I0.*omega_c./(Km_distributed.*Up.f(y2).^2);
F_bar = @(k1) Fbar(k1,kappa1,dkappa1,Int);

%% Cauchy-like term
Cauchy_Cp1 = F_bar(Cp1.k);
Cauchy_Cp2 = F_bar(Cp2.k);
Cauchy_Cm1 = F_bar(Cm1.k);
Cauchy_Cm2 = F_bar(Cm2.k);
Cauchy_CLp = F_bar(C_CLp.k);
Cauchy_CLm = F_bar(C_CLm.k);
disp('Computed forcing from upstream terms');

%% Leading coefficient
A = 1./(1i); % Much simpler than before

%% k1-dependence

% The "steepest descent" contours
K_Cp1 = A*Cauchy_Cp1./(DZ1p.*Kernel_p_p1);
K_Cp2 = A*Cauchy_Cp2./(DZ2p.*Kernel_p_p2);
K_Cm1 = A*Cauchy_Cm1.*DZ2m.*Kernel_m_m1./(Dw1m);
K_Cm2 = A*Cauchy_Cm2.*DZ1m.*Kernel_m_m2./(Dw2m);

% The "critical-layer" contours
K_CLp1 = A*Cauchy_CLp.*DZ2_CLp.*Kernel_m_CLp./(Dw_CLp);
K_CLm1 = A*Cauchy_CLm.*DZ2_CLm.*Kernel_m_CLm./(Dw_CLm);
K_CLp2 = A*Cauchy_CLp.*DZ1_CLp.*Kernel_m_CLp./(Dw_CLp);
K_CLm2 = A*Cauchy_CLm.*DZ1_CLm.*Kernel_m_CLm./(Dw_CLm);

disp('Constructed integrands');

%% Integration: acoustic

disp('Computing steepest descent integral');
I_generic = struct('phi',[],'dphi',[],'p',[],'v',[]);
physical_variables = fieldnames(I_generic);
physical_sign = [1,-1,1,-1];
% For x2 < 0, some of the signs need to be flipped. This does it.
% If any fields are deleted, delete their sign, too.

I_acoustic_Cp1 = I_generic;
I_acoustic_Cp2 = I_generic;
I_acoustic_Cm1 = I_generic;
I_acoustic_Cm2 = I_generic;
I_acoustic_1 = I_generic;
I_acoustic_2 = I_generic;
I_acoustic = I_generic;
for j = 1:numel(physical_variables)
    s = physical_variables{j};
    I_acoustic_Cp1.(s) = Integral_Nov2018(Cp1.wk,K_Cp1,Xp1,Yp1.(s));
    I_acoustic_Cp2.(s) = Integral_Nov2018(Cp2.wk,K_Cp2,Xp2,Yp2.(s));
    I_acoustic_Cm1.(s) = Integral_Nov2018(Cm1.wk,K_Cm1,Xm1,Ym1.(s));
    I_acoustic_Cm2.(s) = Integral_Nov2018(Cm2.wk,K_Cm2,Xm2,Ym2.(s));
    I_acoustic_1.(s) = [I_acoustic_Cp1.(s);I_acoustic_Cm1.(s)];
    I_acoustic_2.(s) = physical_sign(j)*[I_acoustic_Cp2.(s)(:,end:-1:1);I_acoustic_Cm2.(s)(:,end:-1:1)];
    I_acoustic.(s) = [I_acoustic_2.(s),I_acoustic_1.(s)];
end

%% Integration: critical-layer
disp('Computing critical-layer integrals');

I_CLp1 = I_generic;
I_CLp2 = I_generic;
I_CLm1 = I_generic;
I_CLm2 = I_generic;
I_CL1 = I_generic;
I_CL2 = I_generic;
I_CL = I_generic;

for j = 1:numel(physical_variables)
    s = physical_variables{j};
    I_CLp1.(s) = Integral_Nov2018(C_CLp.wk,K_CLp1,X_CLp,Y_CLp1.(s));
    I_CLp2.(s) = Integral_Nov2018(C_CLp.wk,K_CLp2,X_CLp,Y_CLp2.(s));
    I_CLm1.(s) = Integral_Nov2018(C_CLm.wk,K_CLm1,X_CLm,Y_CLm1.(s));
    I_CLm2.(s) = Integral_Nov2018(C_CLm.wk,K_CLm2,X_CLm,Y_CLm2.(s));
    I_CL1.(s) = I_CLp1.(s)+I_CLm1.(s);
    I_CL2.(s) = physical_sign(j)*I_CLp2.(s)(:,end:-1:1)+physical_sign(j)*I_CLm2.(s)(:,end:-1:1);
    I_CL.(s) = [I_CL2.(s),I_CL1.(s)];
end


%% Integration: modal contributions
disp('Determining modal contributions')


% Downstream wake modes
N_wake_m = numel(zeros_wake_m);

%[Y1_wake_m] = Y_Continuous_Scattering_Nov2018_Mod(x2(x2>=0),zeros_wake_m,k3,omega,Up,cp,deltap,N_bl);
%[Y2_wake_m] = Y_Continuous_Scattering_Nov2018_Mod(sort(-x2(x2<0)),zeros_wake_m,k3,omega,Um,cm,deltam,N_bl);

%X_wake_m = Fourier_Nov2018(x1(x1>=0),zeros_wake_m);

K1_wake_m = -2*pi*1i*A*DZ2(zeros_wake_m,0).*K_factorised.m(zeros_wake_m).*F_bar(zeros_wake_m)./(dDw_wake_m);
K2_wake_m = -2*pi*1i*A*DZ1(zeros_wake_m,0).*K_factorised.m(zeros_wake_m).*F_bar(zeros_wake_m)./(dDw_wake_m);

I_wake_m1 = cell(N_wake_m,1);
I_wake_m2 = cell(N_wake_m,1);
I_wake_m = cell(N_wake_m,1);

I_wake_m_total = I_generic;
for j = 1:numel(physical_variables)
    s = physical_variables{j};
    I_wake_m_total.(s) = zeros(numel(x1(x1>=0)),numel(x2));
end

for k = 1:N_wake_m
    I_wake_m1{k} = I_generic;
    I_wake_m2{k} = I_generic;
    I_wake_m{k} = I_generic;
    for j = 1:numel(physical_variables)
        s = physical_variables{j};
        I_wake_m1{k}.(s) = K1_wake_m(k)*X_wake_m(:,k)*Y1_wake_m.(s)(:,k).';
        I_wake_m2{k}.(s) = physical_sign(j)*K2_wake_m(k)*X_wake_m(:,k)*Y2_wake_m.(s)(:,k).';
        I_wake_m{k}.(s) = [I_wake_m2{k}.(s)(:,end:-1:1),I_wake_m1{k}.(s)];
        I_wake_m_total.(s) = I_wake_m_total.(s) + I_wake_m{k}.(s);
    end
end


%% Combinations

Z11 = zeros(numel(x1(x1<0)),numel(x2(x2>=0)));
Z12 = zeros(numel(x1(x1<0)),numel(x2(x2<0)));
Z21 = zeros(numel(x1(x1>=0)),numel(x2(x2>=0)));
Z22 = zeros(numel(x1(x1>=0)),numel(x2(x2<0)));

% Firstly: extending various things to cover all x1, x2
I_modal = I_generic;
I_crit = I_generic;
I_upstream_extended = I_generic;

for j = 1:numel(physical_variables)
    s = physical_variables{j};
    I_modal.(s) = [ [Z12,Z11];I_wake_m_total.(s)] ;
    I_crit.(s) = [ [Z12,Z11];I_CL.(s)];
    I_upstream_extended.(s) = [ [Z12;Z22],I_upstream.(s)];
end

% Secondly: adding things together in a variety of ways
I_scattered = I_generic;
I_total = I_generic;

for j = 1:numel(physical_variables)
    s = physical_variables{j};
    I_scattered.(s) = I_modal.(s) + I_crit.(s) + I_acoustic.(s) ;
    I_total.(s) = I_scattered.(s) + I_upstream_extended.(s);
end
