%% Computes the wall pressure due to a point vortex at some location y2 above a hard wall

% Base case
U_points_base = [0.6,0.7,1];
delta_points_base = [0,0.1,1];
U_base = Struct_U_PWParabolic(U_points_base,delta_points_base);

% Sheared case
U_points_sheared = [0.2,0.75,1];
delta_points_sheared = [0,0.3,1];
U_sheared = Struct_U_PWParabolic(U_points_sheared,delta_points_sheared);


y2 = 0.1:0.1:0.5;
shift_y2 = 0.005;

delta = delta_points_base(end); % should hopefully be same both base and sheared
Uinf = U_base.f(delta);
cinf = 35;

omega = logspace(0,3,30);
f = 10*omega./(2*pi*0.02);

% Various computational parameters
N_bl = 3000;
k3 = 0;
c0 = Struct_c_Const(cinf);
q0 = 1;

Full_Pressure_base = zeros(numel(omega),numel(y2));
Full_Pressure_sheared = zeros(numel(omega),numel(y2));

for j = 1:numel(omega)
    Omega = omega(j);
    C0_base = @(k1) 1i*(Omega - U_base.f(0)*k1);
    dC1_base = @(k1) -1i*U_base.df(0)*k1;
    V1_base.l1 = @(k1) -C0_base(k1).^2;
    V1_base.l2 = @(k1) -3*C0_base(k1).*dC1_base(k1);
    L1_base = V1_base;    
    
    C0_sheared = @(k1) 1i*(Omega - U_sheared.f(0)*k1);
    dC1_sheared = @(k1) -1i*U_sheared.df(0)*k1;
    V1_sheared.l1 = @(k1) -C0_sheared(k1).^2;
    V1_sheared.l2 = @(k1) -3*C0_sheared(k1).*dC1_sheared(k1);
    L1_sheared = V1_sheared;
    
    [WallPressure_base] = Distributed_WallPressure(y2,Omega./U_base.f(y2),shift_y2,Omega,k3,U_base,c0,delta,N_bl,L1_base.l1,L1_base.l2,-1);
    Pressure_Upstream_base = -q0*WallPressure_base.*c0.f(y2).^2./U_base.df(y2);
    Full_Pressure_base(j,:) = Pressure_Upstream_base;
    
    [WallPressure_sheared] = Distributed_WallPressure(y2,Omega./U_sheared.f(y2),shift_y2,Omega,k3,U_sheared,c0,delta,N_bl,L1_sheared.l1,L1_sheared.l2,-1);
    Pressure_Upstream_sheared = -q0*WallPressure_sheared.*c0.f(y2).^2./U_sheared.df(y2);
    Full_Pressure_sheared(j,:) = Pressure_Upstream_sheared;
    
    PlotRI(y2,[Full_Pressure_base(j,:);Full_Pressure_sheared(j,:)]);
    disp([num2str(j),' of ',num2str(numel(omega))]);
end