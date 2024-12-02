%% Computes the wall pressure due to a point vortex at some location y2 above a hard wall

% Base case
U_points_base = [0.6,0.7,1];
delta_points_base = [0,0.1,1];
U_base = Struct_U_PWParabolic(U_points_base,delta_points_base);

% Sheared case
U_points_sheared = [0.2,0.75,1];
delta_points_sheared = [0,0.3,1];
U_sheared = Struct_U_PWParabolic(U_points_sheared,delta_points_sheared);


y2 = linspace(0, 1, 40);
shift_y2 = 0.005;

delta = delta_points_base(end); % should hopefully be same both base and sheared
Uinf = U_base.f(delta);
cinf = 35;

omega = logspace(-1,3,40);
f = 10*omega./(2*pi*0.02);

% Various computational parameters
N_bl = 3000;
k3 = 0;
c0 = Struct_c_Const(cinf);
q0 = 1;

FF_Pressure_base = zeros(numel(omega),numel(y2));
FF_Pressure_sheared = zeros(numel(omega),numel(y2));

for j = 1:numel(omega)
    Omega = omega(j);
    Y2 = y2;
    U = U_base;
    WallPressure_PointSource_FarField;
    FF_Pressure_base(j,:) = FF1p.p(:);
%     PlotRI(y2,FF_Pressure_base(j,:));
    U = U_sheared;
    WallPressure_PointSource_FarField;
    FF_Pressure_sheared(j,:) = FF1p.p(:);
%     PlotRI(y2,FF_Pressure_sheared(j,:));
    disp([num2str(j),' of ',num2str(numel(omega))]);
end