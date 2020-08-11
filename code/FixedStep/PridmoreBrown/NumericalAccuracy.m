k3 = 0;
omega = 1;
delta = 1;

c = Struct_c_Const(2);
U = Struct_U_Parabolic(0.5,1,delta);

k1 = [-1, 0, 1i, 3, 10, 1.5+0.1i, 1.5+0.05i, 1.5+0.01i]; % This is the most changeable parameter.

%N_step = 10;
N_count = 100;
N_logmin = 0.3;
N_max = 10000 ;
N_logmax = 4;

%N_extra = 10*N_max;
N_extra = 10^(N_logmax + 1);

%N = N_step:N_step:N_max;
N = logspace(N_logmin,N_logmax,N_count);
N = round(N);
D = zeros(numel(k1),numel(N));

D_max = zeros(numel(k1),1);

h = waitbar(0,'Recalculating masses');

for i = 1:numel(k1)
    for j = 1:numel(N)
        [dec] = Dec_Aug18A(k1(i),k3,omega,U,c,delta,N(j));
        D(i,j) = dec.v(1); % i.e. dispersion relation, evaluated on the wall
        waitbar(j./numel(N),h,[num2str(round(j./numel(N)*100,1)),'% complete']);
    end
    [dec] = Dec_Aug18A(k1(i),k3,omega,U,c,delta,N_extra);
    D_max(i) = dec.v(1);
    disp([num2str(i),' of ',num2str(numel(k1)),' k1 evaluated']);
end
close(h)

l = loglog(N,abs(D-D_max));
for j = 1:numel(l)
    l(j).LineWidth = 2;
end

xlabel('Integration points, $N = 1/(\delta x_2)$','Interpreter','latex');
ylabel('Relative error, $|D_h(k_1,N) - D_h(k_1,N_\infty)|$','Interpreter','latex');

m = cell(size(k1));
for j = 1:numel(m)
    m{j} = ['$k1 = ',num2str(k1(j)),'$'];
end

leg = legend(m,'Interpreter','latex','location','northeast');