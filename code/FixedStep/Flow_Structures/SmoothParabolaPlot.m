% Smooth parabolic profile plot

U0 = 0.5:-0.1:0;

N = numel(U0);
U = cell(1,N);

for j = 1:N
    U{j} = Struct_U_Parabolic(U0(j),1,1);
end

x2a = 0:0.01:1;
x2b = 1.01:0.01:2;

Ua = zeros([numel(x2a),N]);
Ub = zeros([numel(x2b),N]);

for j = 1:N
    Ua(:,j) = U{j}.f(x2a);
    Ub(:,j) = U{j}.f(1).*ones(size(x2b));
end

U_comb = [Ua;Ub];


h = plot(U_comb,[x2a,x2b]);

for j = 1:numel(h)
    h(j).LineWidth = 2;
end

col = get(gca,'colororder');
for j = 2:numel(h)
    h(j).LineStyle = '--';
    h(j).Color = col(1,:);
end

axis([ 0 1.2 0 2]);

xlabel('$U(x_2)$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');