files = dir(('Data structures/Prachtig'));
%files = dir(('Data structures/XXZ_1D_delta_1_sizeD'));
l = length(files);
deltas = zeros(1, 5);
epsilon1s = zeros(1, 5);
variances = zeros(1, 5);
bond_dims = zeros(1, 5);
energies = zeros(1, 5);
k = 1;
for i = 3:l
    file = files(i);
    name = file.name;
    name_l = length(name);
    if strcmp(name(name_l-9:name_l), '_final.mat')
        disp(name);
        load(name);
        %dimensions = dims(gs_mps.AL(1));
        %bond_dim = dimensions(1) + dimensions(3);
        %bond_dims(k) = bond_dim;
        [V, D] = transfereigs(gs_mps, gs_mps, 10);
        epsilons = zeros(1,10);
        for j = 1:10
            epsilons(j) = -log(norm(D(j,j)));
        end
        disp(epsilons);
        corr_length = 1/epsilons(2);
        deltas(k) = epsilons(3) - epsilons(2);
        epsilon1s(k) = epsilons(2);
        %{
        AL1 = gs_mps.AL(1);
        AC2 = gs_mps.AC(2);
        W1 = mpo_joint{1};
        W2 = mpo_joint{2};
        [GL, GR] = environments(H1, gs_mps, gs_mps);
        [P, ~, Q] = size(W1);

        %exp = expectation_value(gs_mps, H1, gs_mps);
        top_str = contract(AL1, [-1 -2 2], AC2, [2 -3 -4]);
        hoera = contract(top_str, [1 3 5 -5], GL{1}, [-1 2 1], W1, [2 -2 4 3], W2(:,:,Q,:), [4 -3 -4 5]);
        exp = contract(hoera, [1 2 3 -1 4], conj(AL1), [1 2 5], conj(AC2), [5 3 4]);
        disp(exp);
        energies(k) = exp.var.var;

        top_str = contract(AL1, [-1 1 2], AC2, [2 3 -4], H, [1 3 -2 -3]);
        hoera = contract(top_str, [1 3 5 -5], GL{1}, [-1 2 1], W1, [2 -2 4 3], W2(:,:,Q,:), [4 -3 -4 5]);
        var_A = contract(hoera, [1 2 3 -1 4], conj(AL1), [1 2 5], conj(AC2), [5 3 4]);
        disp(var_A.var.var);
        variances(k) = var_A.var.var;
        %}

        k = k+1;
    end
end

%%

scatter(deltas, epsilon1s);
hold on
p1 = polyfit(deltas, epsilon1s, 1);
x = linspace(0,max(deltas)*1.1,500);
y = p1(2) + p1(1)*x;
plot(x, y, "red");
xlabel('$\delta = \epsilon_2 - \epsilon_1$', 'interpreter', 'latex');
ylabel('$\epsilon_1 = 1/\xi$', 'interpreter', 'latex')
title('Extrapolation of $\epsilon_1$ for the Heisenberg XXX model', 'interpreter', 'latex')
hold off

pol = arrayfun(@(x) polynomial(x, p1(1), p1(2)), deltas);
RMS1 = sqrt(mean((pol-epsilon1s).^2));

%%
rec_bond_dims = 1./bond_dims;
scatter(rec_bond_dims, epsilon1s);
hold on
p2 = polyfit(rec_bond_dims, epsilon1s, 1);
x = linspace(0,max(rec_bond_dims)*1.1,500);
y = p2(2) + p2(1)*x;
plot(x, y, "red");

xlabel('$1/D$', 'interpreter', 'latex');
ylabel('$\epsilon_1 = 1/\xi$', 'interpreter', 'latex')
title('Extrapolation of $\epsilon_1$ for the Heisenberg XXX model', 'interpreter', 'latex')
hold off

pol = arrayfun(@(x) polynomial(x, p2(1), p2(2)), deltas);
RMS2 = sqrt(mean((pol-epsilon1s).^2));

%%
p3 = polyfit(variances, epsilon1s, 1);

function y = polynomial(x, a, b)
    y = a*x+b;
end


