folder = 'Data structures/Hubbard_t_1_24_dec/';
folder = 'Data structures/Hubbard_1D_U_6_new/';
%folder = 'Data structures/Hubbard_1D_U_8_doping/';

files = dir((folder));
l = length(files);
doPath;
disp('started');
%{
deltas = zeros(1, l-2);
epsilon1s = zeros(1, l-2);
variances = zeros(1, l-2);
bond_dims = zeros(1, l-2);
energies = zeros(1, l-2);
stag_magn = zeros(1, l-2);
inv_corr = zeros(1, l-2);
Us = zeros(1, l-2);
%}

k = 1;
for i = 3:l
    file = files(i);
    naam = file.name;
    name_l = length(naam);
    if name_l == 1
        naam = naam{1};
        name_l = length(naam);
    end
    disp(naam);
    U = str2double(naam(31:length(naam)-38));
    disp(U);
    disp(length(naam));
    if strcmp(naam(name_l-9:name_l), '_final.mat')
%        gs_mps = canonicalize(mps, 'Order', 'rl');
%    end
%    if 0 == 0
        load(strcat(folder, naam))%, 'gs_mps', 'gs_energy');
        disp(naam);
        disp(naam(length(naam)-8:length(naam)));
        %load(name);
        dimensions = dims(gs_mps.AL(1).var);
        bond_dim = dimensions(1) + dimensions(3);
        bond_dims{k} = bond_dim;
        energies{k} = gs_energy;
        disp(U);
        Us{k} = U;
        etas{k} = eta;

        %{
        filling{k} = P/Q;
        if mod(P,2) == 0
            unit_cell{k} = Q;
        else
            unit_cell{k} = 2*Q;
        end
        U = 8;
        H_one_site = get_hamiltonian('Hubbard_one_site', pspace, trivspace, U);
        single_site_energies = zeros(1, unit_cell{k});
        for site = 1:unit_cell{k}
            AC = gs_mps.AC(site);
            single_site_energies(site) = contract(AC, [1 2 3], H_one_site, [-1 4 -2 2], twist(conj(AC),3), [1 4 3]).var.var;
        end
        self_interaction_energies{k} = mean(single_site_energies);
        %}

        [V, D] = transfereigs(gs_mps, gs_mps, 3);
        epsilons = zeros(1,3);
        for j = 1:3
            epsilons(j) = -log(norm(D(j,j)));
        end
        all_sectors{k} = D;
        deltas{k} = epsilons(3) - epsilons(2);
        epsilon1s{k} = epsilons(2);
        disp('0 sector done');
        %{
        [V, D] = transfereigs(gs_mps, gs_mps, 3, 'Charge', ProductCharge(U1(2*Q), U1(0), fZ2(0)));
        epsilons = zeros(1,3);
        for j = 1:3
            epsilons(j) = -log(norm(D(j,j)));
        end
        disp(epsilons(1));
        charge_sector{k} = D;
        deltas_charge{k} = epsilons(2) - epsilons(1);
        epsilon1s_charge{k} = epsilons(1);
        disp('charge sector done');
        %}
        %{
        [V, D] = transfereigs(gs_mps, gs_mps, 3, 'Charge', ProductCharge(U1(0), U1(2), fZ2(0)));
        epsilons = zeros(1,3);
        for j = 1:3
            epsilons(j) = -log(norm(D(j,j)));
        end
        spin_sector{k} = D;
        deltas_spin{k} = epsilons(2) - epsilons(1);
        epsilon1s_spin{k} = epsilons(1);
        disp('Spin sector done');
        %}

        %{
        stag_m = get_magnetisation('XXZ', gs_mps, pspace, trivspace, 1.5, false, true);
        fprintf('Staggered magnetisation is %s \n', stag_m)
        stag_magn(k) = stag_m;
        %}
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
x_data = cell2mat(deltas_charge);
y_data = cell2mat(epsilon1s_charge);

resolution = 500;
scatter(x_data, y_data);

hold on
p1 = polyfit(x_data, y_data, 1);
x = linspace(0,max(x_data)*1.1,resolution);
y = p1(2) + p1(1)*x;

%plot(x, y, "red");
xlabel('$\delta = \epsilon_2 - \epsilon_1$', 'interpreter', 'latex');
ylabel('$\epsilon_1 = 1/\xi$', 'interpreter', 'latex')
title('Extrapolation of $\epsilon_1$ for the 1D Hubbard model - Charge sector', 'interpreter', 'latex')
hold off

pol = arrayfun(@(x) x*p1(1) + p1(2), x_data);
RMS1 = sqrt(mean((pol-y_data).^2));

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


p3 = polyfit(deltas(1:11), stag_magn(1:11), 1);


scatter(bond_dims(1:11)./2, abs(stag_magn(1:11)));
xlabel('Bond dimension')
ylabel('spontaneous staggered magnetisation')
title('Staggered magnetisation as finite bond dimension effect')



%%

scatter(bond_dims./2, abs(stag_magn));
xlabel('Bond dimension')
ylabel('spontaneous staggered magnetisation')
title('Helix, N = 5')

%%

scatter(Us_array./(Us_array-4), cell2mat(energies));
xlim([0,10]);

%%
plot(x, y);
hold on
scatter(1./bond_dims, abs(stag_magn));
hold off
xlabel('1/D', 'interpreter', 'latex');
ylabel('spontaneous staggered magnetisation');
title('Staggered magnetisation as finite bond dimension effect');

%%

scatter(Us_a, cell2mat(energies)/2);
xlim([0 10]);

%%

scatter(Us_a./(Us_a+4), cell2mat(energies)/2);
%%
function y = polynomial(x, a, b)
    y = a*x+b;
end

