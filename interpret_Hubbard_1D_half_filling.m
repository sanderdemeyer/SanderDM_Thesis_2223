
[pspace, vspace, trivspace, fusion_trees] = get_spaces('Hubbard', false, 1, 1, 12, 3);

%%

% For atomic separation = 2å of 1D hydrogen chain: 
% t = 1.5 en abs(U/t) = 6
t = 1.5;
U = 9;
mu = 0;
h_field = 0;

H = get_hamiltonian('Hubbard_external_field', fusion_trees, pspace, t, U, mu, h_field);
H = tpermute(H, [3 4 1 2], [2 2]);
mpo = get_mpo(H, 0);

H1 = InfJMpo({mpo mpo});

%mps = UniformMps.randnc(pspace, vspace);

disp('initialization correct');
cuts = [4 5 6];


%%
length = 4;
bond_dimensions = zeros(1,length);
deltas = zeros(1,length);
energies = zeros(1,length);
corr_lengths = zeros(1,length);
epsilon1s = zeros(1,length);
%name0 = 'Hubbard_1D_half_filling_idmrg2_t_1_U_8_cut_';
%names = {'4_final', '5_final', '7'};
cuts = [4 5 7];
%names = {'4_final', '5_final', '6', '7'};
for i = 1:length
    disp(i);
    %name = 'Hubbard_1D_half_filling_VUMPS_' + string(i+3);
    %name = 'Hubbard_1D_half_filling_idmrg2_t_0.625_U_10.8_cut_5_final';
    %name = strcat(name0, names{i});
    name = 'Hubbard_1D_half_filling_VUMPS_t_1_U_8_cut_' + string(cuts(i)) + '_final.mat';
    %name = strcat('Hubbard_1D_half_filling_idmrg2_t_1_U_6_cut_', names{i});
    disp(name);
    load(name);
    %gs_mps = canonicalize(mps, 'Order', 'rl');
    plot_entanglementspectrum(gs_mps);
    AC1 = gs_mps.AC(1);
    AC2 = gs_mps.AC(2);
    AL1 = gs_mps.AL(1);
    AL2 = gs_mps.AL(2);
    disp(AL1.dims);
    bond_dimensions(i) = max(AL1.dims);
    E1 = contract(AL1, [1 2 3], AC2, [3 4 5], conj(AL1), [1 6 7], conj(twist(AC2,3)), [7 8 5], H, [2 4 6 8]);
    E2 = contract(AL2, [1 2 3], AC1, [3 4 5], conj(AL2), [1 6 7], conj(twist(AC1,3)), [7 8 5], H, [2 4 6 8]);
    E = (E1+E2)/2;
    disp(E);
    energies(i) = E;
    [V, D] = transfereigs(gs_mps, gs_mps, 5);
    epsilons = zeros(1,5);
    for j = 1:5
        epsilons(j) = -log(norm(D(j,j)));
    end
    disp(epsilons);
    epsilon1s(i) = epsilons(2);
    corr_lengths(i) = 1/epsilons(2);
    delta = epsilons(3)-epsilons(2);
    deltas(i) = delta;
end

%%

p1 = polyfit(deltas, epsilon1s, 1);
a = p1(2);
b = p1(1);
scatter(deltas,epsilon1s)
title('Extrapolation', 'interpreter', 'latex');
%xlabel('$\delta$', 'interpreter', 'latex');
%ylabel('$\epsilon$ = $\xsi$', 'interpreter', 'latex');
hold on
x = linspace(0,0.02,500);
y = p1(2) + p1(1)*x;
plot(x, y);
hold off
pol = arrayfun(@(x) polynomial(x, a, b), deltas);
RMS = sqrt(mean((pol-epsilon1s).^2));

function y = polynomial(x, a, b)
    y = a + b*x;
end
