% Use underlying code to interpret the XXX model.

N = 0;
delta = 1;
stag_h_field = 0;
h_field = 0;

charges = U1([1 -1]);
pspace = GradedSpace.new(charges, [1 1], false);
trivspace = GradedSpace.new(U1(0), 1, false);

H = get_hamiltonian('XXZ', delta, pspace);
H = tpermute(H, [3 4 1 2], [2 2]);

%%
length = 4;
bond_dimensions = zeros(1,length);
deltas = zeros(1,length);
energies = zeros(1,length);
corr_lengths = zeros(1,length);
epsilon1s = zeros(1,length);
%name0 = 'Hubbard_1D_half_filling_idmrg2_t_1_U_8_cut_';
%names = {'4_final', '5_final', '7'};
%names = {'3_4_0_final.matVUMPS_final', '3_5_0_final.matVUMPS_final', '3_6_0.matVUMPS_final', '3_7_0.matVUMPS_final'};
names = {'3_4_0_final.mat', '3_5_0_final.mat', '3_6_0.mat', '3_7_0.mat'};
finals = [1 1 0 0];
%names = {'4_final', '5_final', '6', '7'};
for i = 1:length
    disp(i);
    %name = 'Hubbard_1D_half_filling_VUMPS_' + string(i+3);
    %name = 'Hubbard_1D_half_filling_idmrg2_t_0.625_U_10.8_cut_5_final';
    %name = strcat(name0, names{i});
    name = strcat('XXZ_Cylinder_', string(names{i}));
    %name = strcat('Hubbard_1D_half_filling_idmrg2_t_1_U_6_cut_', names{i});
    disp(name);
    load(name);
    if ~finals(i)
        gs_mps = canonicalize(mps, 'Order', 'rl');
    end
    %plot_entanglementspectrum(gs_mps);
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
scatter(deltas,epsilon1s, "red")
title('Extrapolation for XXX Heisenberg model on a helix with N = 3 with zero magnetic field.', 'interpreter', 'latex');
%xlabel('$\delta$', 'interpreter', 'latex');
%ylabel('$\epsilon$ = $\xsi$', 'interpreter', 'latex');
hold on
x = linspace(0,max(deltas)*1.1,500);
y = p1(2) + p1(1)*x;
plot(x, y, "red");
hold on
pol = arrayfun(@(x) polynomial(x, a, b), deltas);
RMS = sqrt(mean((pol-epsilon1s).^2));

scatter(deltas_new, epsilon1s_new, "blue")
hold on


p1_new = polyfit(deltas_new, epsilon1s_new, 1);
y_new = p1_new(2) + p1_new(1)*x;
plot(x,y_new, "blue");

p1_new_spec = polyfit(deltas_new(2:4), epsilon1s_new(2:4), 1);
y_new_spec = p1_new_spec(2) + p1_new_spec(1)*x;
plot(x, y_new_spec, "black");

legend("After IDmrg2: Data points", "After IDmrg2: Extrapolations", "After IDmrg2 and VUMPS: Data points", "After IDmrg2 and VUMPS: Extrapolation", "After IDmrg2 and VUMPS: Extrapolation with 3 best points");
xlabel('$\delta$', 'interpreter', 'latex');
ylabel('$\epsilon_1 = \xi^{-1}$', 'interpreter', 'latex');
hold off
function y = polynomial(x, a, b)
    y = a + b*x;
end
