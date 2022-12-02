N = 0;
delta = 2;
stag_h_field = 0;
shear = 0;

charges = U1([1 -1]);
pspace = GradedSpace.new(charges, [1 1], false);
trivspace = GradedSpace.new(U1(0), 1, false);

%%
H = get_hamiltonian('XXZ', delta, pspace);
H = tpermute(H, [3 4 1 2], [2 2]);
H_one_site = get_hamiltonian('one_site_XXZ', pspace, trivspace, stag_h_field, shear);
H_one_site_A = H_one_site{1};
H_one_site_B = H_one_site{2};

mpo = get_mpo(H, N);
mpo_joint = {mpo mpo};

mpo_joint{1}(1, 1, N+3, 1) = tpermute(H_one_site_A, [2 1 4 3], [2 2]);
mpo_joint{2}(1, 1, N+3, 1) = tpermute(H_one_site_B, [2 1 4 3], [2 2]);

H1 = InfJMpo(mpo_joint);
disp('Initialization done');
%%
%errors = zeros(1,6);
deltas = zeros(1,5);
energies = zeros(1,5);
corr_lengths = zeros(1,5);
epsilon1s = zeros(1,5);
for cut = 4:8
    load("XXZ_1D_2_VUMPS_" + string(cut) + "_final.mat");
    disp(cut);
    disp(gs_mps.AL);

    %plot_entanglementspectrum(gs_mps);
    AC1 = gs_mps.AC(1);
    AC2 = gs_mps.AC(2);
    AL1 = gs_mps.AL(1);
    AL2 = gs_mps.AL(2);
    
    E1 = contract(AL1, [1 2 3], AC2, [3 4 5], conj(AL1), [1 6 7], conj(twist(AC2,3)), [7 8 5], H, [2 4 6 8]);
    E2 = contract(AL2, [1 2 3], AC1, [3 4 5], conj(AL2), [1 6 7], conj(twist(AC1,3)), [7 8 5], H, [2 4 6 8]);
    E = (E1+E2)/2;
    disp(E);
    energies(cut-3) = E;
    %{
    H_squared_1 = contract(AL1, [1 2 3], AC2, [3 4 5], conj(AL1), [1 6 7], conj(twist(AC2,3)), [7 8 5], H, [2 4 9 10], H, [9 10 6 8]);
    H_squared_2 = contract(AL2, [1 2 3], AC1, [3 4 5], conj(AL2), [1 6 7], conj(twist(AC1,3)), [7 8 5], H, [2 4 9 10], H, [9 10 6 8]);
    H_squared = (H_squared_1 + H_squared_2)/2;
    
    error = sqrt(H_squared - E^2)/abs(E);
    %get_filling(mps, pspace, trivspace, 2, false);
    errors(power+1) = error;
    %}
    %plot_entanglementspectrum(mps);
    [V, D] = transfereigs(gs_mps, gs_mps, 5);

    epsilons = zeros(1,5);
    for i = 1:5
        epsilons(i) = -log(norm(D(i,i)));
    end
    disp(epsilons);
    delta = epsilons(8)-epsilons(6);
    deltas(cut-3) = delta;
    epsilon1s(cut-3) = epsilons(2);
    corr_lengths(cut-3) = -1/(norm(D(2,2)));
end

disp(deltas);
disp(energies);
disp(corr_lengths);

%%
plot(errors, energies);
xlabel('epsilon spectrum');
ylabel('energy of the calculation')
title('Extrapolation for the energy of the Heisenberg XXZ model')
%%

plot(deltas)
xlabel('cut = 1e(-x-2)')
ylabel('delta')
%%

scatter(deltas, epsilon1s);
hold on
x = linspace(0, 0.4, 100);
plot(x, 0.3798 + 0.5380*x);
xlabel('delta');
ylabel('epsilon1 = 1/correlation length');
hold off