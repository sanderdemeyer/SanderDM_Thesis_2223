load('Hubbard_1D_half_filling_VUMPS_6.mat');

%%
mps = canonicalize(mps, 'Order', 'rl');

[pspace, vspace, trivspace, fusion_trees] = get_spaces('Hubbard', false, 1, 1, 12, 3);

t = 1.5;
U = 9;
mu = 0;
h_field = 0;

H = get_hamiltonian('Hubbard_external_field', fusion_trees, pspace, t, U, mu, h_field);
H = tpermute(H, [3 4 1 2], [2 2]);

AC1 = mps.AC(1);
AC2 = mps.AC(2);
AL1 = mps.AL(1);
AL2 = mps.AL(2);

E1 = contract(AL1, [1 2 3], AC2, [3 4 5], conj(AL1), [1 6 7], conj(twist(AC2,3)), [7 8 5], H, [2 4 6 8]);
E2 = contract(AL2, [1 2 3], AC1, [3 4 5], conj(AL2), [1 6 7], conj(twist(AC1,3)), [7 8 5], H, [2 4 6 8]);
E = (E1+E2)/2;

H_squared_1 = contract(AL1, [1 2 3], AC2, [3 4 5], conj(AL1), [1 6 7], conj(twist(AC2,3)), [7 8 5], H, [2 4 9 10], H, [9 10 6 8]);
H_squared_2 = contract(AL2, [1 2 3], AC1, [3 4 5], conj(AL2), [1 6 7], conj(twist(AC1,3)), [7 8 5], H, [2 4 9 10], H, [9 10 6 8]);
H_squared = (H_squared_1 + H_squared_2)/2;

error = sqrt(H_squared - E^2)/abs(E);
%get_filling(mps, pspace, trivspace, 2, false);

plot_entanglementspectrum(mps);
%%
max_dist = 20;
number_corr = get_number_correlation(mps, pspace, trivspace, fusion_trees, 2, 1/2, true, max_dist);

%%
%{
values = zeros(1,5);
for i = 1:20
    values(i) = norm(D(i,i));
end
disp(values);
%}

plot(number_corr)
xlabel('Number of sites between i and j')
ylabel('<n(i,up) n(j,down)> - <n(i,up)> <n(j,down)>')
title('Correlation')
