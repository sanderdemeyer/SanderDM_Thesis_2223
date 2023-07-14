%[gs_mps, gs_energy] = Heisenberg_XXZ_1D(1, 4, -2, [20 20 500], 2.75, 2.75 + 1.5, 0, 0);
for trunc = [3 3.25 3.5 3.75 4 4.25 4.5 4.75 5 5.25 5.5 5.75 6 6.25 6.5 6.75 7]
    prev_name = 'Heisenberg_XXZ_1D_delta_1_trunctotdim_' + string(trunc-0.25) + '_final.mat';
    [gs_mps, gs_energy] = Heisenberg_XXZ_1D(1, 4, -2, [20 20 5.0], trunc, trunc + 1.5, prev_name, 2);
end