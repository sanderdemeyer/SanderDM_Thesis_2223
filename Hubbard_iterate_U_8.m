disp('Started running');
Hubbard_1D_external_field_half_filling(1, 6, 300, [75 5 2], 4, -1, 0, 0, 0, 0)

for D = [400 500 600 700 800]
    fprintf('Started for D = %d \n', D);
    name = 'Hubbard_1D_half_filling_t_1_U_6_trunctotdim_' + string(D-100) + '_h_0_redef_0_final.mat';
    Hubbard_1D_external_field_half_filling(1, 8, D, [75 5 2], 4, -1, 0, 0, name, 2)
end
