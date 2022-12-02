for cut = [4 5 7]
    disp(cut);
    if cut == 4
        name = 'Hubbard_1D_half_filling_idmrg2_t_1_U_8_cut_4_final.mat';
    elseif cut == 5
        name = 'Hubbard_1D_half_filling_idmrg2_t_1_U_8_cut_5_final.mat';
    elseif cut == 7
        name = 'Hubbard_1D_half_filling_idmrg2_t_1_U_8_cut_7.mat';
    end
    Hubbard_1D_external_field_half_filling_VUMPS(1, 8, cut, name)
end