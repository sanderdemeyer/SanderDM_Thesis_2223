function do_vanalles(U)
    gs_mps = Hubbard_1D_external_field_half_filling(1, U, 7, 0);
    for cut = 6:-1:4
        gs_mps = Hubbard_1D_external_field_half_filling(1, U, cut, 1, gs_mps);
    end
end