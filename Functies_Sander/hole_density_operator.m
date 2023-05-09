function O = hole_density_operator()
    % This assumes symmetry None_SU2

    pspace = get_spaces_Hubbard_None_SU2();
    tens_one_site = Tensor(pspace, pspace);

    O = fill_tensor(tens_one_site, {reshape([1 0; 0 -1], 2, 2), 0});
end