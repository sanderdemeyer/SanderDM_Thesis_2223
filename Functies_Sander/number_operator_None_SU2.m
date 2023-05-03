function O = number_operator_None_SU2()
    % This assumes symmetry None_SU2

    pspace = get_spaces_Hubbard_None_SU2();
    tens_one_site = Tensor(pspace, pspace);

    O = fill_tensor(tens_one_site, {reshape([0 0; 0 2], 2, 2), 1});
end