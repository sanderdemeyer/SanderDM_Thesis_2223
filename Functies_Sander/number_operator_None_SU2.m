function H = number_operator_None_SU2(pspace)
    % This assumes symmetry None_SU2
    
    tens_one_site = Tensor(pspace, pspace);

    tblocks = {reshape([0 0; 0 2], 2, 2), 1};
    
    H = fill_tensor(tens_one_site, tblocks);
end