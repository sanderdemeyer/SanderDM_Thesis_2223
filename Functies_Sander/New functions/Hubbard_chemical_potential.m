function H = Hubbard_chemical_potential(pspace, trivspace, mu)
    % This assumes symmetry None_SU2
    
    tens_one_site = Tensor([trivspace pspace ], [pspace trivspace]);

    tblocks = {reshape([0 0; 0 2*mu], 2, 1, 2), mu};
    
    H = fill_tensor(tens_one_site, tblocks);
end