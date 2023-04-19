function H = Hubbard_Onesite_Hamiltonian_None_SU2(pspace, trivspace, U)
    % arguments are pspace, trivspace, U
    tens_one_site = Tensor([pspace trivspace], [trivspace pspace]);

    tblocks = {reshape([0 0; 0 U], 2, 1, 2), 0};
    
    H = fill_tensor(tens_one_site, tblocks);
    H = tpermute(H, [2 1 4 3], [2 2]);
end