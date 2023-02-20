function H = Hubbard_Onesite_Hamiltonian_SU2(pspace, trivspace, U)
    % arguments are pspace, trivspace, U
    tens_one_site = Tensor([pspace trivspace], [trivspace pspace]);

    tblocks = num2cell([0 0 1]*U);
    
    H = fill_tensor(tens_one_site, tblocks);
    H = tpermute(H, [2 1 4 3], [2 2]);
end