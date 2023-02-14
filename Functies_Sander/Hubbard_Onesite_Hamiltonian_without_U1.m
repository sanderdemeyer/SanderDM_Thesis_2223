function H = Hubbard_Onesite_Hamiltonian_without_U1(pspace, trivspace, U, mu)
    % arguments are pspace, trivspace, U
    tens_one_site = Tensor([pspace trivspace], [trivspace pspace]);

    tblocks = {-mu; reshape([0 0;0 U - 2*mu], 2, 1, 2); -mu};
    
    H = fill_tensor(tens_one_site, tblocks);
    H = tpermute(H, [2 1 4 3], [2 2]);
end