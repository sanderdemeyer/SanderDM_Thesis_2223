function H = Hubbard_delta_Hamiltonian_None_SU2(pspace, trivspace, delta)
    % arguments are pspace, trivspace, delta
    % Gives an extra energy equal to delta per electron.
    tens_one_site = Tensor([trivspace pspace], [pspace trivspace]);

    tblocks = {reshape([0 0; 0 2], 1, 2, 1, 2), 1};
    
    H = delta*fill_tensor(tens_one_site, tblocks);
end