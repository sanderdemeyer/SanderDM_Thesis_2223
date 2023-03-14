function H = Hubbard_delta_Hamiltonian(pspace, trivspace, delta)
    % arguments are pspace, trivspace, delta
    % Gives an extra energy equal to delta per electron.
    tens_one_site = Tensor([trivspace pspace], [pspace trivspace]);

    tblocks = num2cell([0 1 2]*delta);
    
    H = fill_tensor(tens_one_site, tblocks);
end