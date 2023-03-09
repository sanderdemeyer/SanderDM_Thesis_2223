function H = Hubbard_delta_Hamiltonian(pspace, trivspace, delta_pd)
    % To be placed only on the Cu atoms in the three-band model.
    arguments
        pspace
        trivspace
        delta_pd
    end
    warning('Not sure about this one');
    tens_one_site = Tensor([pspace trivspace], [trivspace pspace]);

    tblocks = num2cell(-delta_pd*[0 1 2]);
    
    H = fill_tensor(tens_one_site, tblocks);
    H = tpermute(H, [2 1 4 3], [2 2]);
end