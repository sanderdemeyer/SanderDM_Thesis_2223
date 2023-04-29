function kinetic_energy = get_kinetic_energy(gs_mps, kwargs)
    arguments
        gs_mps
        kwargs.t = 1
        kwargs.P = 1
        kwargs.Q = 1
        kwargs.D = 1
        kwargs.symmetries = 'U1_SU2'
        kwargs.convention = 'conventional'
    end
    if strcmp(kwargs.symmetries, 'U1_SU2')
        [pspace, ~, trivspace] = get_spaces_Hubbard_SU2(kwargs.P, kwargs.Q, 'D', kwargs.D);
        Hopping_t = Hubbard_Hopping_Hamiltonian_SU2(kwargs.t, kwargs.P, kwargs.Q, 'convention', kwargs.convention);
    elseif strcmp(kwargs.symmetries, 'None_SU2')
        [pspace, ~, trivspace] = get_spaces_Hubbard_None_SU2('D1', kwargs.D, 'D2', kwargs.D);
        Hopping_t = Hubbard_Hopping_Hamiltonian_None_SU2(kwargs.t, 'convention', kwargs.convention);
    else
        error('Symmetry %s not implemented!', kwargs.symmetries)
    end

    AC = gs_mps.AC;
    AR = gs_mps.AR;
    w = period(gs_mps);
    kinetic_energies = cell(1, w);
    for b = 1:w
        kinetic_energies{b} = contract(AC(b), [1 2 3], AR(loop(b,1,w)), [3 4 5], conj(AC(b)), [1 6 7], twist(conj(AR(loop(b,1,w))),3), [7 8 5], Hopping_t, [6 8 4 2]);
    end
    kinetic_energy = mean(cell2mat(kinetic_energies));
end