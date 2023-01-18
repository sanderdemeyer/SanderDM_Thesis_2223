function H = get_Hubbard_JMpo(t, U, kwargs)
    arguments
        t
        U
        kwargs.P = 1;
        kwargs.Q = 1;
        kwargs.system = {'1D'}
        kwargs.D = 1
    end

    if kwargs.P == 1
        [pspace, ~, trivspace] = get_spaces_Hubbard_symmetric(kwargs.P, kwargs.Q, 'D1', kwargs.D, 'D2', kwargs.D);
    else
        [pspace, ~, trivspace] = get_spaces_Hubbard_asymmetric(kwargs.P, kwargs.Q, 'D', kwargs.D);
    end

    H_twosite = Hubbard_Hopping_Hamiltonian(t, kwargs.P, kwargs.Q);
    H_onesite = get_hamiltonian('Hubbard_one_site', pspace, trivspace, U);
    
    if strcmp(kwargs.system{1}, '1D')
        mpo = get_mpo_1D(H_twosite, H_onesite);
    elseif strcmp(kwargs.system{1}, 'Helix')
        mpo = get_mpo_helix(H_twosite, H_onesite, kwargs.system{2});
    elseif strcmp(kwargs.system{1}, 'Cylinder')
        mpo = get_mpo_cylinder(H_twosite, H_onesite, kwargs.system{2}, 1);
    elseif strcmp(kwargs.system{1}, 'DoubleCylinder')
        mpo = get_mpo_cylinder(H_twosite, H_onesite, kwargs.system{2}, 2);
    else
        error('System %s not implemented', kwargs.system{1});
    end

    H = InfJMpo(mpo);
end