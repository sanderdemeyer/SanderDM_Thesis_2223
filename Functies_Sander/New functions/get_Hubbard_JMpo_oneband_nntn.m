function H = get_Hubbard_JMpo_oneband_nntn(t, t2, t3, U, V, kwargs)
    arguments
        t
        t2
        t3
        U
        V
        kwargs.P = 1
        kwargs.Q = 1
        kwargs.system = {'1D'}
        kwargs.D = 1
        kwargs.len = 0
        kwargs.symmetries = 'U1_U1'
        kwargs.mu = 0
        kwargs.convention = 'conventional'
    end
    assert(strcmp(kwargs.symmetries, 'U1_SU2'), 'TBA symmetry');
    assert(strcmp(kwargs.system{1}, 'Cylinder_multiple_rungs'), 'TBA system');

    if strcmp(kwargs.symmetries, 'U1_SU2')
        [pspace, ~, trivspace] = get_spaces_Hubbard_SU2(kwargs.P, kwargs.Q, 'D', kwargs.D);
        Hopping_t = Hubbard_Hopping_Hamiltonian_SU2(t, kwargs.P, kwargs.Q, 'convention', kwargs.convention);
        Hopping_t2 = Hubbard_Hopping_Hamiltonian_SU2(t2, kwargs.P, kwargs.Q, 'convention', kwargs.convention);
        Hopping_t3 = Hubbard_Hopping_Hamiltonian_SU2(t3, kwargs.P, kwargs.Q, 'convention', kwargs.convention);
        V_term = Hubbard_V_Hamiltonian(pspace, trivspace, V, 'convention', kwargs.convention);
        U_term = Hubbard_Onesite_Hamiltonian_SU2(pspace, trivspace, U);
    else
        error('Invalid symmetry')
    end    

    if strcmp(kwargs.system{1}, 'Cylinder_multiple_rungs')
        mpo = get_mpo_cylinder_oneband_nntn(U_term, Hopping_t + V_term, Hopping_t2, Hopping_t3, kwargs.system{2}, kwargs.system{3});
    else
        error('TBA system');
    end

    H = InfJMpo(mpo);

end