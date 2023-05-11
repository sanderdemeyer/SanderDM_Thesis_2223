function H = get_Hubbard_JMpo_oneband(t, t2, U, V, kwargs)
    arguments
        t
        t2
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
    assert(strcmp(kwargs.system{1}, 'Cylinder_multiple_rungs'), 'TBA system');

    if strcmp(kwargs.symmetries, 'U1_SU2')
        [pspace, ~, trivspace] = get_spaces_Hubbard_SU2(kwargs.P, kwargs.Q, 'D', kwargs.D);
        Hopping_t = Hubbard_Hopping_Hamiltonian_SU2(t, kwargs.P, kwargs.Q, 'convention', kwargs.convention);
        Hopping_t2 = Hubbard_Hopping_Hamiltonian_SU2(t2, kwargs.P, kwargs.Q, 'convention', kwargs.convention);
        V_term = Hubbard_V_Hamiltonian(pspace, trivspace, V, 'convention', kwargs.convention);
        U_term = Hubbard_Onesite_Hamiltonian_SU2(pspace, trivspace, U);
    elseif strcmp(kwargs.symmetries, 'None_SU2')
        [pspace, ~, trivspace] = get_spaces_Hubbard_None_SU2('D1', kwargs.D, 'D2', kwargs.D);
        Hopping_t = Hubbard_Hopping_Hamiltonian_None_SU2(t, 'convention', kwargs.convention);
        Hopping_t2 = Hubbard_Hopping_Hamiltonian_None_SU2(t2, 'convention', kwargs.convention);
        V_term = Hubbard_V_Hamiltonian_None_SU2(pspace, trivspace, V, 'convention', kwargs.convention);
        U_term_base = Hubbard_Onesite_Hamiltonian_None_SU2(pspace, trivspace, U);
        H_mu_term = Hubbard_chemical_potential(pspace, trivspace, kwargs.mu);
        U_term = U_term_base + H_mu_term;
    else
        error('Invalid symmetry')
    end

    if strcmp(kwargs.system{1}, 'Cylinder_multiple_rungs')
        mpo = get_mpo_cylinder_oneband(U_term, Hopping_t + V_term, Hopping_t2, kwargs.system{2}, kwargs.system{3});
    else
        error('TBA system');
    end

    H = InfJMpo(mpo);

end