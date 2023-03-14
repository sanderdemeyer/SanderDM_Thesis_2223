function H = get_Hubbard_JMpo_threeband(param, kwargs)
    arguments
        param
        kwargs.P = 1
        kwargs.Q = 1
        kwargs.system = {'1D'}
        kwargs.D = 1
        kwargs.len = 0
        kwargs.symmetries = 'U1_U1'
        kwargs.mu = 0
        kwargs.convention = 'conventional'
    end
    
    assert(strcmp(kwargs.symmetries, 'U1_SU2'), 'TBA symmetry %s', kwargs.symmetries);
    assert(strcmp(kwargs.system{1}, 'Cylinder_multiple_rungs'), 'TBA system %s', kwargs.system{1});
    if ~strcmp(kwargs.convention, 'conventional')
        warning('It is strongly recommended to work with convention = conventional.');
    end

    if strcmp(kwargs.symmetries, 'U1_SU2')
        [pspace, ~, trivspace] = get_spaces_Hubbard_SU2(kwargs.P, kwargs.Q, 'D', kwargs.D);
        Hopping_t_dp = Hubbard_Hopping_Hamiltonian_SU2(param.t_dp, kwargs.P, kwargs.Q, 'convention', kwargs.convention);
        Hopping_t_pp = Hubbard_Hopping_Hamiltonian_SU2(param.t_pp, kwargs.P, kwargs.Q, 'convention', kwargs.convention);
        delta_dp_term = Hubbard_delta_Hamiltonian(pspace, trivspace, param.delta_dp); 
        warning('check dp or pd!');
        V_pp_term = Hubbard_V_Hamiltonian(pspace, trivspace, param.V_pp, 'convention', kwargs.convention);
        V_dp_term = Hubbard_V_Hamiltonian(pspace, trivspace, param.V_dp, 'convention', kwargs.convention);
        V_dd_term = Hubbard_V_Hamiltonian(pspace, trivspace, param.V_dd, 'convention', kwargs.convention);
        U_dd_term = Hubbard_Onesite_Hamiltonian_SU2(pspace, trivspace, param.U_dd);
        U_pp_term = Hubbard_Onesite_Hamiltonian_SU2(pspace, trivspace, param.U_pp);

        H1_Cu = U_dd_term + delta_dp_term;
        H1_O = U_pp_term;
        H2_OO = Hopping_t_pp + V_pp_term;
        H2_CuCu = V_dd_term;
        H2_CuO = Hopping_t_dp + V_dp_term;
    else
        error('Invalid symmetry')
    end

    if strcmp(kwargs.system{1}, 'Cylinder_multiple_rungs')
        mpo = get_mpo_cylinder_threeband(H1_Cu, H1_O, H2_CuCu, H2_CuO, H2_OO, kwargs.system{2}, kwargs.system{3});
    else
        error('TBA system');
    end

    H = InfJMpo(mpo);

end