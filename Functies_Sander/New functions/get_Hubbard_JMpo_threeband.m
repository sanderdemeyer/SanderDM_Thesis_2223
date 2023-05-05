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
    
    assert(strcmp(kwargs.system{1}, 'Cylinder_multiple_rungs'), 'TBA system %s', kwargs.system{1});
    if ~strcmp(kwargs.convention, 'conventional')
        warning('It is strongly recommended to work with convention = conventional.');
    end

    if strcmp(kwargs.symmetries, 'U1_SU2')
        % minus sign in Hubbard_Hopping_Hamiltonian_SU2 denote the fact
        % that those functions are defined with a negative sign, whereas
        % here the negative sign will be included in the construction of
        % the mpo. The definition of the various terms can be found in the
        % paper Ground-state phase diagram of the three-band Hubbard model 
        % from density matrix embedding theory - Zhi-Hao Cui 2020.

        [pspace, ~, trivspace] = get_spaces_Hubbard_SU2(kwargs.P, kwargs.Q, 'D', kwargs.D);
        Hopping_t_dp = Hubbard_Hopping_Hamiltonian_SU2(param.t_dp, kwargs.P, kwargs.Q, 'convention', kwargs.convention);
        Hopping_t_pp = Hubbard_Hopping_Hamiltonian_SU2(param.t_pp, kwargs.P, kwargs.Q, 'convention', kwargs.convention);

        %delta dp term will be put with a positive sign on the d-orbitals.
        delta_dp_term = Hubbard_delta_Hamiltonian(pspace, trivspace, param.delta_dp); 
        warning('check dp or pd!');
        V_pp_term = Hubbard_V_Hamiltonian(pspace, trivspace, param.V_pp, 'convention', kwargs.convention);
        V_dp_term = Hubbard_V_Hamiltonian(pspace, trivspace, param.V_dp, 'convention', kwargs.convention);
        V_dd_term = Hubbard_V_Hamiltonian(pspace, trivspace, param.V_dd, 'convention', kwargs.convention);
        U_dd_term = Hubbard_Onesite_Hamiltonian_SU2(pspace, trivspace, param.U_dd);
        U_pp_term = Hubbard_Onesite_Hamiltonian_SU2(pspace, trivspace, param.U_pp);

        % V interactions are always positive, hopping terms can differ in
        % sign in the 3-band model. Hence, pos is used whenever the hopping
        % term is positive, neg is used whenever the hopping term is
        % negative. This way, it is ensured the V-terms are always
        % positive.

        
    elseif strcmp(kwargs.symmetries, 'None_SU2')

        [pspace, ~, trivspace] = get_spaces_Hubbard_None_SU2('D1', kwargs.D, 'D2', kwargs.D);
        Hopping_t_dp = Hubbard_Hopping_Hamiltonian_None_SU2(param.t_dp, 'convention', kwargs.convention);
        Hopping_t_pp = Hubbard_Hopping_Hamiltonian_None_SU2(param.t_pp, 'convention', kwargs.convention);

        delta_dp_term = Hubbard_delta_Hamiltonian_None_SU2(pspace, trivspace, param.delta_dp); 
        V_pp_term = Hubbard_V_Hamiltonian_None_SU2(pspace, trivspace, param.V_pp, 'convention', kwargs.convention);
        V_dp_term = Hubbard_V_Hamiltonian_None_SU2(pspace, trivspace, param.V_dp, 'convention', kwargs.convention);
        V_dd_term = Hubbard_V_Hamiltonian_None_SU2(pspace, trivspace, param.V_dd, 'convention', kwargs.convention);
        U_dd_term = Hubbard_Onesite_Hamiltonian_None_SU2(pspace, trivspace, param.U_dd);
        U_pp_term = Hubbard_Onesite_Hamiltonian_None_SU2(pspace, trivspace, param.U_pp);
    else
        error('Invalid symmetry')
    end

    H1_Cu = U_dd_term + delta_dp_term;
    H1_O = U_pp_term;
    H2_OO_pos = Hopping_t_pp + V_pp_term;
    H2_OO_neg = -Hopping_t_pp + V_pp_term;
    H2_CuCu = V_dd_term;
    H2_CuO_pos = Hopping_t_dp + V_dp_term;
    H2_CuO_neg = -Hopping_t_dp + V_dp_term;

    if strcmp(kwargs.system{1}, 'Cylinder_multiple_rungs')
        mpo = get_mpo_cylinder_threeband(H1_Cu, H1_O, H2_CuCu, H2_CuO_pos, H2_CuO_neg, H2_OO_pos, H2_OO_neg, kwargs.system{2}, kwargs.system{3}, 'convention', kwargs.convention);
    else
        error('TBA system');
    end

    H = InfJMpo(mpo);

end