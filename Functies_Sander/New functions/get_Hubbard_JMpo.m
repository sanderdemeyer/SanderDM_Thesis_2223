function H = get_Hubbard_JMpo(t, U, kwargs)
    arguments
        t
        U
        kwargs.t2 = 0
        kwargs.V = 0
        kwargs.P = 1
        kwargs.Q = 1
        kwargs.system = {'1D'}
        kwargs.D = 1
        kwargs.len = 0
        kwargs.symmetries = 'U1_U1'
        kwargs.mu = 0
        kwargs.convention = 'conventional'
        kwargs.bitstring = none
    end

    if strcmp(kwargs.symmetries, 'U1_U1')
        [pspace, ~, trivspace] = get_spaces_Hubbard_asymmetric(kwargs.P, kwargs.Q, 'D', kwargs.D);
        Hopping_t = Hubbard_Hopping_Hamiltonian(t, kwargs.P, kwargs.Q, 'convention', kwargs.convention);
        H_onesite = Hubbard_Onesite_Hamiltonian(pspace, trivspace, U);
        if kwargs.t2 ~= 0
            Hopping_t2 = Hubbard_Hopping_Hamiltonian(kwargs.t2, kwargs.P, kwargs.Q, 'convention', kwargs.convention);
        end
    elseif strcmp(kwargs.symmetries, 'U1_SU2')
        [pspace, ~, trivspace] = get_spaces_Hubbard_SU2(kwargs.P, kwargs.Q, 'D', kwargs.D);
        Hopping_t = Hubbard_Hopping_Hamiltonian_SU2(t, kwargs.P, kwargs.Q, 'convention', kwargs.convention);
        H_onesite = Hubbard_Onesite_Hamiltonian_SU2(pspace, trivspace, U);
        if kwargs.t2 ~= 0
            Hopping_t2 = Hubbard_Hopping_Hamiltonian_SU2(kwargs.t2, kwargs.P, kwargs.Q, 'convention', kwargs.convention);
        end
    elseif strcmp(kwargs.symmetries, 'None_SU2')
        warning('not checked!!')
        [pspace, ~, trivspace] = get_spaces_Hubbard_None_SU2('D1', kwargs.D, 'D2', kwargs.D);
        Hopping_t = Hubbard_Hopping_Hamiltonian_None_SU2(t, 'convention', kwargs.convention, 'bitstring', kwargs.bitstring);
        H_onesite = Hubbard_Onesite_Hamiltonian_None_SU2(pspace, trivspace, U);
        if kwargs.t2 ~= 0
            Hopping_t2 = Hubbard_Hopping_Hamiltonian_None_SU2(kwargs.t2, 'convention', kwargs.convention);
        end
    elseif strcmp(kwargs.symmetries, 'None_U1')
        error('Probably wrong.');
        [pspace, ~, trivspace] = get_spaces_Hubbard_without_U1('D1', kwargs.D, 'D2', kwargs.D);
        Hopping_t = Hubbard_Hopping_Hamiltonian_without_U1(t, 'convention', kwargs.convention);
        H_onesite = Hubbard_Onesite_Hamiltonian_without_U1(pspace, trivspace, U, kwargs.mu);
        if kwargs.t2 ~= 0
            Hopping_t2 = Hubbard_Hopping_Hamiltonian_without_U1(kwargs.t2, kwargs.P, kwargs.Q);
        end
    else
        error('Invalid symmetry')
    end

    if kwargs.V ~= 0
        neighbour_interaction = Hubbard_V_Hamiltonian(pspace, trivspace, kwargs.V, 'convention', kwargs.convention);
        Hopping_t = Hopping_t + neighbour_interaction;
    end
    
    if strcmp(kwargs.system{1}, '1D')
        if kwargs.t2 == 0
            mpo = get_mpo_1D(Hopping_t, H_onesite, 'len', kwargs.len, 'convention', kwargs.convention);
        else
            warning('convention not fixed')
            mpo = get_mpo_1D_ntn(Hopping_t, Hopping_t2, H_onesite, 'len', kwargs.len);
        end
    elseif strcmp(kwargs.system{1}, 'Helix')
        if kwargs.t2 == 0
            mpo = get_mpo_helix(Hopping_t, H_onesite, kwargs.system{2}, 'len', kwargs.len, 'convention', kwargs.convention);
        else
            warning('convention not fixed');
            mpo = get_mpo_helix_ntn(Hopping_t, Hopping_t2, H_onesite, kwargs.system{2});
        end
    elseif strcmp(kwargs.system{1}, 'Cylinder')
        warning('convention not fixed');
        mpo = get_mpo_cylinder(Hopping_t, H_onesite, kwargs.system{2}, 1);
    elseif strcmp(kwargs.system{1}, 'Helix_multiple_rungs')
        assert(kwargs.t2 == 0, 't2 != 0 to be added');
        mpo = get_mpo_helix(Hopping_t, H_onesite, kwargs.system{2}, kwargs.system{3}, 'convention', kwargs.convention);
    elseif strcmp(kwargs.system{1}, 'Cylinder_multiple_rungs')
        assert(kwargs.t2 == 0, 't2 != 0 to be added');
        mpo = get_mpo_cylinder(Hopping_t, H_onesite, kwargs.system{2}, kwargs.system{3}, 'convention', kwargs.convention);
    elseif strcmp(kwargs.system{1}, 'DoubleCylinder')
        warning('convention not fixed');
        mpo = get_mpo_cylinder(Hopping_t, H_onesite, kwargs.system{2}, 2);
    else
        error('System %s not implemented', kwargs.system{1});
    end

    H = InfJMpo(mpo);
end