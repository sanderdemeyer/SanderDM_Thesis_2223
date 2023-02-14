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
        kwargs.U1 = true
        kwargs.mu = 0
    end

   % if kwargs.P == 1000 % change!
   %     [pspace, ~, trivspace] = get_spaces_Hubbard_symmetric(kwargs.P, kwargs.Q, 'D1', kwargs.D, 'D2', kwargs.D);
   % else
    if kwargs.U1
        [pspace, ~, trivspace] = get_spaces_Hubbard_asymmetric(kwargs.P, kwargs.Q, 'D', kwargs.D);
        Hopping_t = Hubbard_Hopping_Hamiltonian(t, kwargs.P, kwargs.Q);
        H_onesite = Hubbard_Onesite_Hamiltonian(pspace, trivspace, U);
    else
        [pspace, ~, trivspace] = get_spaces_Hubbard_without_U1('D1', kwargs.D, 'D2', kwargs.D);
        Hopping_t = Hubbard_Hopping_Hamiltonian_without_U1(t);
        H_onesite = Hubbard_Onesite_Hamiltonian_without_U1(pspace, trivspace, U, kwargs.mu);
    end

    if kwargs.V ~= 0
        disp('V not 0');
        neighbour_interaction = kwargs.V*Hubbard_operators('nn', kwargs.P, kwargs.Q);
        Hopping_t = Hopping_t + neighbour_interaction;
    end
    if kwargs.t2 ~= 0
        disp('t2 not 0');
        Hopping_t2 = Hubbard_Hopping_Hamiltonian(kwargs.t2, kwargs.P, kwargs.Q);
    end
    
    if strcmp(kwargs.system{1}, '1D')
        if kwargs.t2 == 0
            mpo = get_mpo_1D(Hopping_t, H_onesite, 'len', kwargs.len);
        else
            mpo = get_mpo_1D_ntn(Hopping_t, Hopping_t2, H_onesite, 'len', kwargs.len);
        end
    elseif strcmp(kwargs.system{1}, 'Helix')
        if kwargs.t2 == 0
            mpo = get_mpo_helix(Hopping_t, H_onesite, kwargs.system{2});
        else
            mpo = get_mpo_helix_ntn(Hopping_t, Hopping_t2, H_onesite, kwargs.system{2});
        end
    elseif strcmp(kwargs.system{1}, 'Cylinder')
        mpo = get_mpo_cylinder(Hopping_t, H_onesite, kwargs.system{2}, 1);
    elseif strcmp(kwargs.system{1}, 'DoubleCylinder')
        mpo = get_mpo_cylinder(Hopping_t, H_onesite, kwargs.system{2}, 2);
    else
        error('System %s not implemented', kwargs.system{1});
    end

    H = InfJMpo(mpo);
end