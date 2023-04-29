function cdc = Hubbard_Hopping_Hamiltonian_SU2_nogeenkeer(t, P, Q, kwargs)
    arguments
        t
        P
        Q
        kwargs.convention = 'conventional'
    end
        
    %[pspace, vspaces, trivspace, prodspace, fusion_trees] = get_spaces('Hubbard', false, P, Q, 12, 3);
    pspace = get_spaces_Hubbard_SU2(P, Q);


    singlet_charge = ProductCharge(U1(P), SU2(2), fZ2(1));
    singlet_space = GradedSpace.new(singlet_charge, 1, false);

    c_dagger = Tensor(pspace, [pspace singlet_space]);
    c = Tensor([singlet_space pspace], pspace);

    %c_dagger_data = {1 sqrt(2)};
    %c_data = {1 -sqrt(2)};

    c_dagger_data = {1 sqrt(2)};
    c_data = {1 sqrt(2)};

    c_dagger = fill_tensor(c_dagger, c_dagger_data);
    c = fill_tensor(c, c_data);

    c_dagger_c = contract(c_dagger, [-1 1 -4], c, [1 -2 -3]);
    cdc = -t*(c_dagger_c + tpermute(conj(c_dagger_c), [4 3 2 1]));

end

