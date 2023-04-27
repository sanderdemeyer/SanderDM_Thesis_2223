function cdc = Hubbard_Hopping_Hamiltonian_None_SU2(t, kwargs)
    arguments
        t
        kwargs.bitstring = 0
        kwargs.convention = 'conventional'
    end
        
    % based on comparison with H_double([1, 3, 4, 2], [1, 3, 4, 2])

    assert(strcmp(kwargs.convention, 'conventional'), 'Only convention = conventional is implemented.');

    pspace = get_spaces_Hubbard_None_SU2();

    singlet_charge = ProductCharge(SU2(2), fZ2(1));
    singlet_space = GradedSpace.new(singlet_charge, 1, false);

    c_dagger = Tensor(pspace, [pspace singlet_space]);
    c = Tensor([singlet_space pspace], pspace);

    c_dagger_data = {reshape([0 sqrt(2)], 2, 1), reshape([1 0], 1, 1, 2)};
    c_data = {reshape([0 -sqrt(2)], 1, 1, 2), [1 0]};

    c_dagger = fill_tensor(c_dagger, c_dagger_data);
    c = fill_tensor(c, c_data);

    c_dagger_c = contract(c_dagger, [-1 1 -4], c, [1 -2 -3]);
    cdc = -t*(c_dagger_c + tpermute(conj(c_dagger_c), [4 3 2 1]));
    return
end