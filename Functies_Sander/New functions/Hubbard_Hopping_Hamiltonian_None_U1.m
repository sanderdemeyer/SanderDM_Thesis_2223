function cdc = Hubbard_Hopping_Hamiltonian_None_U1(t, kwargs)
    arguments
        t
        kwargs.convention = 'conventional'
    end

    %[pspace, vspaces, trivspace, prodspace, fusion_trees] = get_spaces('Hubbard', false, P, Q, 12, 3);
    [pspace, ~, trivspace] = get_spaces_Hubbard_without_U1();
    up_charge = ProductCharge(U1(1), fZ2(1));
    down_charge = ProductCharge(U1(-1), fZ2(1));

    up_space = GradedSpace.new(up_charge, 1, false);
    down_space = GradedSpace.new(down_charge, 1, false);

    c_dagger_up = Tensor(pspace, [pspace up_space]);
    c_dagger_down = Tensor(pspace, [pspace down_space]);

    c_dagger_up = fill_tensor(c_dagger_up, {[0;1];reshape([1 0], 1, 1, 2)});
    c_dagger_down = fill_tensor(c_dagger_down, {[1;0];reshape([0 -1], 1, 1, 2)});
    
    c_up = conj(c_dagger_up);
    c_up = tpermute(c_up, [3 2 1]);
    c_down = conj(c_dagger_down);
    c_down = tpermute(c_down, [3 2 1]);

    if strcmp(kwargs.convention, 'first')
        cdc_up = contract(c_dagger_up, [-3 1 -1], c_up, [-4 1 -2]);
        cdc_down = contract(c_dagger_down, [-3 1 -1], c_down, [-4 1 -2]);
    
        cdc_base = cdc_up + cdc_down;
        cdc = cdc_base + tpermute(conj(cdc_base), [3 4 1 2]);
        cdc = -t*cdc;
    elseif strcmp(kwargs.convention, 'conventional')
        cdc_up = contract(c_dagger_up, [-1 1 -4], c_up, [-2 1 -3]);
        cdc_down = contract(c_dagger_down, [-1 1 -4], c_down, [-2 1 -3]);
    
        cdc_base = cdc_up + cdc_down;
        cdc = cdc_base + tpermute(conj(cdc_base), [4 3 2 1]);
        cdc = -t*cdc;
    else
        error('Convention must be either first or conventional.')
    end

end

