function cdc_e = Hubbard_operators_SU2(type, P, Q)
    [pspace, ~, trivspace] = get_spaces_Hubbard_SU2(P, Q);
    
    e_charge = ProductCharge(U1(Q), SU2(2), fZ2(1));
    e_space = GradedSpace.new(e_charge, 1, false);

    c_dagger_e = Tensor(pspace, [pspace e_space]);
    c_dagger_e = fill_tensor(c_dagger_e, num2cell([1 1]));

    c_e = conj(c_dagger_e);
    c_e = tpermute(c_e, [3 2 1]);
    cdc_e = contract(c_dagger_e, [-1 1 -4], c_e, [-2 1 -3]);
end