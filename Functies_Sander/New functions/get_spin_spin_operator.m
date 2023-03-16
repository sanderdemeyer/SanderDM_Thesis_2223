function O = get_spin_spin_operator(P, Q, kwargs)
    arguments
        P
        Q
        kwargs.convention = 'conventional'
        kwargs.symmetries = 'U1_SU2'
    end
    assert(strcmp(kwargs.convention, 'conventional'), 'Only implemented for convention = conventional. Other conventions are strongly discouraged.');

    if strcmp(kwargs.symmetries, 'U1_SU2')
        pspace = get_spaces_Hubbard_SU2(P, Q);
        O_vars = num2cell([0 0 0 0 0 0 0 0 0 3/4 0 0 0 0 1/4 0 0 0 0 0]);
    elseif strcmp(kwargs.symmetries, 'U1_U1')
        pspace = get_spaces_Hubbard_asymmetric(P, Q);
        O_vars = num2cell([0 0 0 0 0 0 0 0 0 1/4 0 0 0 0 0 1/2 -1/4 0 0 -1/4 1/2 0 0 0 0 0 1/4 0 0 0 0 0 0 0 0 0]);
    end
    O = Tensor([pspace pspace], [pspace pspace]);
    O = fill_tensor(O, O_vars);
end