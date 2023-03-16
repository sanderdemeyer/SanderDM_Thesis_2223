function [O_alpha, O_beta] = get_number_number_operators(P, Q, kwargs)
    arguments
        P
        Q
        kwargs.convention = 'conventional'
        kwargs.symmetries = 'U1_SU2'
    end
    assert(strcmp(kwargs.convention, 'conventional'), 'Only implemented for convention = conventional. Other conventions are strongly discouraged.');

    if strcmp(kwargs.symmetries, 'U1_SU2')
        [pspace, ~, trivspace] = get_spaces_Hubbard_SU2(P, Q);
    elseif strcmp(kwargs.symmetries, 'U1_U1')
        [pspace, ~, trivspace] = get_spaces_Hubbard_asymmetric(P, Q);
    end

    number_data = num2cell([0 1 2]);
    number_L = Tensor(pspace, [pspace trivspace]);
    number_R = Tensor([pspace trivspace], pspace);

    O_alpha = fill_tensor(number_L, number_data);
    O_beta = fill_tensor(number_R, number_data);
end