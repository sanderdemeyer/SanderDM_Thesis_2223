function nn = Hubbard_V_Hamiltonian_None_SU2(pspace, trivspace, V, kwargs)
    arguments
        pspace
        trivspace
        V
        kwargs.convention = 'conventional'
    end
    number_data = {reshape([0 0; 0 2], 2, 1, 2), 1};

    number_L = Tensor(pspace, [pspace trivspace]);
    number_R = Tensor([pspace trivspace], pspace);

    number_L = fill_tensor(number_L, number_data);
    number_R = fill_tensor(number_R, number_data);

    if strcmp(kwargs.convention, 'first')
        nn = V*contract(number_L, [-3 1 -1], number_R, [-4 1 -2]);
    elseif strcmp(kwargs.convention, 'conventional')
        nn = V*contract(number_L, [-1 1 -4], number_R, [-2 1 -3]);
        %nn2 = contract(number_L, [-1 1 -4], number_R2, [1 -2 -3]);
    else
        error('Convention must be either first or conventional.')
    end
end