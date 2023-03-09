function nn = Hubbard_V_Hamiltonian(pspace, trivspace, V, kwargs)
    arguments
        pspace
        trivspace
        V
        kwargs.convention = 'conventional'
    end
    number_data = num2cell([0 1 2]);
    number_L = Tensor(pspace, [pspace trivspace]);
    number_R = Tensor([pspace trivspace], pspace);
    number_R2 = Tensor([trivspace pspace], pspace);

    number_L = fill_tensor(number_L, number_data);
    number_R = fill_tensor(number_R, number_data);
    %number_R2 = fill_tensor(number_R2, number_data);

    if strcmp(kwargs.convention, 'first')
        nn = V*contract(number_L, [-3 1 -1], number_R, [-4 1 -2]);
    elseif strcmp(kwargs.convention, 'conventional')
        nn = V*contract(number_L, [-1 1 -4], number_R, [-2 1 -3]);
        %nn2 = contract(number_L, [-1 1 -4], number_R2, [1 -2 -3]);
    else
        error('Convention must be either first or conventional.')
    end

    %{

    if strcmp(kwargs.symmetries, 'U1_U1')
        if strcmp(kwargs.convention, 'first')
            warning('use with care')
            cdc_L = Tensor(pspace', [pspace' trivspace]);
            cdc_R = Tensor([pspace' trivspace], pspace');
            number_L = fill_tensor(cdc_L, num2cell([0 1 1 2]));
            number_R = fill_tensor(cdc_R, num2cell([0 1 1 2]));
            nn = V*contract(number_L, [-1 1 -3], number_R, [-2 1 -4]);

        elseif strcmp(kwargs.convention, 'conventional')
            cdc_L = Tensor(pspace, [pspace trivspace]);
            cdc_R = Tensor([pspace trivspace], pspace);
            number_L = fill_tensor(cdc_L, num2cell([0 1 1 2]));
            number_R = fill_tensor(cdc_R, num2cell([0 1 1 2]));
            nn = V*contract(number_L, [-1 1 -4], number_R, [-2 1 -3]);
        else
            error('Invalid symmetry');
        end
        %{   
        warning('not checked, use with care.')
        n_L = fill_tensor(cdc_L, num2cell([0 1 1 -2]));
        n_R = fill_tensor(cdc_R, num2cell([0 1 1 -2]));
        nn2 = contract(n_L, [-1 1 -3], n_R, [-2 1 -4]);
        %}
    elseif strcmp(kwargs.symmetries, 'U1_SU2')
        nn = Tensor([pspace' pspace'], [pspace' pspace']);
        %data = num2cell([0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 3 5 0 0 7 0 0 0 0 0 0 9 0 11 13 0 0 15 17 0 19]);
        data = num2cell(V*[0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 -1 0 0 -1 0 0 0 0 0 0 -1 0 2 2 0 0 2 2 0 4]);
        nn = fill_tensor(nn, data);
    elseif strcmp(symmetries, 'None_U1')
        error('TBA');
    else
        error('Invalid symmetry')
    end

    %}
end