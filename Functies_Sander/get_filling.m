function number = get_filling(gs_mps, kwargs)
    arguments
        gs_mps
        kwargs.symmetries = 'None_U1'
        kwargs.D = 1
    end
    AC = gs_mps.AC;
    % SU2 is a boolean and indicates whether spin up and spin down are the
    % same
    w = period(gs_mps);

    if strcmp(kwargs.symmetries, 'None_U1')
        [pspace, ~, trivspace] = get_spaces_Hubbard_None_SU2('D1', kwargs.D, 'D2', kwargs.D);
        O = Tensor(pspace, pspace);
        O_data = {[0 0; 0 2], 1};
        O = fill_tensor(O, O_data);
    else
        error('TBA');
    end
    number_list = cell(1, w);
    for b = 1:w
        number_list{b} = contract(AC(b), [1 2 3], O, [4 2], twist(conj(AC(b)),3), [1 4 3]);
    end
    number = mean(cell2mat(number_list));
end