function cdc = Hubbard_Hopping_Hamiltonian(t, P, Q, kwargs)
    arguments
        t
        P
        Q
        kwargs.convention = 'conventional'
    end
    % Convention: first was used prior to 16 february 2023. 
    % Conventional is what is should be, with indices 1 2 3 4 meaning
    % bottom left, bottom right, top right, top left, respectively
    % P and Q signify filling f = P/Q

    %[pspace, vspaces, trivspace, prodspace, fusion_trees] = get_spaces('Hubbard', false, P, Q, 12, 3);
    %[pspace, vspaces, trivspace, prodspace, fusion_trees] = get_spaces('Hubbard_asymmetric', P, Q);
    [pspace, vspaces, trivspace, prodspace, fusion_trees] = get_spaces_Hubbard_asymmetric(P, Q);    
    up_charge = ProductCharge(U1(Q), U1(1), fZ2(1));
    down_charge = ProductCharge(U1(Q), U1(-1), fZ2(1));

    up_space = GradedSpace.new(up_charge, 1, false);
    down_space = GradedSpace.new(down_charge, 1, false);

    c_dagger_up = Tensor(pspace, [pspace up_space]);
    c_dagger_down = Tensor(pspace, [pspace down_space]);

    c_dagger_up = fill_tensor(c_dagger_up, num2cell([1 1]));
    c_dagger_down = fill_tensor(c_dagger_down, num2cell([1 -1]));
    
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
    return
    %{
    cdc_up_correct = contract(c_dagger_up, [-1 1 -4], c_up, [-2 1 -3]);
    cdc_down_correct = contract(c_dagger_down, [-1 1 -4], c_down, [-2 1 -3]);
    cdc_base_correct = cdc_up_correct + cdc_down_correct;
    cdc_correct = cdc_base_correct + tpermute(conj(cdc_base_correct), [4 3 2 1]);
    cdc_correct = -t*cdc_correct;
    %}
    %{
    ttest = Tensor([pspace', pspace'], [pspace', pspace']);
    var = num2cell([0 1 0 0 1 1 0 0 1 0 0 1 -1 0 -1 0 0 -1 1 0 0 1 0 1 -1 0 0 -1 0 0 -1 -1 0 0 -1 0]);
    H = fill_tensor(ttest, var);


    space1 = GradedSpace.new(up_charge, 1, false);
    space2 = GradedSpace.new(up_charge, 2, false);
    space3 = GradedSpace.new(up_charge, 3, false);
    space4 = GradedSpace.new(up_charge, 4, false);
    testje = Tensor([space1', space2'], [space4', space3']);
    testje_conj = conj(testje);
    %}

    %{
    %prodspace = pspace * pspace;
    c_dagger = Tensor(pspace, [prodspace pspace]);
    c = Tensor([prodspace pspace], pspace);
    
    c_dagger = tpermute(c_dagger, [2 3 1], [1 2]);
    c = tpermute(c, [3 1 2], [2 1]);
    
    data_dagger_up =  num2cell([0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0]);
    c_dagger_up = fill_tensor(c_dagger, data_dagger_up);
    
    data_dagger_down = num2cell([0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0]);
    c_dagger_down = fill_tensor(c_dagger, data_dagger_down);
    
    data_down = num2cell([0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0]);
    c_down = fill_tensor(c, data_down);
    
    data_up = num2cell([0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0]);
    c_up = fill_tensor(c, data_up);
    
    %%
    cdc_up_ij = contract(c_dagger_up, [-1 1 -3], c_up, [-2 1 -4]);
    cdc_up_ji = contract(c_dagger_up, [-2 1 -4], c_up, [-1 1 -3]);
    cdc_down_ij = contract(c_dagger_down, [-1 1 -3], c_down, [-2 1 -4]);
    cdc_down_ji = contract(c_dagger_down, [-2 1 -4], c_down, [-1 1 -3]);

    %{
    % Probably wrong
    cdc_up_ij = contract(c_dagger_up, [-1 1 -4], c_up, [-2 1 -3]);
    cdc_up_ji = contract(c_dagger_up, [-2 1 -3], c_up, [-1 1 -4]);
    cdc_down_ij = contract(c_dagger_down, [-1 1 -4], c_down, [-2 1 -3]);
    cdc_down_ji = contract(c_dagger_down, [-2 1 -3], c_down, [-1 1 -4]);
    %}
    
    cdc = -t*(cdc_up_ij + cdc_up_ji + cdc_down_ij + cdc_down_ji)/2;

    %%
    %{
    c_test = Tensor([pspace' pspace'], [pspace' pspace']);
    c_up_ij_data = num2cell([0 0 0 0 0 1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0]);
    c_test = fill_tensor(c_test, c_up_ij_data);
    %}
    
    %%
    %{
    ctest2 = Tensor([pspace' pspace'], [pspace' pspace']);
    ctestt2 = fill_tensor(ctest2, num2cell(-(t/2)*[0 1 0 0 1 1 0 0 1 0 0 1 1 0 1 0 0 1 1 0 0 1 0 1 1 0 0 1 0 0 1 1 0 0 1 0]));
    disp(reshape(double(ctestt2), 16, 16));
    cdc = ctestt2;
    return
    %}

    %%
    ttest = Tensor([pspace', pspace'], [pspace', pspace']);
    %var = num2cell(1:36);
    var = num2cell([0 1 0 0 1 1 0 0 1 0 0 1 -1 0 -1 0 0 -1 1 0 0 1 0 1 -1 0 0 -1 0 0 -1 -1 0 0 -1 0]);
    ttest = fill_tensor(ttest, var);
    [tc, tb] = tensorblocks(ttest);

    %}
end

