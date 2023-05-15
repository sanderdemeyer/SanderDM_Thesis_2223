function cdc = Hubbard_Hopping_Hamiltonian_SU2(t, P, Q, kwargs)
    arguments
        t
        P
        Q
        kwargs.convention = 'conventional'
    end
        
    pspace = get_spaces_Hubbard_SU2(P, Q);

    %{
    one_charge = ProductCharge(U1(Q), SU2(3), fZ2(1));
    half_charge = ProductCharge(U1(Q), SU2(2), fZ2(1));
    zero_charge = ProductCharge(U1(Q), SU2(1), fZ2(1));

    one_space = GradedSpace.new(one_charge, 1, false);
    half_space = GradedSpace.new(half_charge, 1, false);
    zero_space = GradedSpace.new(zero_charge, 1, false);

    c_dagger_up = Tensor(pspace, [pspace half_space]);
    c_dagger_down = Tensor(pspace, [pspace half_space]);

    c_dagger_up = fill_tensor(c_dagger_up, num2cell([20 1]));
    c_dagger_down = fill_tensor(c_dagger_down, num2cell([1 -1]));
    
    c_up = conj(c_dagger_up);
    c_up = tpermute(c_up, [3 2 1]);
    c_down = conj(c_dagger_down);
    c_down = tpermute(c_down, [3 2 1]);

    cdc_up = contract(c_dagger_up, [-3 1 -1], c_up, [-4 1 -2]);
    cdc_down = contract(c_dagger_down, [-3 1 -1], c_down, [-4 1 -2]);

    cdc_base = cdc_up + cdc_down;
    cdc = cdc_base + tpermute(conj(cdc_base), [3 4 1 2]);
    cdc = -t*cdc;

    cdc_up_correct = contract(c_dagger_up, [-1 1 -4], c_up, [-2 1 -3]);
    cdc_down_correct = contract(c_dagger_down, [-1 1 -4], c_down, [-2 1 -3]);
    cdc_base_correct = cdc_up_correct + cdc_down_correct;
    cdc_correct = cdc_base_correct + tpermute(conj(cdc_base_correct), [4 3 2 1]);
    cdc_correct = -t*cdc_correct;
    %}


    %cdc_old = Tensor([pspace' pspace'], [pspace' pspace']);
    cdc = Tensor([pspace' pspace' pspace pspace], []);


    %cdc = fill_tensor(cdc, num2cell([2 1 4 3 5 6 7 8 9 10 11 12 13 14 16 15 18 20 17 22]));
    cdc = fill_tensor(cdc, num2cell(t*sqrt(2)*[0 -1 0 0 1 0 0 -1 -1 0 -1 0 1 0 0 1 0 0 1 0]));

    %cdc = fill_tensor(cdc, num2cell(-t*[0 1 0 0 1 0 -sqrt(2) 0 -sqrt(2) 0 -sqrt(2) 0 -sqrt(2) 0 0 -1 0 0 -1 0]));
    %cdc = fill_tensor(cdc, num2cell(-t*[0 1 0 0 1 0 -sqrt(2) 0 -sqrt(2) 0 -sqrt(2) 0 -sqrt(2) 0 0 -1 0 0 -1 0]));
    if strcmp(kwargs.convention, 'conventional')
        cdc = tpermute(cdc, [4 3 1 2]);
    elseif strcmp(kwargs.convention, 'first')
        warning('convention: first');
    else
        error('Convention not defined')
    end
    return

    %{
    tens = Tensor([pspace', pspace'], [pspace', pspace']);
    var = num2cell(t*[0 1 0 0 1 1 0 0 1 0 0 1 -1 0 -1 0 0 -1 1 0 0 1 0 1 -1 0 0 -1 0 0 -1 -1 0 0 -1 0]);
    H = fill_tensor(tens, var);


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

