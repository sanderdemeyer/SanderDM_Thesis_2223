function [O, O2, O3, O4] = Hubbard_operators(type, P, Q)
    warning('convention = first');
    
    O2 = 0;
    O3 = 0;
    O4 = 0;
    %if P == 1
    %    [pspace] = get_spaces_Hubbard_symmetric(P, Q);
    %else
    [pspace, ~, trivspace] = get_spaces_Hubbard_asymmetric(P, Q);

    up_charge = ProductCharge(U1(Q), U1(1), fZ2(1));
    down_charge = ProductCharge(U1(Q), U1(-1), fZ2(1));
    up_flip_charge = ProductCharge(U1(Q-P), U1(2), fZ2(0));
    down_flip_charge = ProductCharge(U1(Q-P), U1(-2), fZ2(0));
    updown_charge = ProductCharge(U1(4*Q-2*P), U1(0), fZ2(0));

    up_space = GradedSpace.new(up_charge, 1, false);
    down_space = GradedSpace.new(down_charge, 1, false);
    up_flip_space = GradedSpace.new(up_flip_charge, 1, false);
    down_flip_space = GradedSpace.new(down_flip_charge, 1, false);
    updown_space = GradedSpace.new(updown_charge, 1, false);

    c_dagger_up = Tensor(pspace, [pspace up_space]);
    c_dagger_down = Tensor(pspace, [pspace down_space]);
    up_flip = Tensor(pspace, [pspace up_flip_space]);
    down_flip = Tensor(pspace, [pspace down_flip_space]);
    c_dagger_updown = Tensor(pspace, [pspace updown_space]);

    c_dagger_up = fill_tensor(c_dagger_up, num2cell([1 1]));
    c_dagger_down = fill_tensor(c_dagger_down, num2cell([1 -1]));
    c_dagger_updown = fill_tensor(c_dagger_updown, num2cell(1));
    up_flip = fill_tensor(up_flip, num2cell(1));
    down_flip = fill_tensor(down_flip, num2cell(1));
    up_flip_conj = conj(up_flip);
    up_flip_conj = tpermute(up_flip_conj, [3 2 1]);
    down_flip_conj = conj(down_flip);
    down_flip_conj = tpermute(down_flip_conj, [3 2 1]);

    c_up = conj(c_dagger_up);
    c_up = tpermute(c_up, [3 2 1]);
    c_down = conj(c_dagger_down);
    c_down = tpermute(c_down, [3 2 1]);
    c_updown = conj(c_dagger_updown);
    c_updown = tpermute(c_updown, [3 2 1]);

    cdc_up = contract(c_dagger_up, [-3 1 -1], c_up, [-4 1 -2]);
    cdc_down = contract(c_dagger_down, [-3 1 -1], c_down, [-4 1 -2]);

    cdc_base = cdc_up + cdc_down;
    cdc = cdc_base + tpermute(conj(cdc_base), [3 4 1 2]);

    cdc_L = Tensor(pspace', [pspace' trivspace]);
    cdc_R = Tensor([pspace' trivspace], pspace');

    n_L = fill_tensor(cdc_L, num2cell([0 1 1 -2]));
    n_R = fill_tensor(cdc_R, num2cell([0 1 1 -2]));
    nn2 = contract(n_L, [-1 1 -3], n_R, [-2 1 -4]);

    number_L = fill_tensor(cdc_L, num2cell([0 1 1 2]));
    number_R = fill_tensor(cdc_R, num2cell([0 1 1 2]));
    nn = contract(number_L, [-1 1 -3], number_R, [-2 1 -4]);

    ttest = Tensor([pspace', pspace'], [pspace', pspace']);
    var = num2cell([0 1 0 0 1 1 0 0 1 0 0 1 -1 0 -1 0 0 -1 1 0 0 1 0 1 -1 0 0 -1 0 0 -1 -1 0 0 -1 0]);
    H = fill_tensor(ttest, var);



    if strcmp('c_dagger_up', type)
        O = c_dagger_up;
    elseif strcmp('c_dagger_down', type)
        O = c_dagger_down;
    elseif strcmp('c_up', type)
        O = c_up;
    elseif strcmp('c_down', type)
        O = c_down;
    elseif strcmp('all basic', type)
        O = c_dagger_up;
        O2 = c_dagger_down;
        O3 = c_up;
        O4 = c_down;
    elseif strcmp('cdc_both', type)
        O = cdc_up;
        O2 = cdc_down;
    elseif strcmp('up_flip', type)
        O = up_flip;
    elseif strcmp('down_flip', type)
        O = down_flip;
    elseif strcmp('up_down_flip', type)
        O = contract(up_flip, [-3 1 -1], up_flip_conj, [-4 1 -2]);
    elseif strcmp('down_up_flip', type)
        O = contract(down_flip, [-3 1 -1], down_flip_conj, [-4 1 -2]);
    elseif strcmp('cdc_up', type)
        O = cdc_up;
    elseif strcmp('cdc_down', type)
        O = cdc_down;
    elseif strcmp('cdc_updown', type)
        O = contract(c_dagger_updown, [-3 1 -1], c_updown, [-4 1 -2]);
    elseif strcmp('cdc', type)
        O = cdc;
    elseif strcmp('nn', type)
        O = nn2;
    else
        error('Type does not exist');
    end

    %{

    space1 = GradedSpace.new(up_charge, 1, false);
    space2 = GradedSpace.new(up_charge, 2, false);
    space3 = GradedSpace.new(up_charge, 3, false);
    space4 = GradedSpace.new(up_charge, 4, false);
    testje = Tensor([space1', space2'], [space4', space3']);
    testje_conj = conj(testje);

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

