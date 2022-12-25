function cdc = Hubbard_operators(t)
    [pspace, vspaces, trivspace, prodspace, fusion_trees] = get_spaces('Hubbard', false, 1, 1, 12, 3);
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
    %{
    cdc_up_ij = contract(c_dagger_up, [-1 1 -3], c_up, [-2 1 -4]);
    cdc_up_ji = contract(c_dagger_up, [-2 1 -4], c_up, [-1 1 -3]);
    cdc_down_ij = contract(c_dagger_down, [-1 1 -3], c_down, [-2 1 -4]);
    cdc_down_ji = contract(c_dagger_down, [-2 1 -4], c_down, [-1 1 -3]);
    %}
    cdc_up_ij = contract(c_dagger_up, [-1 1 -4], c_up, [-2 1 -3]);
    cdc_up_ji = contract(c_dagger_up, [-2 1 -3], c_up, [-1 1 -4]);
    cdc_down_ij = contract(c_dagger_down, [-1 1 -4], c_down, [-2 1 -3]);
    cdc_down_ji = contract(c_dagger_down, [-2 1 -3], c_down, [-1 1 -4]);

    
    cdc = -t*(cdc_up_ij + cdc_up_ji + cdc_down_ij + cdc_down_ji)/2;

    %%
    %{
    ctest2 = Tensor([pspace' pspace'], [pspace' pspace']);
    ctestt2 = fill_tensor(ctest2, num2cell(-(t/2)*[0 1 0 0 1 1 0 0 1 0 0 1 1 0 1 0 0 1 1 0 0 1 0 1 1 0 0 1 0 0 1 1 0 0 1 0]));
    disp(reshape(double(ctestt2), 16, 16));
    cdc = ctestt2;
    return
    %}
end