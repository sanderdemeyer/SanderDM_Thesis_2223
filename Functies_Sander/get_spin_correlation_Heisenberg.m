function corr_list = get_spin_correlation_Heisenberg(gs_mps, pspace, max_dist)
    % calculates <up,i up,j> - <up,i><up,j>
    tens_double = Tensor([pspace' pspace'], [pspace' pspace']);
    tens_single = Tensor(pspace', pspace');
    vars_double = [0 0 0 0 0 1];
    vars_single = [0 1];
    H_double = fill_tensor(tens_double, num2cell(vars_double));
    H_single = fill_tensor(tens_single, num2cell(vars_single));
    corr_list = correlation_function(H_double, gs_mps, 2, false, max_dist);
    AC1 = gs_mps.AC(1);
    AC2 = gs_mps.AC(2);
    single = contract(AC1, [1 2 3], conj(AC1), [1 4 3], H_single, [2 4]);
    corr_list = corr_list - single^2;
    return
end