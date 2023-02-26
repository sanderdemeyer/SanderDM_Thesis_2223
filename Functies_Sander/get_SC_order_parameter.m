function Pd = get_SC_order_parameter(mps, N, max_dist, kwargs)
    arguments
        mps
        N
        max_dist
        kwargs.P = 1
        kwargs.Q = 1
    end
    AL = mps.AL;
    AC = mps.AC;
    AR = mps.AR;

    warning('might not be using the correct 4-points correlation functions');
    warning('ri=0 does not have to be multiplied by 2, change this!');
    
    [cdc_up, cdc_down] = Hubbard_operators('cdc_both', kwargs.P, kwargs.Q);
    summation = 0;
    for j = [1, -1, N, -N]
        for k = [1, -1, N, -N]
            min_i = min([N, N + j, N + k]);
            indices = [N-min_i 4; N-min_i 3; N + j-min_i 2; N + k-min_i 1;];
            corr = get_4_site_correlation_function_new(cdc_up, cdc_down, AL, AC, AR, N, indices, max_dist);
            factor_j = 1*(abs(j) == 1) - 2*(abs(j) == N);
            factor_k = 1*(abs(k) == 1) - 2*(abs(k) == N);
            summation = summation - factor_j*factor_k*(corr(1) + 2*sum(corr(2:end))); % Change, this is wrong
        end
    end
    Pd = summation;
end