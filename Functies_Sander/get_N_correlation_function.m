function corr_list = get_N_correlation_function(gs_mps, pspace, trivspace, site_number, varargin)
    AC1 = gs_mps.AC(1);
    AC2 = gs_mps.AC(2);
    AR1 = gs_mps.AR(1);
    AR2 = gs_mps.AR(2);

    corr_list = correlation_function_new(O, AC1, AC2, AR1, AR2, max_dist);
end