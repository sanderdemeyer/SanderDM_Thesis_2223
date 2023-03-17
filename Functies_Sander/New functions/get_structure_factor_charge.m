function s = get_structure_factor_charge(gs_mps, P, Q, max_dist, symmetries, operator_type, kwargs)
    arguments
        gs_mps
        P
        Q
        max_dist
        symmetries = 'U1_SU2'
        operator_type = 'separate'
        kwargs.corr = []
        kwargs.convention = 'conventional'
        kwargs.system = {'1D'}
        kwargs.qx_size = 500
        kwargs.tol_check = false
    end
    if strcmp(operator_type, 'joint')
        warning('This might be more efficient if the 2 number operators are given separately');
    end

    [O_alpha, O_beta] = get_number_number_operators(P, Q, 'convention', kwargs.convention, 'symmetries', symmetries, 'operator_type', operator_type);


    corr = correlation_function({O_alpha, O_beta}, gs_mps, max_dist, operator_type, kwargs.convention, kwargs.tol_check);
    save('corr_list', 'corr');
    disp('corr_list saved!')

    struct = namedargs2cell(kwargs);
    s = get_structure_factor({O_alpha, O_beta}, gs_mps, max_dist, operator_type, struct{:});
end