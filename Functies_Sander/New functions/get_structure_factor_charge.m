function [s, corr] = get_structure_factor_charge(gs_mps, P, Q, max_dist, kwargs)
    arguments
        gs_mps
        P
        Q
        max_dist
        kwargs.symmetries = 'U1_SU2'
        kwargs.operator_type = 'separate'
        kwargs.corr = []
        kwargs.convention = 'conventional'
        kwargs.system = {'1D'}
        kwargs.qx_size = 500
        kwargs.tol_check = false
    end
    if strcmp(kwargs.operator_type, 'joint')
        warning('This might be more efficient if the 2 number operators are given separately');
    end

    [O_alpha, O_beta] = get_number_number_operators(P, Q, 'convention', kwargs.convention, 'symmetries', kwargs.symmetries, 'operator_type', kwargs.operator_type);


    corr = correlation_function({O_alpha, O_beta}, gs_mps, max_dist, kwargs.operator_type, 'convention', kwargs.convention, 'tol', kwargs.tol_check);
    disp('corr_list saved!')

    struct = namedargs2cell(kwargs);
    s = get_structure_factor({O_alpha, O_beta}, gs_mps, max_dist, kwargs.operator_type, struct{5:end}, 'corr', corr);
end