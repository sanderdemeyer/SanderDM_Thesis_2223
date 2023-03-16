function s = get_structure_factor_charge(gs_mps, P, Q, max_dist, symmetries, kwargs)
    arguments
        gs_mps
        P
        Q
        max_dist
        symmetries = 'U1_SU2'
        kwargs.corr = []
        kwargs.convention = 'conventional'
        kwargs.system = {'1D'}
        kwargs.qx_size = 500
        kwargs.tol_check = false
    end
    warning('This might be more efficient if the 2 number operators are given separately');
    
    [O_alpha, O_beta] = get_number_number_operators(P, Q, 'convention', kwargs.convention, 'symmetries', symmetries);

    struct = namedargs2cell(kwargs);
    s = get_structure_factor({O_alpha, O_beta}, gs_mps, max_dist, 'joint', struct{:});
end