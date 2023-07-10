function [s, corr] = get_structure_factor_spin(gs_mps, P, Q, max_dist, kwargs)
    arguments
        gs_mps
        P
        Q
        max_dist
        kwargs.symmetries = 'U1_SU2'
        kwargs.corr = []
        kwargs.convention = 'conventional'
        kwargs.system = {'1D'}
        kwargs.qx_size = 500
        kwargs.tol_check = false
    end

    O = get_spin_spin_operator(P, Q, 'convention', kwargs.convention, 'symmetries', kwargs.symmetries);

    corr = correlation_function(O, gs_mps, max_dist, 'twosite', 'convention', kwargs.convention, 'tol', kwargs.tol_check);
    disp('corr_list saved!')

    struct = namedargs2cell(kwargs);
    s = get_structure_factor(O, gs_mps, max_dist, 'twosite', struct{3:end}, 'corr', corr);
end