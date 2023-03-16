function s = get_structure_factor_spin(gs_mps, P, Q, max_dist, operator_type, symmetries, kwargs)
    arguments
        gs_mps
        P
        Q
        max_dist
        operator_type
        symmetries = 'U1_SU2'
        kwargs.corr = []
        kwargs.convention = 'conventional'
        kwargs.system = {'1D'}
        kwargs.qx_size = 500
        kwargs.tol_check = false
    end

    O = get_spin_spin_operator(P, Q, 'convention', kwargs.convention, 'symmetries', symmetries);

    struct = namedargs2cell(kwargs);
    s = get_structure_factor(O, gs_mps, max_dist, operator_type, struct{:});
end