function [corr_list_density, SC_lists_average] = order_parameters(gs_mps, kwargs)
    arguments
        gs_mps
        kwargs.symmetries = 'None_SU2'
        kwargs.max_dist = 250
        kwargs.P = 0
        kwargs.Q = 0
    end
    assert(strcmp(kwargs.symmetries, 'None_SU2'), 'Only None_SU2 implemented.');

    O_hole = hole_density_operator();
    corr_list_density = correlation_function({O_hole O_hole}, gs_mps, kwargs.max_dist, 'separate');

    [~, SC_lists_average, ~] = get_SC_order_parameter_april(gs_mps, N, kwargs.P, kwargs.Q, kwargs.max_dist, 0, 0, 'symmetries', kwargs.symmetries);

end