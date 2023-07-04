function [exp_SC_list, exp_SC_disconnected_list, exp_density_list] = analysis_SC_vs_density_postprocessing(density_list, phi_yy)
    
    l = length(phi_yy);
    assert(l == length(density_list), 'SC and density correlation functions should have the same length')
    
    exp_SC_list = cell(1,l);
    exp_SC_disconnected_list = cell(1,l);
    exp_density_list = cell(1,l);

    for k = 1:l
        [exp_SC, maximum_x, maximum_corr] = fit_correlation_function(phi_yy{k}, 'connected', true);
        [exp_SC_disconnected, maximum_x_disc, maximum_corr_disc] = fit_correlation_function(phi_yy{k}, 'connected', false);
        [exp_density, maximum_x_dens, maximum_corr_dens] = fit_correlation_function(density_list{k}, 'connected', false);%, [5.5 6]);

        exp_SC_list = exp_SC;
        exp_SC_disconnected_list = exp_SC_disconnected;
        exp_density_list = exp_density;
    end
end