function [exp_SC_list, exp_density_list, phase_list] = analysis_SC_vs_density(mu)
    folder_base = 'Data structures\Superconductivity - None_SU2\oneband_Hg\Getting the Dome\';
    
    max_dist = 300;
    
    O_hole = hole_density_operator();
    load('operators_SC.mat');

    folder = strcat(folder_base, string(mu), '/');
    number_of_files = length(dir(fullfile(folder,'**','*.mat')));
    files = dir((folder));
    l = length(files);
    k = 1;

    check_trunc_list = cell(1,l-3);
    phi_yy_list_list = cell(1, l-3);
    corr_list_density_list = cell(1, l-3);
    exp_SC_list = cell(1, l-3);
    exp_SC_disconnected_list = cell(1, l-3);
    exp_density_list = cell(1, l-3);
    phase_list = cell(1,l-3);
    phase_disc_list = cell(1,l-3);

    for i = 3:l
        file = files(i);
        naam = file.name;
        name_l = length(naam);
        if name_l == 1
            naam = naam{1};
            name_l = length(naam);
        end
    
        string_mu = strfind(naam,'mu');
        string_trunctotdim = strfind(naam,'trunctotdim');
        string_final = strfind(naam,'final');
    
        mu = str2double(naam(string_mu+3:string_trunctotdim-2));
        trunc = str2double(naam(string_trunctotdim+12:string_final-2));
    
        fprintf('Started with mu = %d, trunc = %0.4f \n', mu, trunc);
        if ~strcmp(naam(name_l-9:name_l), '_final.mat')
            load(strcat(folder, naam), 'mps', 'lambda', 'eta');
            gs_mps = canonicalize(mps, 'Order', 'rl');
            gs_energy = lambda;
        else
            load(strcat(folder, naam), 'gs_mps', 'gs_energy', 'eta');
        end
    
        if eta < 10^(-trunc)
            check_trunc_list{k} = 1;

            phi_yy_list = get_SC_phi_yy(gs_mps, 2, 0, 0, max_dist, 'symmetries', 'None_SU2', 'operators', {L1 L2 L3 L4});
            corr_list_density = correlation_function({O_hole O_hole}, gs_mps, 2*max_dist, 'separate');
            scatter(log(1:2*max_dist), log(corr_list_density-corr_list_density(end))); 
            hold on; 
            scatter(log(1:max_dist), log(phi_yy_list)); 
            hold off; 
            legend('density', 'pair-pair');
            xlabel('log of the distance');
            ylabel('log of the correlation function')
            title(sprintf('Comparison between the pair-pair and hole density correlation functions - mu = %d, trunc = %d', mu, round(trunc,2)));
            
            exp_SC = fit_correlation_function(phi_yy_list, 'connected', true);
            exp_SC_disconnected = fit_correlation_function(phi_yy_list, 'connected', false);
            exp_density = fit_correlation_function(corr_list_density, 'connected', false);

            phi_yy_list_list{k} = phi_yy_list;
            corr_list_density_list{k} = corr_list_density;
            exp_SC_list{k} = exp_SC;
            exp_SC_disconnected_list{k} = exp_SC_disconnected;
            exp_density_list{k} = exp_density;

            if exp_SC < exp_density
                phase_list{k} = 'Non_SC';
            else
                phase_list{k} = 'SC';
            end

            if exp_SC_disconnected < exp_density
                phase_disc_list{k} = 'Non_SC';
            else
                phase_disc_list{k} = 'SC';
            end
        else
            check_trunc_list{k} = 0;
        end
    end
    save(sprintf('analysis_mu_%d_trunc_%.4f.mat', mu, trunc))
end