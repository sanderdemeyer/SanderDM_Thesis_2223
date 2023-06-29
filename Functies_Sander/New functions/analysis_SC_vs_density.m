function [exp_SC_list, exp_density_list, phase_list] = analysis_SC_vs_density(mu)
    folder_base = 'Data structures\Superconductivity - None_SU2\oneband_Hg\Getting the Dome\';
    folder_base = '\home\sanddmey\Data_and_analysis_getting_the_dome_2023_06_09/';
    folder_base = 'Data structures\Superconductivity - None_SU2\oneband_Hg\Getting the Dome\';
    
    max_dist = 500;
    
    O_hole = hole_density_operator();
    load('operators_SC.mat');

    folder = strcat(folder_base, string(mu), '/');
    %folder = 'Data structures\Superconductivity - None_SU2\mu_m1.6428/';
    disp(folder);
    %folder = folder_base;
    number_of_files = length(dir(fullfile(folder,'**','*.mat')));
    %files = dir(num2str(mu));
    files = dir(num2str(folder));
    l = length(files);
    fprintf('The number of files in the dirgectory for mu = %d is %d \n', mu, l);
    k = 1;

    check_trunc_list = cell(1,l-3);
    phi_yy_list_list = cell(1, l-3);
    corr_list_density_list = cell(1, l-3);
    exp_SC_list = cell(1, l-3);
    exp_SC_disconnected_list = cell(1, l-3);
    exp_density_list = cell(1, l-3);
    phase_list = cell(1,l-3);
    phase_disc_list = cell(1,l-3);

    tot_bonddim_list = cell(1,l-3);
    filling_list = cell(1,l-3);

    maximum_x_list = cell(1,l-3);
    maximum_corr_list = cell(1,l-3);
    maximum_x_disc_list = cell(1,l-3);
    maximum_corr_disc_list = cell(1,l-3);
    maximum_x_dens_list = cell(1,l-3);
    maximum_corr_dens_list = cell(1,l-3);

    all_sectors{k} = cell(1,l-3);
    deltas{k} = cell(1,l-3);
    deltas_42{k} = cell(1,l-3);
    epsilon1s{k} = cell(1,l-3);

    for i = 3:l
        file = files(i);
        naam = file.name;
        name_l = length(naam);
        if name_l == 1
            naam = naam{1};
            name_l = length(naam);
        end
        disp(naam);
        string_mu = strfind(naam,'mu');
        string_trunctotdim = strfind(naam,'trunctotdim');
        string_final = strfind(naam,'final');
    
        %mu = str2double(naam(string_mu+3:string_trunctotdim-2));
        if isempty(string_final)
            string_mat = strfind(naam, 'mat');
            trunc = str2double(naam(string_trunctotdim+12:string_mat-2));
        else
            trunc = str2double(naam(string_trunctotdim+12:string_final-2));
        end

        fprintf('Started with mu = %d, trunc = %0.4f \n', mu, trunc);
        if ~strcmp(naam(name_l-9:name_l), '_final.mat')
            %load(strcat(folder, naam), 'mps', 'lambda', 'eta');
            load(naam, 'mps', 'lambda', 'eta');
            gs_mps = canonicalize(mps, 'Order', 'rl');
            gs_energy = lambda;
        else
            %load(strcat(folder, naam), 'gs_mps', 'gs_energy', 'eta');
            load(naam, 'gs_mps', 'gs_energy', 'eta');
        end
    
        if eta < 10^(-trunc)% && trunc > 4.25 
            check_trunc_list{k} = 1;

            phi_yy_list = get_SC_phi_yy(gs_mps, 2, 0, 0, max_dist, 'symmetries', 'None_SU2', 'operators', {L1 L2 L3 L4});
            corr_list_density = correlation_function({O_hole O_hole}, gs_mps, max_dist, 'separate');
            scatter((1:max_dist), log(corr_list_density-corr_list_density(end))); 
            hold on; 
            scatter((1:max_dist), log(phi_yy_list)); 
            hold off; 
            legend('density', 'pair-pair');
            xlabel('log of the distance');
            ylabel('log of the correlation function')
            title(sprintf('Comparison between the pair-pair and hole density correlation functions - mu = %d, trunc = %d', mu, round(trunc,2)));
            
            %[exp_SC, maximum_x, maximum_corr] = fit_correlation_function_contour(phi_yy_list, 'connected', true);
            %[exp_SC_disconnected, maximum_x_disc, maximum_corr_disc] = fit_correlation_function_contour(phi_yy_list, 'connected', false);
            %[exp_density, maximum_x_dens, maximum_corr_dens] = fit_correlation_function_contour(corr_list_density, 'connected', false);
            [exp_SC, maximum_x, maximum_corr] = fit_correlation_function(phi_yy_list, 'connected', true);
            [exp_SC_disconnected, maximum_x_disc, maximum_corr_disc] = fit_correlation_function(phi_yy_list, 'connected', false);
            [exp_density, maximum_x_dens, maximum_corr_dens] = fit_correlation_function(corr_list_density, 'connected', false, 'edges', [5.5 6]);

            phi_yy_list_list{k} = phi_yy_list;
            corr_list_density_list{k} = corr_list_density;
            exp_SC_list{k} = exp_SC;
            exp_SC_disconnected_list{k} = exp_SC_disconnected;
            exp_density_list{k} = exp_density;

            tot_bonddim_list{k} = get_tot_bonddim(gs_mps);
            filling_list{k} = get_filling(gs_mps);

            disp('Started with transfereigs');

            [V, D] = transfereigs(gs_mps, gs_mps, 5);
            epsilons = zeros(1,5);
            for j = 1:5
                epsilons(j) = -log(norm(D(j,j)));
            end
            all_sectors{k} = D;
            deltas{k} = epsilons(3) - epsilons(2);
            deltas_42{k} = epsilons(4) - epsilons(2);
            epsilon1s{k} = epsilons(2);

            maximum_x_list{k} = maximum_x;
            maximum_corr_list{k} = maximum_corr;
            maximum_x_disc_list{k} = maximum_x_disc;
            maximum_corr_disc_list{k} = maximum_corr_disc;
            maximum_x_dens_list{k} = maximum_x_dens;
            maximum_corr_dens_list{k} = maximum_corr_dens;

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
    
            save(sprintf('analysis_mu_%.4f_trunc_%.4f.mat', mu, trunc))
        else
            check_trunc_list{k} = 0;
        end
        k = k+1;
    end
    save(sprintf('analysis_mu_%.4f.mat', mu))
end