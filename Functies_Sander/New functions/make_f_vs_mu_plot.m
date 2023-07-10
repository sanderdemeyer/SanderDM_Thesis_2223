function [mu_list, filling_list] = make_f_vs_mu_plot()    
    folder = 'Data structures\Superconductivity - None_SU2\mu_m1.6428/';
    folder = 'Data structures\Superconductivity - None_SU2\simple Hubbard - conclusion - mu -1.6428/';
    folder = 'Data structures\Superconductivity - None_SU2\oneband_Hg\bigger bond dimension/';
    folder = 'Data structures\Superconductivity - None_SU2\oneband_Hg\fillings greater than 0.875\slow increase/';
    %folder = 'Data structures\Superconductivity - None_SU2\oneband_Hg\other fillings/';
    files = dir(num2str(folder));
    l = length(files);
    k = 1;

    mu_list = cell(1,l-3);
    filling_list = cell(1,l-3);

    for i = 3:l
        file = files(i);
        naam = file.name;
        name_l = length(naam);
        if name_l == 1
            naam = naam{1};
            name_l = length(naam);
        end
        if name_l < 20
            mu_list = cell2mat(mu_list);
            filling_list = cell2mat(filling_list);
            return
        end
        disp(naam);
        string_mu = strfind(naam,'mu');
        string_rungs = strfind(naam,'rungs');
        string_trunctotdim = strfind(naam,'trunctotdim');
        string_cut = strfind(naam,'cut');
        string_final = strfind(naam,'final');
        string_mat = strfind(naam, 'mat');
    

        mu = str2double(naam(string_mu+3:string_trunctotdim-2));

        if ~strcmp(naam(name_l-9:name_l), '_final.mat')
            %load(strcat(folder, naam), 'mps', 'lambda', 'eta');
            load(naam, 'mps', 'lambda', 'eta');
            gs_mps = canonicalize(mps, 'Order', 'rl');
            gs_energy = lambda;
        else
            %load(strcat(folder, naam), 'gs_mps', 'gs_energy', 'eta');
            load(naam, 'gs_mps', 'gs_energy', 'eta');
        end

        mu_list{k} = mu;
        filling_list{k} = get_filling(gs_mps);

        k = k+1;
    end
    mu_list = cell2mat(mu_list);
    filling_list = cell2mat(filling_list);
end