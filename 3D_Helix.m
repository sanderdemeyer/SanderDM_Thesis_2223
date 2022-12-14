function [gs_mps, gs_energy] = 3D_Helix(N, delta, trunc, maxiter, vumps_way, stag_h_field, starting_name, final)
    doPath
    disp('Code has started running');
    %N = 3;
    %delta = 2;
    %stag_h_field = 0;
    h_field = 0;
    trunc_tot = ~iscell(trunc);
    if trunc_tot
        trunc_way = {'TruncTotalDim', trunc};
    else
        trunc_way = {'TruncBelow', 10^(-trunc{2}), 'TruncDim', trunc{1}};
    end
    
    charges = U1([1 -1]);
    fusion_trees = FusionTree.new([2 2], charges, false, charges, false, charges, false, charges, false);
    pspace = GradedSpace.new(charges, [1 1], false);
    trivspace = GradedSpace.new(U1(0), 1, false);
    size_D = 2;
    %D = 12;

    if size_D == 3
        D1 = [1 6 30 50 70 50 30 6 1];
        D2 = [1 8 20 40 65 65 40 20 8 1];

        vspace1 = GradedSpace.new(U1([-7 -5 -3 -1 1 3 5 7 9]), D1, false);
        vspace2 = GradedSpace.new(U1([-8 -6 -4 -2 0 2 4 6 8 10]), D2, false);
    elseif size_D == 2
        D1 = [2 16 39 52 39 16 2];
        D2 = [1 9 28 49 49 28 9 1];

        vspace1 = GradedSpace.new(U1([-5 -3 -1 1 3 5 7]), D1, false);
        vspace2 = GradedSpace.new(U1([-6 -4 -2 0 2 4 6 8]), D2, false);
    elseif size_D == 1
        D1 = [12 25 40 25 12];
        D2 = [6 16 32 32 16 6];

        vspace1 = GradedSpace.new(U1([-3 -1 1 3 5]), D1, false);
        vspace2 = GradedSpace.new(U1([-4 -2 0 2 4 6]), D2, false);

    elseif size_D == 0
        D1 = [6 18 30 18 6];
        D2 = [2 12 26 26 12 2];

        vspace1 = GradedSpace.new(U1([-3 -1 1 3 5]), D1, false);
        vspace2 = GradedSpace.new(U1([-4 -2 0 2 4 6]), D2, false);
    end
    %vspace1 = GradedSpace.new(U1([-1 1]), [D D], false);
    %vspace2 = GradedSpace.new(U1([2, 0, -2]), [D/2 D/2 D/2], false);
    
    %%
    H_A = get_hamiltonian('XXZ_stag', pspace, trivspace, stag_h_field, 'AB', delta);
    H_B = get_hamiltonian('XXZ_stag', pspace, trivspace, stag_h_field, 'BA', delta);
    H_one_site = get_hamiltonian('one_site_XXZ', pspace, trivspace, 0, 0);
    H_one_site_A = H_one_site{1};
    H_one_site_B = H_one_site{2};
    
    %mpo_A = get_mpo(H_A, N, '3D_Helix');
    %mpo_B = get_mpo(H_B, N, '3D_Helix');
    mpo = get_mpo(H_A, [P Q], '3D_Helix');

    %mpo_joint = {mpo_A mpo_B};
    mpo_joint = {mpo mpo}; 

    mpo_joint{1}(1, 1, N+3, 1) = tpermute(H_one_site_A, [2 1 4 3], [2 2]);
    mpo_joint{2}(1, 1, N+3, 1) = tpermute(H_one_site_B, [2 1 4 3], [2 2]);
    
    H1 = InfJMpo(mpo_joint);
    
    %%
    %load(starting_name);
    mps = UniformMps.randnc(pspace, [vspace1 vspace2]);

    %{
    if ~final
        mps = Canonicalize(mps, 'Order', 'rl');
    else
        mps = gs_mps;
    end
    %}
    disp('initialization correct');
    %%
    %naam = strcat(starting_name, 'VUMPS')
    if trunc_tot
        naam = 'XXZ_Cylinder_vumps_' + string(N) + '_delta_' + string(delta) + '_trunctotdim_' + string(trunc) + '_stagh_' + string(stag_h_field);
    else
        naam = 'XXZ_Cylinder_vumps_' + string(N) + '_delta_' + string(delta) + '_truncbond_' + string(trunc{1}) + '_cut_' + string(trunc{2}) + '_stagh_' + string(stag_h_field);
    end
    %alg = IDmrg2('dynamical_tols', true, 'which', 'smallestreal', 'trunc', {'TruncBelow', 10^(-cut), 'TruncDim', maxdim}, 'tol', 10^(-5), 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'name', naam, 'doSave', true, 'saveIterations', 1);
    if vumps_way == 1
        alg1 = Vumps('which', 'smallestreal', 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'doSave', true, 'name', strcat(naam, '.mat'), 'tol', 10^(-6), 'doplot', true);
        [gs_mps, gs_energy] = fixedpoint(alg1, H1, mps);
    elseif vumps_way == 2
        alg2 = Vumps2('which', 'smallestreal', 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'doSave', true, 'trunc', trunc_way, 'name', strcat(naam, '.mat'), 'tol', 10^(-6), 'doplot', true);
        [gs_mps, gs_energy] = fixedpoint(alg2, H1, mps);
    else
        maxiter1 = maxiter(1);
        maxiter2 = maxiter(2);
        iterations = maxiter(3);
        alg1 = Vumps('which', 'smallestreal', 'miniter', 1, 'maxiter', maxiter1, 'verbosity', Verbosity.iter, 'doSave', true, 'name', strcat(naam, '.mat'), 'tol', 10^(-6), 'doplot', true);
        alg2 = Vumps2('which', 'smallestreal', 'miniter', 1, 'maxiter', maxiter2, 'verbosity', Verbosity.iter, 'doSave', true, 'trunc', trunc_way, 'name', strcat(naam, '.mat'), 'tol', 10^(-6), 'doplot', true);
        gs_mps = mps;
        for i = 1:iterations
            fprintf('Big iteration %s \n', i);
            [gs_mps, gs_energy] = fixedpoint(alg2, H1, gs_mps);
            [gs_mps, gs_energy] = fixedpoint(alg1, H1, gs_mps);
        end
    end
    %%
    save(strcat(naam, '_final.mat'));
    disp('Done');
    %magn = get_magnetisation(gs_mps, pspace, trivspace, 2, true);
    
    %[prob_up, prob_down] = get_probabilities(gs_mps, pspace, trivspace, 2);
