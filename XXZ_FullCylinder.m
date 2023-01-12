function [gs_mps, gs_energy] = XXZ_FullCylinder(N, trunc, maxiter, vumps_way, stag_h_field, starting_name, finalized)
    warning('You should probably use the file Heisenberg_XXZ_cylinder');
    doPath
    disp('Code started running');
    %N = 3;
    delta = 1;
    %stag_h_field = 0;
    
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
    
    %D1 = [2 16 39 52 39 16 2];
    %D2 = [1 9 28 49 49 28 9 1];
    D1 = [12 12];
    D2 = [8 8 8];

    %vspace1 = GradedSpace.new(U1([-5 -3 -1 1 3 5 7]), D1, false);
    %vspace2 = GradedSpace.new(U1([-6 -4 -2 0 2 4 6 8]), D2, false);
    vspace1 = GradedSpace.new(U1([-1 1]), D1, false);
    vspace2 = GradedSpace.new(U1([-2 0 2]), D2, false);
    

    %%
    H = get_hamiltonian('XXZ', 1, pspace);
    H_one_site = get_hamiltonian('one_site_XXZ', pspace, trivspace, stag_h_field, 0);
    H_one_site_A = H_one_site{1};
    H_one_site_B = H_one_site{2};
    H_one_site_A = tpermute(H_one_site_A, [2 1 4 3], [2 2]);
    H_one_site_B = tpermute(H_one_site_B, [2 1 4 3], [2 2]);

    mpo_joint = get_mpo(H, N, 'FullCylinder_inefficient');
    %mpo_joint = get_mpo(H, N, 'Helix');

    sz = 2*N+2;
    for n = 1:2*N
        if ((mod(n,2) == 0) - 1/2)*((n < N + 1) - 1/2) == 1/4
            mpo_joint{n}(1,1,sz,1) = H_one_site_B;
        else
            mpo_joint{n}(1,1,sz,1) = H_one_site_A;
        end
    end
    H1 = InfJMpo(mpo_joint);
    
    %%
    %mps = UniformMps.randnc(pspace, [vspace1 vspace2]);

    pspaces = repmat(pspace, 1, 2*N);
    vspaces = repmat([vspace1 vspace2], 1, N);
    
    if finalized == 2
        load(starting_name, 'gs_mps');
        mps = gs_mps;
    elseif finalized == 1
        load(starting_name, 'mps');
        mps = Canonicalize(mps, 'Order', 'rl');
    else
        mps = UniformMps.randnc(pspaces , vspaces);
    end
    disp('initialization correct');
    
    %%
    if trunc_tot
        naam = 'XXZ_FullCylinder_vumps_' + string(N) + '_delta_' + string(delta) + '_trunctotdim_' + string(trunc) + '_stagh_' + string(stag_h_field);
    else
        naam = 'XXZ_FullCylinder_vumps_' + string(N) + '_delta_' + string(delta) + '_truncbond_' + string(trunc{1}) + '_cut_' + string(trunc{2}) + '_stagh_' + string(stag_h_field);
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
            fprintf('Big iteration %d \n', i);
            [gs_mps, gs_energy, ~, ~, eta] = fixedpoint(alg2, H1, gs_mps);
            [gs_mps, gs_energy, ~, ~, eta] = fixedpoint(alg1, H1, gs_mps);
            if eta < 10^(-5)
                save(strcat(naam, '_final.mat'));
                disp('Done');
                return
            end
        end
    end



    %%
    save(strcat(naam, '_final.mat'));
    disp('Done');
    %magn = get_magnetisation(gs_mps, pspace, trivspace, 2, true);
    
    %[prob_up, prob_down] = get_probabilities(gs_mps, pspace, trivspace, 2);
