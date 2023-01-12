function XXZ_Cylinder_VUMPS(N, maxdim, cut, maxiter, stag_h_field, starting_name, final)
    disp('Code started running');
    %N = 3;
    delta = 1;
    %stag_h_field = 0;
    h_field = 0;
    
    charges = U1([1 -1]);
    fusion_trees = FusionTree.new([2 2], charges, false, charges, false, charges, false, charges, false);
    pspace = GradedSpace.new(charges, [1 1], false);
    trivspace = GradedSpace.new(U1(0), 1, false);
    
    D = 12;
    vspace1 = GradedSpace.new(U1([-1 1]), [D D], false);
    vspace2 = GradedSpace.new(U1([2, 0, -2]), [D/2 D/2 D/2], false);
    
    %%
    H_A = get_hamiltonian('XXX_stag', pspace, trivspace, stag_h_field, 'AB');
    H_B = get_hamiltonian('XXX_stag', pspace, trivspace, stag_h_field, 'BA');
    H_one_site = get_hamiltonian('one_site_XXZ', pspace, trivspace, 0, 0);
    H_one_site_A = H_one_site{1};
    H_one_site_B = H_one_site{2};
    
    mpo_A = get_mpo(H_A, N);
    mpo_B = get_mpo(H_B, N);

    mpo_joint = {mpo_A mpo_B};
    
    %mpo_joint{1}(1, 1, N+3, 1) = tpermute(H_one_site_A, [4 3 2 1], [2 2]);
    %mpo_joint{2}(1, 1, N+3, 1) = tpermute(H_one_site_B, [4 3 2 1], [2 2]);
    
    %mpo_joint{1}(1, 1, N+3, 1) = tpermute(H_one_site_A, [4 1 2 3], [2 2]);
    %mpo_joint{2}(1, 1, N+3, 1) = tpermute(H_one_site_B, [4 1 2 3], [2 2]);
    
    mpo_joint{1}(1, 1, N+3, 1) = tpermute(H_one_site_A, [2 1 4 3], [2 2]);
    mpo_joint{2}(1, 1, N+3, 1) = tpermute(H_one_site_B, [2 1 4 3], [2 2]);
    
    H1 = InfJMpo(mpo_joint);
    
    %%
    load(starting_name);
    %mps = UniformMps.randnc(pspace, [vspace1 vspace2]);
    %load('XXZ_Cylinder_N_3_cut_5_final.mat');

    if ~final
        mps = Canonicalize(mps, 'Order', 'rl');
    else
        mps = gs_mps;
    end
    disp('initialization correct');
    %%
    %cut = 5;
    %naam = strcat(starting_name, 'VUMPS')
    naam = 'XXZ_Cylinder_idmrg2_' + string(N) + '_maxdim_' + string(maxdim) + '_cut_' + string(cut) + '_stagh_' + string(stag_h_field);
    %alg = IDmrg2('dynamical_tols', true, 'which', 'smallestreal', 'trunc', {'TruncBelow', 10^(-cut), 'TruncDim', maxdim}, 'tol', 10^(-5), 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'name', naam, 'doSave', true, 'saveIterations', 1);
    %alg = Vumps('which', 'smallestreal', 'maxiter', 2, 'verbosity', Verbosity.iter);
    alg = Vumps('which', 'smallestreal', 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'dynamical_tols', false, 'doSave', true, 'name', strcat(naam, '.mat'), 'tol', 10^(-6));
    [gs_mps, gs_energy] = fixedpoint(alg, H1, mps);
    
    
    %%
    save(strcat(naam, '_final.mat'));
    disp('Done');
    %magn = get_magnetisation(gs_mps, pspace, trivspace, 2, true);
    
    %[prob_up, prob_down] = get_probabilities(gs_mps, pspace, trivspace, 2);
