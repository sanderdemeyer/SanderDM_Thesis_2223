function gs_mps = FullCylinder_VUMPS(N, cut, maxiter, stag_h_field, starting_name, finalized)
    doPath
    load(starting_name);
    disp('Code started running');
    %N = 3;
    delta = 1;
    %stag_h_field = 0;
    
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
        mps = gs_mps;
    elseif finalized == 1
        mps = canonicalize(mps, 'Order', 'rl');
    else
        mps = UniformMps.randnc(pspaces , vspaces);
    end
    disp('initialization correct');
    %%
    %cut = 5;
    %naam = strcat(starting_name, 'VUMPS')
    naam = 'XXX_FullCylinder_VUMPS_' + string(N) + '_cut_' + string(cut) + '_stagh_' + string(stag_h_field);
    %alg = IDmrg2('dynamical_tols', true, 'which', 'smallestreal', 'trunc', {'TruncBelow', 10^(-cut), 'TruncDim', maxdim}, 'tol', 10^(-5), 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'name', strcat(naam, '.mat'), 'doSave', true, 'saveIterations', 1);
    %alg = Vumps('which', 'smallestreal', 'maxiter', 2, 'verbosity', Verbosity.iter);
    alg = Vumps2('which', 'smallestreal', 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'trunc', {'TruncBelow', 10^(-cut)}, 'doSave', true, 'name', strcat(naam, '.mat'), 'tol', 10^(-6));
    [gs_mps, gs_energy] = fixedpoint(alg, H1, mps);
    
    %%
    save(strcat(naam, '_final.mat'));
    disp('Done');
    %magn = get_magnetisation(gs_mps, pspace, trivspace, 2, true);
    
    %[prob_up, prob_down] = get_probabilities(gs_mps, pspace, trivspace, 2);
