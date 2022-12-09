function [gs_mps, gs_energy] = XXZ_Cylinder(N, delta, maxdim, cut, maxiter, size_D, stag_h_field, starting_name, final)
    disp('Code has started running');
    %N = 3;
    %delta = 2;
    %stag_h_field = 0;
    h_field = 0;
    
    charges = U1([1 -1]);
    fusion_trees = FusionTree.new([2 2], charges, false, charges, false, charges, false, charges, false);
    pspace = GradedSpace.new(charges, [1 1], false);
    trivspace = GradedSpace.new(U1(0), 1, false);
    
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
    
    mpo_A = get_mpo(H_A, N, 'Helix');
    mpo_B = get_mpo(H_B, N, 'Helix');

    mpo_joint = {mpo_A mpo_B};
        
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
    naam = 'XXZ_Cylinder_vumps_' + string(N) + '_delta_' + string(delta) + '_size_bond_' + string(size_D) + '_cut_' + string(cut) + '_stagh_' + string(stag_h_field);
    %alg = IDmrg2('dynamical_tols', true, 'which', 'smallestreal', 'trunc', {'TruncBelow', 10^(-cut), 'TruncDim', maxdim}, 'tol', 10^(-5), 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'name', naam, 'doSave', true, 'saveIterations', 1);
    %alg = Vumps2('which', 'smallestreal', 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'dynamical_tols', true, 'doSave', true, 'trunc', {'TruncBelow', 10^(-cut)}, 'name', strcat(naam, '.mat'), 'tol', 10^(-6), 'doplot', true);
    alg = Vumps('which', 'smallestreal', 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'dynamical_tols', true, 'doSave', true, 'name', strcat(naam, '.mat'), 'tol', 10^(-6), 'doplot', true);
    [gs_mps, gs_energy] = fixedpoint(alg, H1, mps);
    
    %%
    save(strcat(naam, '_final.mat'));
    disp('Done');
    %magn = get_magnetisation(gs_mps, pspace, trivspace, 2, true);
    
    %[prob_up, prob_down] = get_probabilities(gs_mps, pspace, trivspace, 2);
