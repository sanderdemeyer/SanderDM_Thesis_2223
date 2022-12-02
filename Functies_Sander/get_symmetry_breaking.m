function magn = get_symmetry_breaking(N, stag_h_field)
    delta = 1;
    shear = 0;
    charges = U1([1 -1]);
    fusion_trees = FusionTree.new([2 2], charges, false, charges, false, charges, false, charges, false);
    pspace = GradedSpace.new(charges, [1 1], false);
    %%
    H = get_hamiltonian('XXZ', delta, pspace);
    H_one_site = get_hamiltonian('one_site_XXZ', pspace, stag_h_field, shear);
    H_one_site_A = H_one_site{1};
    H_one_site_B = H_one_site{2};
    
    mpo = get_mpo(H, N);
    mpo_joint = {mpo mpo};
    
    mpo_joint{1}(1, 1, N+3, 1) = tpermute(H_one_site_A, [4 3 2 1], [2 2]);
    mpo_joint{2}(1, 1, N+3, 1) = tpermute(H_one_site_B, [4 3 2 1], [2 2]);
    
    H1 = InfJMpo(mpo_joint);
    
    %%
    D = 150;
    vspace1 = GradedSpace.new(U1([-1 1]), [D D], false);
    vspace2 = GradedSpace.new(U1([2, 0, -2]), [D/2 D/2 D/2], false);
    mps = UniformMps.randnc([pspace' pspace'], [vspace1 vspace2]);
    
    %%
    alg = Vumps('which', 'smallestreal', 'maxiter', 40, 'verbosity', Verbosity.iter);
    [gs_mps, gs_energy] = fixedpoint(alg, H1, mps);
    
    %%
    trivspace = GradedSpace.new(U1(0), 1, false);
    
    magn = get_magnetisation(gs_mps, pspace, trivspace, 2, true);
    return
end
