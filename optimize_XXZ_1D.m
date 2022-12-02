function gs_mps = optimize_XXZ_1D(power)
    doPath
    disp('Code started running');
    N = 0;
    delta = 2;
    stag_h_field = 0;
    shear = 0;
    
    f = 1.1^power;

    charges = U1([1 -1]);
    pspace = GradedSpace.new(charges, [1 1], false);
    trivspace = GradedSpace.new(U1(0), 1, false);
    start_1 = [8 41 94 117 80 25 3];
    start_2 = [2 20 70 114 102 49 11];
    dim_1 = round(f*start_1);
    dim_2 = round(f*start_2);

    vspace1 = GradedSpace.new(U1([-5 -3 -1 1 3 5 7]), dim_1, false);
    vspace2 = GradedSpace.new(U1([-6 -4 -2 0 2 4 6]), dim_2, false);
    
    H = get_hamiltonian('XXZ', delta, pspace);
    H = tpermute(H, [3 4 1 2], [2 2]);
    H_one_site = get_hamiltonian('one_site_XXZ', pspace, trivspace, stag_h_field, shear);
    H_one_site_A = H_one_site{1};
    H_one_site_B = H_one_site{2};
    
    mpo = get_mpo(H, N);
    mpo_joint = {mpo mpo};
        
    mpo_joint{1}(1, 1, N+3, 1) = tpermute(H_one_site_A, [2 1 4 3], [2 2]);
    mpo_joint{2}(1, 1, N+3, 1) = tpermute(H_one_site_B, [2 1 4 3], [2 2]);
    
    H1 = InfJMpo(mpo_joint);

    mps = UniformMps.randnc(pspace, [vspace1 vspace2]);
    
    disp('initialization correct');

    alg = Vumps('which', 'smallestreal', 'maxiter', 200, 'verbosity', Verbosity.iter, 'dynamical_tols', false, 'doSave', true, 'name', 'optimize_XXZ_1D_' + string(power));
    [gs_mps, gs_energy] = fixedpoint(alg, H1, mps);

    save('optimize_XXZ_1D_' + string(power) + '_final');
    disp('Done');
end