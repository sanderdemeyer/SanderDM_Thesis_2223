function [gs_mps, gs_energy, eta] = Heisenberg_XXZ_1D(delta, D1, vumps_way, maxiter, trunc, tol, finalized, starting_name, kwargs)
    arguments
        delta
        D1
        vumps_way
        maxiter
        trunc
        tol
        finalized
        starting_name
        kwargs.len = 2
        kwargs.convention = 'conventional'
    end

    charges = U1([1 -1]);
    fusion_trees = FusionTree.new([2 2], charges, false, charges, false, charges, false, charges, false);
    pspace = GradedSpace.new(charges, [1 1], false);
    trivspace = GradedSpace.new(U1(0), 1, false);
    prodspace = 0;
    vspace1 = GradedSpace.new(U1([-1 1]), [D1 D1], false);
    vspace2 = GradedSpace.new(U1([2, 0, -2]), [D1 D1 D1], false);
    vspace = [vspace1 vspace2];

    H = OLD_get_hamiltonian('XXZ', delta, pspace);
    H = tpermute(H, [3 4 2 1]);

    tens_one_site = Tensor([trivspace pspace], [pspace trivspace]);
    tblocks = num2cell([0 0 0 0]);
    
    H0 = fill_tensor(tens_one_site, tblocks);
    
    mpo = get_mpo_1D(H, H0, 'len', kwargs.len, 'convention', kwargs.convention);
    H1 = InfJMpo(mpo);

    if finalized == 2
        load(starting_name, 'gs_mps');
        mps = gs_mps;
    elseif finalized == 1
        load(starting_name, 'mps')
        mps = Canonicalize(mps, 'Order', 'rl');
    else
        mps = UniformMps.randnc(pspace, vspace1, pspace, vspace2);
    end

    name = 'Heisenberg_XXZ_1D_delta_' + string(delta);

    [gs_mps, gs_energy, eta] = doVumps(H1, mps, vumps_way, maxiter, trunc, tol, name);

end