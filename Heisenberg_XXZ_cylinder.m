function [gs_mps, gs_energy] = Heisenberg_XXZ_cylinder(N, delta, trunc, maxiter, tol, vumps_way, starting_name, finalized)
    disp('Code started running');
    charges = U1([1 -1]);
    fusion_trees = FusionTree.new([2 2], charges, false, charges, false, charges, false, charges, false);
    pspace = GradedSpace.new(charges, [1 1], false);
    trivspace = GradedSpace.new(U1(0), 1, false);
        
    D1 = [2 16 39 52 39 16 2];
    D2 = [1 9 28 49 49 28 9 1];

    vspace1 = GradedSpace.new(U1([-5 -3 -1 1 3 5 7]), D1, false);
    vspace2 = GradedSpace.new(U1([-6 -4 -2 0 2 4 6 8]), D2, false);
    vspaces = [vspace1 vspace2];

    trunc_tot = ~iscell(trunc);
    if trunc_tot
        trunc_way = {'TruncTotalDim', trunc};
    else
        trunc_way = {'TruncBelow', 10^(-trunc{2}), 'TruncDim', trunc{1}};
    end

    stag_h_field = 0;
    H = get_hamiltonian('XXZ', 1, pspace);
    H_one_site = get_hamiltonian('one_site_XXZ', pspace, trivspace, stag_h_field, 0);

    mpo_way = 'FullCylinder';
    mpo_joint = get_mpo(H, N, mpo_way);        

    if strcmp('FullCylinder', mpo_way)
        place = N+4;
    elseif strcmp('FullCylinder_inefficient', mpo_way)
        place = 2*N+2;
    else
        error('Mpo_way not implemented')
    end

    for i = 1: 2*N
        mpo_joint{i}(1, 1, place, 1) = H_one_site;
    end

    H1 = InfJMpo(mpo_joint);
    if finalized == 2
        load(starting_name, 'gs_mps');
        mps = gs_mps;
    elseif finalized == 1
        load(starting_name, 'mps');
        mps = canonicalize(mps, 'Order', 'rl');
    else
        if length(vspaces) == 1
            mps = UniformMps.randnc(pspace, vspaces);
        elseif length(vspaces) == 2
            if N == 2
                mps = UniformMps.randnc(pspace, vspaces(1), pspace, vspaces(2), pspace, vspaces(1), pspace, vspaces(2));                
            elseif N == 4
                mps = UniformMps.randnc(pspace, vspaces(1), pspace, vspaces(2), pspace, vspaces(1), pspace, vspaces(2), pspace, vspaces(1), pspace, vspaces(2), pspace, vspaces(1), pspace, vspaces(2));
            elseif N == 6
                mps = UniformMps.randnc(pspace, vspaces(1), pspace, vspaces(2), pspace, vspaces(1), pspace, vspaces(2), pspace, vspaces(1), pspace, vspaces(2), pspace, vspaces(1), pspace, vspaces(2), pspace, vspaces(1), pspace, vspaces(2), pspace, vspaces(1), pspace, vspaces(2));
            else
                error('TBA');
            end
        else
            error('check length of vspaces');
        end
    end
    disp('initialization correct');
    
    if trunc_tot
        name = 'Heisenberg_cylinder_delta_' + string(delta) + '_trunctotdim_' + string(trunc);
    else
        name = 'Heisenberg_cylinder_delta' + string(delta) + '_truncbond_' + string(trunc{1}) + '_cut_' + string(trunc{2});
    end

    if vumps_way == 1
        alg1 = Vumps('which', 'smallestreal', 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'doSave', true, 'name', strcat(name, '.mat'), 'tol', 10^(-6), 'doplot', true);
        [gs_mps, gs_energy] = fixedpoint(alg1, H1, mps);
    elseif vumps_way == 2
        alg2 = Vumps2('which', 'smallestreal', 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'doSave', true, 'trunc', trunc_way, 'name', strcat(name, '.mat'), 'tol', 10^(-6), 'doplot', true);
        [gs_mps, gs_energy] = fixedpoint(alg2, H1, mps);
    else
        [gs_mps, gs_energy, eta] = doVumps(H1, mps, name, maxiter, tol, trunc_way);
    end
    
    save(strcat(name, '_final.mat'));
    disp('Done, Hooray!');

end