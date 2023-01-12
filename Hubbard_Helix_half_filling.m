function [gs_mps, gs_energy] = Hubbard_helix_half_filling(N, t, U, trunc, maxiter, tol, vumps_way, starting_name, finalized)
    disp('Code started running');
    P = 1;
    Q = 2;
    [pspace, vspaces, trivspace, fusion_trees] = get_spaces('Hubbard_asymmetric', false, P, Q, 12, 3);
        
    trunc_tot = ~iscell(trunc);
    if trunc_tot
        trunc_way = {'TruncTotalDim', trunc};
    else
        trunc_way = {'TruncBelow', 10^(-trunc{2}), 'TruncDim', trunc{1}};
    end

    mu = 0;

    H = Hubbard_Hamiltonian(t);
    H_one_site = get_hamiltonian('Hubbard_one_site', pspace, trivspace, U);

    mpo = get_mpo(H, N, 'Helix');
    mpo_joint = {mpo mpo};
        
    mpo_joint{1}(1, 1, N+3, 1) = H_one_site;
    mpo_joint{2}(1, 1, N+3, 1) = H_one_site;

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
            mps = UniformMps.randnc(pspace, vspaces(1), pspace, vspaces(2));
        else
            error('check length of vspaces');
        end
    end
    disp('initialization correct');
    
    if trunc_tot
        name = 'Hubbard_Helix_half_filling_N_' + string(N) + '_t_' + string(t) + '_U_' + string(U) + '_trunctotdim_' + string(trunc);
    else
        name = 'Hubbard_Helix_half_filling_N_' + string(N) + '_t_' + string(t) + '_U_' + string(U) + '_truncbond_' + string(trunc{1}) + '_cut_' + string(trunc{2});
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