function [gs_mps, gs_energy] = Hubbard_cylinder_half_filling(N, t, U, trunc, maxiter, tol, vumps_way, starting_name, finalized)
    disp('Code started running');
    P = 1;
    Q = 1;
    [pspace, vspaces, trivspace, prodspace, fusion_trees] = get_spaces('Hubbard', false, P, Q, 12, 3);
%    [pspace, vspaces, trivspace] = get_spaces('Heisenberg XXZ');
        
    trunc_tot = ~iscell(trunc);
    if trunc_tot
        trunc_way = {'TruncTotalDim', trunc};
    else
        trunc_way = {'TruncBelow', 10^(-trunc{2}), 'TruncDim', trunc{1}};
    end

    mu = 0;

    H = Hubbard_Hamiltonian(t, P, Q);
    H_one_site = get_hamiltonian('Hubbard_one_site', pspace, trivspace, U);

%    H = get_hamiltonian('XXZ', 1, pspace);
%    H_one_site = get_hamiltonian('one_site_XXZ', pspace, trivspace, 0, 0);


    mpo_way = 'FullCylinder';
    mpo_joint = get_mpo(H, N, mpo_way);        

    if strcmp('FullCylinder', mpo_way)
        one_site_place = N+4;
    elseif strcmp('FullCylinder_inefficient', mpo_way)
        one_site_place = 2*N+2;
    elseif strcmp('1D', mpo_way)
        one_site_place = 3;
    else
        error('Mpo_way not implemented')
    end

    for i = 1: max(2*N,2)
        mpo_joint{i}(1, 1, one_site_place, 1) = H_one_site;
    end

    H1 = InfJMpo(mpo_joint);
    if finalized == 2
        load(starting_name, 'gs_mps');
        mps = gs_mps;
    elseif finalized == 1
        load(starting_name, 'mps');
        mps = canonicalize(mps, 'Order', 'rl');
    else
        args = cell(2, max(2*N,2));
        for i = 1:max(2*N,2)
            args{1,i} = pspace;
            args{2,i} = vspaces((mod(i-1,2)+1));
        end
        mps = UniformMps.randnc(args{:});

    end
    disp('initialization correct');
    
    if trunc_tot
        name = 'Hubbard_FullCylinder_half_filling_N_' + string(N) + '_t_' + string(t) + '_U_' + string(U) + '_trunctotdim_' + string(trunc);
    else
        name = 'Hubbard_FullCylinder_half_filling_N_' + string(N) + '_t_' + string(t) + '_U_' + string(U) + '_truncbond_' + string(trunc{1}) + '_cut_' + string(trunc{2});
    end

    if vumps_way == 1
        alg1 = Vumps('which', 'smallestreal', 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'doSave', true, 'name', strcat(name, '.mat'), 'tol', 10^(-6), 'doplot', true);
        [gs_mps, gs_energy] = fixedpoint(alg1, H1, mps);
    elseif vumps_way == 2
        alg2 = Vumps2('which', 'smallestreal', 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'doSave', true, 'trunc', trunc_way, 'name', strcat(name, '.mat'), 'tol', 10^(-6), 'doplot', true);
        [gs_mps, gs_energy] = fixedpoint(alg2, H1, mps);
    elseif vumps_way == 3
        cut = 5;
        maxdim = 70;
        alg = IDmrg2('dynamical_tols', true, 'which', 'smallestreal', 'trunc', {'TruncBelow', 10^(-cut), 'TruncDim', maxdim}, 'tol', 10^(-tol), 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'name', strcat(name, '.mat'), 'doSave', true, 'saveIterations', 1);
        [gs_mps, gs_energy] = fixedpoint(alg, H1, mps);
    else
        [gs_mps, gs_energy, eta] = doVumps(H1, mps, name, maxiter, tol, trunc_way);
    end
    
    save(strcat(name, '_final.mat'));
    disp('Done, Hooray!');

end