function gs_mps = Hubbard_1D_external_field_half_filling(t, U, maxdim, cut, maxiter, redefined, previous, varargin)
    doPath
    disp('Code started running');
    [pspace, vspace, trivspace, fusion_trees] = get_spaces('Hubbard', false, 1, 1, 12, 3);
        
    % For atomic separation = 2å of 1D hydrogen chain: 
    % t = 1.5 en abs(U/t) = 6
    %t = 1.5;
    %U = 9;
    mu = 0;
    h_field = 0;
    N = 0;

    %H = get_hamiltonian('Hubbard_external_field', fusion_trees, pspace, t, 0, 0, 0);
    H = get_hamiltonian('Hubbard_two_site', fusion_trees, pspace, t, mu, true);
    %H = tpermute(H, [3 4 1 2], [2 2]);
    if redefined
        H_one_site = get_hamiltonian('Hubbard_one_site_redefined', pspace, trivspace, U);
    else
        H_one_site = get_hamiltonian('Hubbard_one_site', pspace, trivspace, U);
    end
    mpo = get_mpo(H, 0, 'Helix');
    mpo_joint = {mpo mpo};
        
    mpo_joint{1}(1, 1, N+3, 1) = H_one_site;
    mpo_joint{2}(1, 1, N+3, 1) = H_one_site;

    H1 = InfJMpo(mpo_joint);
    if previous
        mps = varargin{1};
    else
        mps = UniformMps.randnc(pspace, vspace);
    end
    disp('initialization correct');
    
    %alg = Vumps('which', 'smallestreal', 'miniter', 1, 'maxiter', 10, 'verbosity', Verbosity.iter);
    %[gs_mps, gs_energy] = fixedpoint(alg, H1, mps);
    %cut = 4;
    naam = 'Hubbard_1D_half_filling_idmrg2_t_' + string(t) + '_U_' + string(U) + '_maxdim_' +string(maxdim) + '_cut_' + string(cut) + '_redef_' + string(redefined);
    %{'TruncBelow', 10^(-cut)} can be used.
    %alg = IDmrg2('dynamical_tols', true, 'which', 'smallestreal', 'trunc', {'TruncBelow', 10^(-cut), 'TruncDim', maxdim}, 'tol', 10^(-6), 'maxiter', maxdim, 'verbosity', Verbosity.iter, 'name', naam, 'doSave', true, 'saveIterations', 1);
    alg = Vumps2('which', 'smallestreal', 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'trunc', {'TruncBelow', 10^(-cut), 'TruncDim', maxdim}, 'dynamical_tols', true, 'doSave', true, 'name', naam, 'tol', 10^(-6), 'doplot', true);
    
    [gs_mps, gs_energy] = fixedpoint(alg, H1, mps);
    
    save(naam + '_final');
    disp('Done, Hooray!');

end