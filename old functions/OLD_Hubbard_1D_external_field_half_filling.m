function [gs_mps, gs_energy] = Hubbard_1D_external_field_half_filling(t, U, trunc, maxiter, tol, vumps_way, starting_name, finalized)
    disp('Code started running');
    P = 2;
    Q = 5;

    %[pspace, vspaces, trivspace, fusion_trees] = get_spaces('Hubbard_asymmetric_again', P, Q);
    D1 = 1;
    D2 = 2;
    [pspace, vspaces, trivspace] = get_spaces_Hubbard_asymmetric(P, Q, D1, D2);
        
    trunc_tot = ~iscell(trunc);
    if trunc_tot
        trunc_way = {'TruncTotalDim', trunc};
    else
        trunc_way = {'TruncBelow', 10^(-trunc{2}), 'TruncDim', trunc{1}};
    end

    % For atomic separation = 2å of 1D hydrogen chain: 
    % t = 1.5 en abs(U/t) = 6
    %t = 1.5;
    %U = 9;
    mu = 0;
    h_field = 0;
    N = 0;

    H = Hubbard_Hamiltonian(t, P, Q);
    H_one_site = get_hamiltonian('Hubbard_one_site', pspace, trivspace, U);

    mpo = get_mpo(H, 0, 'Helix');
    %mpo_joint = {mpo mpo};
    mpo_joint = cell(1, length(vspaces));
    mpo_joint(:) = {mpo};
    %mpo_joint = {mpo mpo mpo mpo mpo mpo};
    
    for i = 1: length(vspaces)
        mpo_joint{i}(1, 1, N+3, 1) = H_one_site;
    end

    %{
    mpo_joint{1}(1, 1, N+3, 1) = H_one_site;
    mpo_joint{2}(1, 1, N+3, 1) = H_one_site;
    mpo_joint{3}(1, 1, N+3, 1) = H_one_site;
    mpo_joint{4}(1, 1, N+3, 1) = H_one_site;
    mpo_joint{5}(1, 1, N+3, 1) = H_one_site;
    mpo_joint{6}(1, 1, N+3, 1) = H_one_site;
    %}
    H1 = InfJMpo(mpo_joint);
    if finalized == 2
        load(starting_name, 'gs_mps');
        mps = gs_mps;
    elseif finalized == 1
        load(starting_name, 'mps');
        mps = canonicalize(mps, 'Order', 'rl');
    else
        args = cell(2, length(vspaces));
        for i = 1:length(vspaces)
            args{1,i} = pspace;
            args{2,i} = vspaces(i);
        end

        mps = UniformMps.randnc(args{:});
        %{
        %mps = UniformMps.randnc(pspace, {vspace1 vspace2});
        if length(vspaces) == 1
            mps = UniformMps.randnc(pspace, vspaces);
        elseif length(vspaces) == 2
            mps = UniformMps.randnc(pspace, vspaces(1), pspace, vspaces(2));
        elseif length(vspaces) == 4
            mps = UniformMps.randnc(pspace, vspaces(1), pspace, vspaces(2), pspace, vspaces(3), pspace, vspaces(4));
        elseif length(vspaces) == 6
            mps = UniformMps.randnc(pspace, vspaces(1), pspace, vspaces(2), pspace, vspaces(3), pspace, vspaces(4), pspace, vspaces(5), pspace, vspaces(6));
        else
            error('check length of vspaces');
        end
        %}
    end
    disp('initialization correct');
    
    if trunc_tot
        name = 'Hubbard_1D_half_filling_t_' + string(t) + '_U_' + string(U) + '_trunctotdim_' + string(trunc) + '_h_' + string(h_field) + '_redef_' + string(redefined);
    else
        name = 'Hubbard_1D_half_filling_t_' + string(t) + '_U_' + string(U) + '_truncbond_' + string(trunc{1}) + '_cut_' + string(trunc{2}) + '_h_' + string(h_field) + '_redef_' + string(redefined);
    end

    if vumps_way == 1
        alg1 = Vumps('which', 'smallestreal', 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'doSave', true, 'name', strcat(name, '.mat'), 'tol', 10^(-tol), 'doplot', true);
        [gs_mps, gs_energy] = fixedpoint(alg1, H1, mps);
    elseif vumps_way == 2
        alg2 = Vumps2('which', 'smallestreal', 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'doSave', true, 'trunc', trunc_way, 'name', strcat(name, '.mat'), 'tol', 10^(-tol), 'doplot', true);
        [gs_mps, gs_energy] = fixedpoint(alg2, H1, mps);
    else
        [gs_mps, gs_energy, eta] = doVumps(H1, mps, name, maxiter, tol, trunc_way);
    end
    
    save(strcat(name, '_final.mat'));
    disp('Done, Hooray!');

end