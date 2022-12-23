function [gs_mps, gs_energy] = XXZ_Cylinder_nosymm(N, delta, trunc, maxiter, vumps_way, stag_h_field, starting_name, final)
    % Does exactly the same as XXZ_Cylinder, but without symmetries.
    doPath
    disp('Code has started running');
    h_field = 0;
    trunc_tot = ~iscell(trunc);
    if trunc_tot
        trunc_way = {'TruncTotalDim', trunc};
        D = trunc;
    else
        trunc_way = {'TruncBelow', 10^(-trunc{2}), 'TruncDim', trunc{1}};
        D = 20;
    end
    
    pspace = ComplexSpace.new(2, false);
    trivspace = ComplexSpace.new(1, false);
    vspace = ComplexSpace.new(D, false);
    
    %%
    H = Tensor([pspace' pspace'], [pspace' pspace']);
    var = zeros(2,2,2,2);
    var(1,1,1,1) = delta;
    var(1,2,1,2) = -delta;
    var(2,1,2,1) = -delta;
    var(2,2,2,2) = delta;
    var(1,2,2,1) = 2;
    var(2,1,1,2) = 2;
    H = fill_tensor(H, var);
    
%%

    [U_base, S, V_base] = tsvd(H, [1 3], [2 4], 'TruncBelow', 1e-12);
    U_base = tpermute(U_base, [2 3 1]);
    dims_U_base = dims(U_base);
    U_new_var = zeros([1, dims_U_base]);
    U_new_var(1,:,:,:) = U_base.var.var;
    cspace = ComplexSpace.new(dims_U_base(2), false);
    U = Tensor([trivspace pspace], [pspace cspace]);
    U = fill_tensor(U, U_new_var);
    L = contract(U, [-1 -2 1 -4], S, [1 -3]);

    V_base = tpermute(V_base, [1 3 2]);
    dims_V_base = dims(V_base);
    V_new_var = zeros([dims_V_base(1:2), 1, dims_V_base(3)]);
    V_new_var(:,:,1,:) = V_base.var.var;
    cspace = ComplexSpace.new(dims_V_base(1), false);
    V = Tensor([cspace pspace], [pspace trivspace]);
    R = fill_tensor(V, V_new_var);
    
    L = MpoTensor(L);
    R = MpoTensor(R);
        
    mpo = MpoTensor.zeros(3, 1, 3, 1);
    mpo(3, 1, 3, 1) = MpoTensor(1);
    mpo(1,1,1,1) = MpoTensor(1);
    mpo(1, 1, 2, 1) = L;
    mpo(2,1,3,1) = R;
    
    H1 = InfJMpo({mpo mpo});
    
    %%
    %load(starting_name);
    mps = UniformMps.randnc(pspace, [vspace vspace]);

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
    if trunc_tot
        naam = 'XXZ_Cylinder_nosymm_vumps_' + string(N) + '_delta_' + string(delta) + '_trunctotdim_' + string(trunc) + '_stagh_' + string(stag_h_field);
    else
        naam = 'XXZ_Cylinder_vumps_nosymm_' + string(N) + '_delta_' + string(delta) + '_truncbond_' + string(trunc{1}) + '_cut_' + string(trunc{2}) + '_stagh_' + string(stag_h_field);
    end
    if vumps_way == 1
        alg1 = Vumps('which', 'smallestreal', 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'doSave', true, 'name', strcat(naam, '.mat'), 'tol', 10^(-6), 'doplot', true);
        [gs_mps, gs_energy] = fixedpoint(alg1, H1, mps);
    elseif vumps_way == 2
        alg2 = Vumps2('which', 'smallestreal', 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'doSave', true, 'trunc', trunc_way, 'name', strcat(naam, '.mat'), 'tol', 10^(-6), 'doplot', true);
        [gs_mps, gs_energy] = fixedpoint(alg2, H1, mps);
    elseif vumps_way == 3
        alg = IDmrg2('dynamical_tols', true, 'which', 'smallestreal', 'trunc', {'TruncBelow', 10^(-trunc{2}), 'TruncDim', trunc{1}}, 'tol', 10^(-5), 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'name', naam, 'doSave', true, 'saveIterations', 1);
        [gs_mps, gs_energy] = fixedpoint(alg, H1, mps);
    else
        maxiter1 = maxiter(1);
        maxiter2 = maxiter(2);
        iterations = maxiter(3);
        alg1 = Vumps('which', 'smallestreal', 'miniter', 1, 'maxiter', maxiter1, 'verbosity', Verbosity.iter, 'doSave', true, 'name', strcat(naam, '.mat'), 'tol', 10^(-6), 'doplot', true);
        alg2 = Vumps2('which', 'smallestreal', 'miniter', 1, 'maxiter', maxiter2, 'verbosity', Verbosity.iter, 'doSave', true, 'trunc', trunc_way, 'name', strcat(naam, '.mat'), 'tol', 10^(-6), 'doplot', true);
        gs_mps = mps;
        for i = 1:iterations
            fprintf('Big iteration %s \n', i);
            [gs_mps, gs_energy] = fixedpoint(alg2, H1, gs_mps);
            [gs_mps, gs_energy] = fixedpoint(alg1, H1, gs_mps);
        end
    end
    %%
    save(strcat(naam, '_final.mat'));
    disp('Done');