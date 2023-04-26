function [gs_mps, gs_energy] = Hubbard_cylinder_threeband(N, model, P, Q, rungs, trunc, maxiter, tol, vumps_way, starting_name, finalized, kwargs)
    % N is the number of 3-site unit cells that is used. The radius is thus equal to 3N
    arguments
        N
        model
        P
        Q
        rungs
        trunc
        maxiter
        tol
        vumps_way
        starting_name
        finalized
        kwargs.symmetries = 'U1_SU2' % Symmetry of the charge and spin sector. Fermionic parity is always implemented.
        kwargs.convention = 'conventional'
        kwargs.trunc_method = 'TruncTotalDim'
    end
    warning('Minus signs not fixed!!! See From Real Materials to Model Hamiltonians With Density Matrix Downfolding');
    if strcmp(model, 'Hg1')
        disp('Model is HgBa2CuO4_Hirayama_2018');
        scale = 1.257;
        param.t_dp = 1;
        param.t_pp = 0.751/scale;
        param.delta_dp = 2.416/scale;
        param.U_dd = 8.84/scale;
        param.U_pp = 5.31/scale;
        param.V_dd = 0.8/scale;
        param.V_pp = 1.21/scale;
        param.V_dp = 1.99/scale;
    elseif strcmp(model, 'La1')
        disp('Model is La2CuO4_Hirayama_2018');
        scale = 1.369;
        param.t_dp = 1;
        param.t_pp = 0.754/scale;
        param.delta_dp = 3.699/scale;
        param.U_dd = 9.61/scale;
        param.U_pp = 6.13/scale;
        param.V_dd = 1.51/scale;
        param.V_pp = 1.86/scale;
        param.V_dp = 2.68/scale;
    else
        error('model %s not implemented', model);
    end

    assert(param.t_dp == 1, 'work with t_dp = 1 as normalisation. The real value of t_dp can be used as energy scale')
    if mod(P, 2) == 0
        len = Q;
    else
        len = 2*Q;
    end

    assert(mod(3*rungs*N, len) == 0, 'Trying to create a unit cell of 3 x %d (N) x %d (rungs) = %d, while the period of the mps has to be a multiple of %d due to filling %d / %d \n', N, rungs, 3*N*rungs, len, P, Q);
    disp('Code started running');

    H1 = get_Hubbard_JMpo_threeband(param, 'P', P, 'Q', Q, 'system', {'Cylinder_multiple_rungs', N, rungs}, 'convention', kwargs.convention, 'symmetries', kwargs.symmetries);

    if finalized == 4
        load(starting_name, 'gs_mps');
        w = period(gs_mps);
        assert(mod(w, N) == 0)
        mps = multiply_mps(gs_mps, rungs*N/w);
    elseif finalized == 3
        load(starting_name, 'mps');
        w = period(mps);
        assert(mod(w, N) == 0)
        mps = multiply_mps(canonicalize(mps, 'Order', 'rl'), rungs*N/w);
    elseif finalized == 2
        load(starting_name, 'gs_mps');
        mps = gs_mps;
    elseif finalized == 1
        load(starting_name, 'mps');
        mps = canonicalize(mps, 'Order', 'rl');
    else
        if strcmp(kwargs.symmetries, 'U1_U1')
            mps = get_Hubbard_mps(P, Q, 'system', {'Cylinder_multiple_rungs', 3*N, rungs});
        elseif strcmp(kwargs.symmetries, 'U1_SU2')
            mps = get_Hubbard_mps(P, Q, 'system', {'Cylinder_multiple_rungs', 3*N, rungs}, 'symmetries', 'U1_SU2');
        elseif strcmp(kwargs.symmetries, 'None_U1')
            error('Using None_U1');
            mps = get_Hubbard_mps_without_U1();
        else
            error('Invalid symmetry')
        end
    end
    disp('initialization correct');

    name = 'Hubbard_FullCylinder_threeband_N_' + string(N) + '_model_' + string(model) + '_P_' + string(P) + '_Q_' + string(Q) + '_rungs_' + string(rungs);

    [gs_mps, gs_energy, eta] = doVumps(H1, mps, vumps_way, maxiter, trunc, tol, name, 'trunc_method', kwargs.trunc_method);
end