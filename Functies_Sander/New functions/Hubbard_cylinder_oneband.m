function [gs_mps, gs_energy] = Hubbard_cylinder_oneband(N, model, P, Q, rungs, trunc, maxiter, tol, vumps_way, starting_name, finalized, kwargs)
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
        kwargs.mu = 0
        kwargs.convention = 'conventional'
        kwargs.trunc_method = 'TruncTotalDim'
    end

    if mod(P, 2) == 0
        len = Q;
    else
        len = 2*Q;
    end

    if strcmp(kwargs.symmetries, 'None_SU2')
        len = 2;
    end

    if strcmp(model, 'Hg_oneband1')
        t = 1;
        U = 9.48;
        t2 = -0.26;
        V = 2.364;
    elseif strcmp(model, 'La_oneband1')
        scale = 0.485;
        t = 1;
        U = 5/scale;
        t2 = -0.073/scale;
        V = 1.11/scale;
    else
        error('Model %s not implemented', model);
    end

    assert(mod(rungs*N, len) == 0, 'Trying to create a unit cell of %d (N) x %d (rungs) = %d, while the period of the mps has to be a multiple of %d due to filling %d / %d \n', N, rungs, N*rungs, len, P, Q);
    disp('Code started running');

    H1 = get_Hubbard_JMpo_oneband(t, t2, U, V, 'P', P, 'Q', Q, 'system', {'Cylinder_multiple_rungs', N, rungs}, 'convention', kwargs.convention, 'symmetries', kwargs.symmetries, 'mu', kwargs.mu);

    if iscell(finalized)
        if strcmp(finalized{1}, 'multiply')
            if finalized{2} == 2
                load(starting_name, 'gs_mps');
                w = period(gs_mps);
                assert(mod(w, N) == 0);
                mps = multiply_mps(gs_mps, rungs*N/w);
            elseif finalized{2} == 1
                load(starting_name, 'mps');
                w = period(mps);
                assert(mod(w, N) == 0);
                mps = multiply_mps(canonicalize(mps, 'Order', 'rl'), rungs*N/w);
            end
        elseif strcmp(finalized{1}, 'divide')
            if finalized{2} == 2
                load(starting_name, 'gs_mps');
                w = period(gs_mps);
                assert(w >= rungs*N);
                mps = divide_mps(gs_mps, rungs*N);
            elseif finalized{2} == 1
                load(starting_name, 'mps');
                w = period(mps);
                assert(w >= rungs*N);
                mps = divide_mps(canonicalize(mps, 'Order', 'rl'), rungs*N);
            end
        end
    elseif finalized == 2
        load(starting_name, 'gs_mps');
        mps = gs_mps;
    elseif finalized == 1
        load(starting_name, 'mps');
        mps = canonicalize(mps, 'Order', 'rl');
    else
        if strcmp(kwargs.symmetries, 'U1_U1')
            mps = get_Hubbard_mps(P, Q, 'system', {'Cylinder_multiple_rungs', N, rungs});
        elseif strcmp(kwargs.symmetries, 'U1_SU2')
            mps = get_Hubbard_mps(P, Q, 'system', {'Cylinder_multiple_rungs', N, rungs}, 'symmetries', 'U1_SU2');
        elseif strcmp(kwargs.symmetries, 'None_SU2')
            mps = get_Hubbard_mps(0, 0, 'system', {'Cylinder_multiple_rungs', N, rungs}, 'symmetries', 'None_SU2');
        elseif strcmp(kwargs.symmetries, 'None_U1')
            error('Using None_U1. This is strongly discouraged!');
            mps = get_Hubbard_mps_without_U1();
        else
            error('Invalid symmetry')
        end
    end
    disp('initialization correct');

    name = 'Hubbard_FullCylinder_oneband_' + string(model) + '_N_' + string(N) + '_P_' + string(P) + '_Q_' + string(Q) + '_rungs_' + string(rungs) + '_mu_' + string(kwargs.mu);

    [gs_mps, gs_energy, eta] = doVumps(H1, mps, vumps_way, maxiter, trunc, tol, name, 'trunc_method', kwargs.trunc_method);
end