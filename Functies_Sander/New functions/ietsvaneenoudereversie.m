function [gs_mps, gs_energy] = Hubbard_cylinder(N, t, U, P, Q, rungs, trunc, maxiter, tol, vumps_way, starting_name, finalized, kwargs)
    arguments
        N
        t
        U
        P
        Q
        rungs
        trunc
        maxiter
        tol
        vumps_way
        starting_name
        finalized
        kwargs.D = 1
        kwargs.t2 = 0
        kwargs.V = 0
        kwargs.symmetries = 'U1_U1' % Symmetry of the charge and spin sector. Fermionic parity is always implemented.
        kwargs.mu = 0
        kwargs.convention = 'conventional'
    end

    if mod(P, 2) == 0
        len = Q;
    else
        len = 2*Q;
    end

    assert(mod(rungs*N, len) == 0, 'Trying to create a unit cell of %d (N) x %d (rungs) = %d, while the period of the mps has to be a multiple of %d \n', N, rungs, N*rungs, len);
    disp('Code started running');
    mu = 0;

    H1 = get_Hubbard_JMpo(t, U, 'P', P, 'Q', Q, 'system', {'Cylinder_multiple_rungs', N, rungs}, 'D', kwargs.D, 'convention', kwargs.convention, 'symmetries', kwargs.symmetries);

    if finalized == 2
        load(starting_name, 'gs_mps');
        mps = gs_mps;
    elseif finalized == 1
        load(starting_name, 'mps');
        mps = canonicalize(mps, 'Order', 'rl');
    else
        if strcmp(kwargs.symmetries, 'U1_U1')
            mps = get_Hubbard_mps(P, Q, 'system', {'Cylinder_multiple_rungs', N, rungs}, 'D', kwargs.D);
        elseif strcmp(kwargs.symmetries, 'U1_SU2')
            mps = get_Hubbard_mps(P, Q, 'system', {'Cylinder_multiple_rungs', N, rungs}, 'D', kwargs.D, 'symmetries', 'U1_SU2');
        elseif strcmp(kwargs.symmetries, 'None_U1')
            error('go');
            mps = get_Hubbard_mps_without_U1();
        else
            error('Invalid symmetry')
        end
    end
    disp('initialization correct');

    name = 'Hubbard_cylinder_N_' + string(N) + '_t_' + string(t) + '_U_' + string(U) + '_P_' + string(P) + '_Q_' + string(Q) + '_rungs_' + string(rungs) + '_trunc_' + string(trunc) + '_t2_' + string(kwargs.t2) + '_V_' + string(kwargs.V);
    [gs_mps, gs_energy, eta] = doVumps(H1, mps, vumps_way, maxiter, trunc, tol, name);
end