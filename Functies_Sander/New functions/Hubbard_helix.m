function [gs_mps, gs_energy] = Hubbard_helix(N, t, U, P, Q, rungs, trunc, maxiter, tol, vumps_way, starting_name, finalized, kwargs)
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
        kwargs.len = []
        kwargs.t2 = 0;
        kwargs.V = 0;
        kwargs.symmetries = 'U1_U1'
        kwargs.convention = 'conventional'
    end

    if isempty(kwargs.len)
        if mod(P, 2) == 0
            len = Q;
        else
            len = 2*Q;
        end
    else
        len = kwargs.len;
    end

    assert(mod(rungs*N, len) == 0, 'Trying to create a unit cell of %d (N) x %d (rungs) = %d, while the period of the mps has to be a multiple of %d \n', N, rungs, N*rungs, len);
    disp('Code started running');

    H1 = get_Hubbard_JMpo(t, U, 'P', P, 'Q', Q, 'system', {'Helix_multiple_rungs', N, rungs}, 't2', kwargs.t2, 'V', kwargs.V, 'symmetries', kwargs.symmetries, 'convention', kwargs.convention);

    if finalized == 2
        load(starting_name, 'gs_mps');
        mps = gs_mps;
    elseif finalized == 1
        load(starting_name, 'mps');
        mps = canonicalize(mps, 'Order', 'rl');
    else
        if strcmp(kwargs.symmetries, 'U1_U1')
            mps = get_Hubbard_mps(P, Q, 'system', {'Cylinder_multiple_rungs', N, rungs});
        elseif strcmp(kwargs.symmetries, 'U1_SU2')
            %pspace = get_spaces_Hubbard_SU2(P, Q);
%            mps = UniformMps.randnc(pspace, pspace, pspace, pspace*pspace);
            mps = get_Hubbard_mps(P, Q, 'system', {'Cylinder_multiple_rungs', N, rungs}, 'symmetries', 'U1_SU2');
        elseif strcmp(kwargs.symmetries, 'None_U1')
            error('Using None_U1')
            mps = get_Hubbard_mps_without_U1();
        else
            error('Invalid symmetry')
        end
    end
    disp('initialization correct');
    name = 'Hubbard_Helix_N_' + string(N) + '_t_' + string(t) + '_U_' + string(U) + '_P_' + string(P) + '_Q_' + string(Q) + '_rungs_' + string(rungs) + '_t2_' + string(kwargs.t2) + '_V_' + string(kwargs.V);
    [gs_mps, gs_energy, eta] = doVumps(H1, mps, vumps_way, maxiter, trunc, tol, name);
end