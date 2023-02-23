function [gs_mps, gs_energy] = Hubbard_helix(N, t, U, P, Q, trunc, maxiter, tol, vumps_way, starting_name, finalized, kwargs)
    arguments
        N
        t
        U
        P
        Q
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

    disp('Code started running');
    if isempty(kwargs.len)
        if mod(P, 2) == 0
            len = Q;
        else
            len = 2*Q;
        end
    else
        len = kwargs.len;
    end
    H1 = get_Hubbard_JMpo(t, U, 'P', P, 'Q', Q, 'system', {'Helix', N}, 't2', kwargs.t2, 'V', kwargs.V, 'len', len, 'symmetries', kwargs.symmetries, 'convention', kwargs.convention);

    if finalized == 2
        load(starting_name, 'gs_mps');
        mps = gs_mps;
    elseif finalized == 1
        load(starting_name, 'mps');
        mps = canonicalize(mps, 'Order', 'rl');
    else
        if strcmp(kwargs.symmetries, 'U1_U1')
            mps = get_Hubbard_mps(P, Q, 'system', {'Helix', N}, 'len', len);
        elseif strcmp(kwargs.symmetries, 'U1_SU2')
            %pspace = get_spaces_Hubbard_SU2(P, Q);
%            mps = UniformMps.randnc(pspace, pspace, pspace, pspace*pspace);
            mps = get_Hubbard_mps(P, Q, 'system', {'Helix', N}, 'symmetries', 'U1_SU2', 'len', len);
        elseif strcmp(kwargs.symmetries, 'None_U1')
            error('Using None_U1')
            mps = get_Hubbard_mps_without_U1();
        else
            error('Invalid symmetry')
        end
    end
    disp('initialization correct');
    name = 'Hubbard_Helix_N_' + string(N) + '_t_' + string(t) + '_U_' + string(U) + '_P_' + string(P) + '_Q_' + string(Q) + '_t2_' + string(kwargs.t2) + '_V_' + string(kwargs.V);
    [gs_mps, gs_energy, eta] = doVumps(H1, mps, vumps_way, maxiter, trunc, tol, name);
end