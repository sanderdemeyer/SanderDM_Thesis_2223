function [gs_mps, gs_energy] = Hubbard_helix(N, t, U, trunc, maxiter, tol, vumps_way, starting_name, finalized)
    disp('Code started running');
    P = 1;
    Q = 1;
    mu = 0;

    H1 = get_Hubbard_JMpo(t, U, 'P', P, 'Q', Q, 'system', {'Helix', N});

    if finalized == 2
        load(starting_name, 'gs_mps');
        mps = gs_mps;
    elseif finalized == 1
        load(starting_name, 'mps');
        mps = canonicalize(mps, 'Order', 'rl');
    else
        mps = get_Hubbard_mps(P, Q, 'system', {'Helix', N});
    end
    disp('initialization correct');
    
    name = 'Hubbard_helix_N_' + string(N) + '_t_' + string(t) + '_U_' + string(U) + '_filling_' + string(P) + '_' + string(Q);
    [gs_mps, gs_energy, eta] = doVumps(H1, mps, vumps_way, maxiter, trunc, tol, name);
end