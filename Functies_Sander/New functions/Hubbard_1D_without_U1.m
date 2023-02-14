function [gs_mps, gs_energy] = Hubbard_1D_without_U1(t, U, trunc, maxiter, tol, vumps_way, starting_name, finalized, kwargs)
    arguments
        t
        U
        mu
        trunc
        maxiter
        tol
        vumps_way
        starting_name
        finalized
        kwargs.t2 = 0;
        kwargs.V = 0;
    end
    disp('Code started running');
        
    % For atomic separation = 2å of 1D hydrogen chain: 
    % t = 1.5 en abs(U/t) = 6
    %t = 1.5;
    %U = 9;

    mu = 0;
    h_field = 0;
    N = 0;

    if mod(P, 2) == 0
        len = Q;
    else
        len = 2*Q;
    end

    H1 = get_Hubbard_JMpo(t, U, 'P', P, 'Q', Q, 'system', {'1D'}, 't2', kwargs.t2, 'V', kwargs.V, 'len', len);

    if finalized == 2
        load(starting_name, 'gs_mps');
        mps = gs_mps;
    elseif finalized == 1
        load(starting_name, 'mps');
        mps = canonicalize(mps, 'Order', 'rl');
    else
        if kwargs.U1
            mps = get_Hubbard_mps(P, Q);
        else
            mps = get_H
    end
    disp('initialization correct');

    name = 'Hubbard_1D_t_' + string(t) + '_U_' + string(U) + '_filling_' + string(P) + '_' + string(Q);
    [gs_mps, gs_energy, eta] = doVumps(H1, mps, vumps_way, maxiter, trunc, tol, name);
end