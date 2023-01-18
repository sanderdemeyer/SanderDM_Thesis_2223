function [gs_mps, gs_energy] = Hubbard_cylinder(N, t, U, trunc, maxiter, tol, vumps_way, starting_name, finalized, kwargs)
    arguments
        N
        t
        U
        trunc
        maxiter
        tol
        vumps_way
        starting_name
        finalized
        kwargs.D = 1
    end

    % Unit cell of 1 rung
    disp('Code started running');
    P = 1;
    Q = 1;
    mu = 0;

    H1 = get_Hubbard_JMpo(t, U, 'P', P, 'Q', Q, 'system', {'Cylinder', N}, 'D', kwargs.D);

    if finalized == 2
        load(starting_name, 'gs_mps');
        mps = gs_mps;
    elseif finalized == 1
        load(starting_name, 'mps');
        mps = canonicalize(mps, 'Order', 'rl');
    else
        mps = get_Hubbard_mps(P, Q, 'system', {'Cylinder', N}, 'D', kwargs.D);
    end
    disp('initialization correct');

    name = 'Hubbard_FullCylinder_half_filling_N_' + string(N) + '_t_' + string(t) + '_U_' + string(U);
    [gs_mps, gs_energy, eta] = doVumps(H1, mps, vumps_way, maxiter, trunc, tol, name);
end