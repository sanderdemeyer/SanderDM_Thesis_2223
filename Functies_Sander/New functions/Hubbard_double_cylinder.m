function [gs_mps, gs_energy] = Hubbard_double_cylinder(N, t, U, trunc, maxiter, tol, vumps_way, starting_name, finalized)
    % Unit cell of 2 rungs
    disp('Code started running');
    P = 1;
    Q = 1;
    mu = 0;

    H1 = get_Hubbard_JMpo(t, U, 'P', P, 'Q', Q, 'system', {'DoubleCylinder', N});

    if finalized == 2
        load(starting_name, 'gs_mps');
        mps = gs_mps;
    elseif finalized == 1
        load(starting_name, 'mps');
        mps = canonicalize(mps, 'Order', 'rl');
    else
        mps = get_Hubbard_mps(P, Q, 'system', {'DoubleCylinder', N});
    end
    disp('initialization correct');

    name = 'Hubbard_FullCylinder_half_filling_N_' + string(N) + '_t_' + string(t) + '_U_' + string(U);
    [gs_mps, gs_energy, eta] = doVumps(H1, mps, vumps_way, maxiter, trunc, tol, name);
end