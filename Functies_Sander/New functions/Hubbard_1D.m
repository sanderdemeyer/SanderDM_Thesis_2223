function [gs_mps, gs_energy] = Hubbard_1D(t, U, P, Q, trunc, maxiter, tol, vumps_way, starting_name, finalized)
    disp('Code started running');
        

    % For atomic separation = 2å of 1D hydrogen chain: 
    % t = 1.5 en abs(U/t) = 6
    %t = 1.5;
    %U = 9;

    mu = 0;
    h_field = 0;
    N = 0;

    H1 = get_Hubbard_JMpo(t, U, 'P', P, 'Q', Q, 'system', '1D');

    if finalized == 2
        load(starting_name, 'gs_mps');
        mps = gs_mps;
    elseif finalized == 1
        load(starting_name, 'mps');
        mps = canonicalize(mps, 'Order', 'rl');
    else
        mps = get_Hubbard_mps(P, Q);
    end
    disp('initialization correct');

    name = 'Hubbard_1D_t_' + string(t) + '_U_' + string(U) + '_filling_' + string(P) + '_' + string(Q);
    [gs_mps, gs_energy, eta] = doVumps(H1, mps, vumps_way, maxiter, trunc, tol, name);
end