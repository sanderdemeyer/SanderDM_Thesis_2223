function [gs_mps, gs_energy] = Hubbard_1D_HPC(t, U, P, Q, trunc, maxiter, tol, vumps_way, starting_name, finalized)
    arguments
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
    end
    disp('Code started running');

    % For atomic separation = 2å of 1D hydrogen chain: 
    % t = 1.5 en abs(U/t) = 6
    %t = 1.5;
    %U = 9;

    h_field = 0;
    N = 0;

    if mod(P, 2) == 0
        len = Q;
    else
        len = 2*Q;
    end

    H1 = get_Hubbard_JMpo(t, U, 'P', P, 'Q', Q, 'system', {'1D'}, 't2', 0, 'V', 0, 'len', len, 'symmetries', 'U1_U1', 'mu', 0, 'convention', 'conventional');
        
    if finalized == 2
        load(starting_name, 'gs_mps');
        mps = gs_mps;
    elseif finalized == 1
        load(starting_name, 'mps');
        mps = canonicalize(mps, 'Order', 'rl');
    else
        if strcmp('U1_U1', 'U1_U1')
            mps = get_Hubbard_mps(P, Q);
        elseif strcmp('U1_U1', 'U1_SU2')
            warning('Only for half filling 1D. TBA')
            pspace = get_spaces_Hubbard_SU2(P, Q, 'half_filling', true);
            mps = UniformMps.randnc(pspace, pspace, pspace, pspace*pspace);
        elseif strcmp('U1_U1', 'None_U1')
            mps = get_Hubbard_mps_without_U1();
        else
            error('Invalid symmetry')
        end
    end
    disp('initialization correct');

    name = 'Hubbard_1D_t_' + string(t) + '_U_' + string(U) + '_filling_' + string(P) + '_' + string(Q);
    [gs_mps, gs_energy, eta] = doVumps(H1, mps, vumps_way, maxiter, trunc, tol, name);
end