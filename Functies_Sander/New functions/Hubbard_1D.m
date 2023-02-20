function [gs_mps, gs_energy] = Hubbard_1D(t, U, P, Q, trunc, maxiter, tol, vumps_way, starting_name, finalized, kwargs)
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
        kwargs.t2 = 0
        kwargs.V = 0
        kwargs.symmetries = 'U1_U1' % Symmetry of the charge and spin sector. Fermionic parity is always implemented.
        kwargs.mu = 0
        kwargs.convention = 'conventional'
    end
    disp('Code started running');

    % For atomic separation = 2å of 1D hydrogen chain: 
    % t = 1.5 en abs(U/t) = 6
    %t = 1.5;
    %U = 9;

    assert(~strcmp(kwargs.symmetries, 'None_U1'), 'Wrong results for symmetries None x U1.')
    h_field = 0;
    N = 0;

    if mod(P, 2) == 0
        len = Q;
    else
        len = 2*Q;
    end

    H1 = get_Hubbard_JMpo(t, U, 'P', P, 'Q', Q, 'system', {'1D'}, 't2', kwargs.t2, 'V', kwargs.V, 'len', len, 'symmetries', kwargs.symmetries, 'mu', kwargs.mu, 'convention', kwargs.convention);
        
    if finalized == 2
        load(starting_name, 'gs_mps');
        mps = gs_mps;
    elseif finalized == 1
        load(starting_name, 'mps');
        mps = canonicalize(mps, 'Order', 'rl');
    else
        if strcmp(kwargs.symmetries, 'U1_U1')
            mps = get_Hubbard_mps(P, Q);
        elseif strcmp(kwargs.symmetries, 'U1_SU2')
            %pspace = get_spaces_Hubbard_SU2(P, Q);
%            mps = UniformMps.randnc(pspace, pspace, pspace, pspace*pspace);
            mps = get_Hubbard_mps(P, Q, 'symmetries', 'U1_SU2');
        elseif strcmp(kwargs.symmetries, 'None_U1')
            mps = get_Hubbard_mps_without_U1();
        else
            error('Invalid symmetry')
        end
    end
    disp('initialization correct');

    name = 'Hubbard_1D_t_' + string(t) + '_U_' + string(U) + '_filling_' + string(P) + '_' + string(Q);
    [gs_mps, gs_energy, eta] = doVumps(H1, mps, vumps_way, maxiter, trunc, tol, name);
end