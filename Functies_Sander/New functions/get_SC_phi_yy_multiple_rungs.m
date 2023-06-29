function phi_yy_list_final = get_SC_phi_yy_multiple_rungs(gs_mps, N, P, Q, max_dist, kwargs)
    arguments
        gs_mps
        N
        P
        Q
        max_dist
        kwargs.symmetries = 'U1_SU2'
        kwargs.operators = []
    end
    if isempty(kwargs.operators)
        if strcmp(kwargs.symmetries, 'U1_SU2')
            [L1, L2, L3, L4] = Hubbard_operators_superconductivity_ordered(P, Q);
        elseif strcmp(kwargs.symmetries, 'None_SU2')
            [L1, L2, L3, L4] = Hubbard_operators_superconductivity_ordered_None_SU2();
        else
            error('invalid symmetry');
        end
    else
        operators = kwargs.operators();
        L1 = operators{1};
        L2 = operators{2};
        L3 = operators{3};
        L4 = operators{4};
    end

    disp('Operators are defined');
    AC = gs_mps.AC;
    AR = gs_mps.AR;
    w = period(gs_mps);

    assert(mod(w, N) == 0);
    rungs = w/N;

    phi_yy_list_sum = zeros(1, max_dist);

    for rung = 0:rungs-1
        phi_yy_list = zeros(1, max_dist);
        L_values = cell(0, N);
    
        for b = 1:N-1
            fprintf('Started for rung = %d and distance b = %d \n', rung, b);
            L_value_b = contract(AC(mod1(rung*N+b,w)), [1 2 3], AR(mod1(rung*N+b+1,w)), [3 4 -1], L1, [5 6 2], L2, [6 9 -3 4], conj(AC(mod1(rung*N+b,w))), [1 5 7], conj(AR(mod1(rung*N+b+1,w))), [7 9 -2]);
            L_values{b} = few_steps(L_value_b, AR, mod1(rung*N+b+2,w), N-2, w);
        end
        L_N_value = contract(AC(mod1(rung*N+1,w)), [1 2 -1], L2, [-4 3 -3 2], conj(AC(mod1(rung*N+1,w))), [1 3 -2]);
        L_N_value = few_steps_4(L_N_value, AR, mod1(rung*N+2,w), N-2, w);
        L_values{N} = contract(L_N_value, [1 2 -3 4], AR(mod1(rung*N+N,w)), [1 3 -1], L1, [5 4 3], conj(AR(mod1(rung*N+N,w))), [2 5 -2]);
    
        for i = 0:max_dist-1
            fprintf('Started for distance i = %d \n', i+1);
            phi_yy_list_i = zeros(1,N);
            
            for b = 1:N-1
                phi_yy_list_i(b) = contract(L_values{b}, [1 2 3], AR(mod1((rung+i)*N+b+N,w)), [1 4 5], AR(mod1((rung+i)*N+b+N+1,w)), [5 10 9], L3, [3 6 8 4], L4, [8 11 10], conj(AR(mod1((rung+i)*N+b+N,w))), [2 6 7], twist(conj(AR(mod1((rung+i)*N+b+N+1,w))),3), [7 11 9]);
            end
            R_N_value = contract(L_values{N}, [1 2 -3], AR(mod1((rung+i)*N+N+1,w)), [1 3 -1], L4, [-4 4 3], conj(AR(mod1((rung+i)*N+N+1,w))), [2 4 -2]);
            R_N_value = few_steps_4(R_N_value, AR, mod1((rung+i)*N+N+2,w), N-2, w);
            phi_yy_list_i(N) = contract(R_N_value, [1 2 3 4], AR(mod1((rung+i)*N+2*N,w)), [1 6 5], L3, [3 7 4 6], twist(conj(AR(mod1((rung+i)*N+2*N,w))),3), [2 7 5]);
    
            phi_yy_list(i+1) = mean(phi_yy_list_i);
    
            for b = 1:N-1
                L_values{b} = few_steps(L_values{b}, AR, mod1((rung+i)*N+N+b,w), N, w);
            end
            L_values{N} = few_steps(L_values{N}, AR, mod1((rung+i)*N+N+1,w), N, w);
        end
        phi_yy_list_sum = phi_yy_list_sum + phi_yy_list;
    end

    phi_yy_list_final = phi_yy_list_sum/rungs;

    function O = few_steps(O, AR, i_start, steps, w)
        for step = 0:steps-1
            O = contract(O, [1 2 -3], AR(loop(i_start, step, w)), [1 3 -1], conj(AR(loop(i_start, step, w))), [2 3 -2]);
        end
    end
    function O = few_steps_4(O, AR, i_start, steps, w)
        for step = 0:steps-1
            O = contract(O, [1 2 -3 -4], AR(loop(i_start, step, w)), [1 3 -1], conj(AR(loop(i_start, step, w))), [2 3 -2]);
        end
    end
end