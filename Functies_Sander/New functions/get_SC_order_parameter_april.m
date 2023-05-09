function [all_lists, lists_average] = get_SC_order_parameter_april(gs_mps, N, P, Q, max_dist, file, bonddim, kwargs)
    arguments
        gs_mps
        N
        P
        Q
        max_dist
        file
        bonddim
        kwargs.symmetries = 'U1_SU2'
    end
    
    if strcmp(kwargs.symmetries, 'U1_SU2')
        [L1, L2, L3, L4] = Hubbard_operators_superconductivity_ordered(P, Q);
    elseif strcmp(kwargs.symmetries, 'None_SU2')
        [L1, L2, L3, L4] = Hubbard_operators_superconductivity_ordered_None_SU2();
    else
        error('invalid symmetry');
    end
%    load('operators.mat');
    disp('Operators are defined');
    AC = gs_mps.AC;
    AR = gs_mps.AR;
    w = period(gs_mps);

    warning('fix minus signs!!');

    vspaces = zeros(1, w);
    for j = 1:w
        space = AC(j).space.dims;
        vspaces(j) = space(1);
    end

    all_lists = cell(1, w);
    for b = 1:w

        option_l = -contract(AC(loop(b,-N,w)), [1 2 -1], L2, [-3 3 -4 2], conj(AC(loop(b,-N,w))), [1 3 -2]);
        option_l = few_steps_4(option_l, AR, loop(b,-N+1,w), N-1, w);
        option_l = contract(option_l, [1 5 3 -3], AC(loop(b,0,w)), [1 2 -1], conj(AC(loop(b,0,w))), [5 4 -2], L1, [4 3 2]);
        option_l = few_steps(option_l, AR, loop(b,1,w), 2*N-1, w);
    
        base_br = contract(AC(loop(b,0,w)), [1 2 -1], L1, [4 -3 2], conj(AC(loop(b,0,w))), [1 4 -2]);
    
        if mod(b,N) == 1
            option_t = few_steps(base_br, AR, loop(b,1,w), N-2, w);
            option_t = contract(option_t, [1 2 3], AR(loop(b,N-1,w)), [1 4 -1], L2, [3 5 -3 4], conj(AR(loop(b,N-1,w))), [2 5 -2]);
            option_t = few_steps(option_t, AR, loop(b,N,w), N, w);
        else
            option_t = contract(AC(loop(b,-1,w)), [1 2 -1], L2, [-3 3 -4 2], conj(AC(loop(b,-1,w))), [1 3 -2]);
            option_t = contract(option_t, [1 5 3 -3], AC(loop(b,0,w)), [1 2 -1], conj(AC(loop(b,0,w))), [5 4 -2], L1, [4 3 2]);
            option_t = few_steps(option_t, AR, loop(b,1,w), 2*N-1, w);
        end
    
    
        if mod(b,N) == 0
            option_b = contract(AC(loop(b,-N+1,w)), [1 2 -1], L2, [-3 3 -4 2], conj(AC(loop(b,-N+1,w))), [1 3 -2]);
            option_b = few_steps_4(option_b, AR, loop(b,-N+2,w), N-2, w);
            option_b = contract(option_b, [1 5 3 -3], AC(loop(b,0,w)), [1 2 -1], conj(AC(loop(b,0,w))), [5 4 -2], L1, [4 3 2]);        
            option_b = few_steps(option_b, AR, loop(b,1,w), 2*N-1, w);
        else
            option_b = contract(base_br, [1 2 3], AR(loop(b,1,w)), [1 4 -1], L2, [3 5 -3 4], conj(AR(loop(b,1,w))), [2 5 -2]);
            option_b = few_steps(option_b, AR, loop(b,2,w), 2*N-2, w);
        end
        
        option_r = few_steps(base_br, AR, loop(b,1,w), N-1, w);
        option_r = -contract(option_r, [1 2 3], AR(loop(b,N,w)), [1 4 -1], L2, [3 5 -3 4], conj(AR(loop(b,N,w))), [2 5 -2]);
        option_r = few_steps(option_r, AR, loop(b,N+1,w), N-1, w);
    
        left_side = option_l + option_t + option_b + option_r;
    
        for i = 1:max_dist
            fprintf('Started with iteration %d of %d for b = %d. For file %d, bonddim %d \n', i, max_dist, b, file, bonddim);
            
            option_l = -contract(left_side, [1  2 -4], AC(loop(b,N+N*i,w)), [1 3 -1], L4, [-3 4 3], conj(AC(loop(b,N+N*i,w))), [2 4 -2]);
            option_l = few_steps_4(option_l, AR, loop(b, N+1+N*i, w), N-1, w);
            option_l = contract(option_l, [1 6 4 3], AC(loop(b,2*N+N*i,w)), [1 2 7], L3, [3 5 4 2], twist(conj(AC(loop(b,2*N+N*i,w))),3), [6 5 7]);
    
            base = few_steps(left_side, AR, loop(b,N+N*i,w), N, w);
            base = contract(base, [1 5 3], AC(loop(b,2*N+N*i,w)), [1 2 -1], L3, [3 4 -3 2], conj(AC(loop(b,2*N+N*i,w))), [5 4 -2]);
    
            if mod(b,N) == 1
                option_t = few_steps(base, AR, loop(b, 2*N+1 + N*i, w), N-2, w);
                option_t = contract(option_t, [1 5 3], AC(loop(b,3*N-1+N*i,w)), [1 2 6], L4, [3 4 2], twist(conj(AC(loop(b,3*N-1+N*i,w))),3), [5 4 6]);
            else
                option_t = few_steps(left_side, AR, loop(b, N + N*i, w), N-2, w);
                option_t = contract(option_t, [1 2 -4], AC(loop(b,2*N-2+N*i,w)), [1 3 -1], L4, [-3 4 3], conj(AC(loop(b,2*N-2+N*i,w))), [2 4 -2]);
                option_t = contract(option_t, [1 5 4 3], AC(loop(b,2*N-1+N*i,w)), [1 2 7], L3, [3 6 4 2], twist(conj(AC(loop(b,2*N-1+N*i,w))),3), [5 6 7]);
            end
    
    
            if mod(b,N) == 0
                option_b = contract(left_side, [1  2 -4], AC(loop(b,N+N*i,w)), [1 3 -1], L4, [-3 4 3], conj(AC(loop(b,N+N*i,w))), [2 4 -2]);
                option_b = few_steps_4(option_b, AR, loop(b,N+1+N*i,w), N-2, w);
                option_b = contract(option_b, [1 6 4 3], AC(loop(b,2*N-1+N*i,w)), [1 2 7], L3, [3 5 4 2], twist(conj(AC(loop(b,2*N-1+N*i,w))),3), [6 5 7]);
            else
                option_b = contract(base, [1 5 3], AC(loop(b,2*N+1+N*i,w)), [1 2 6], L4, [3 4 2], twist(conj(AC(loop(b,2*N+1+N*i,w))),3), [5 4 6]);
            end
    
            option_r = few_steps(base, AR, loop(b,2*N+1+N*i,w), N-2, w);
            option_r = -contract(option_r, [1 5 3], AC(loop(b,3*N-1+N*i,w)), [1 2 6], L4, [3 4 2], twist(conj(AC(loop(b,3*N-1+N*i,w))),3), [5 4 6]);
    
            list_tot(i) = option_l + option_t + option_b + option_r;
    
            left_side = few_steps(left_side, AR, loop(b,N+N*i,w), N, w);
        end

        all_lists{b} = list_tot;
    end

    lists_average = zeros(1, max_dist);
    for i = 1:max_dist
       som = 0;
       for b = 1:w
           som = som + all_lists{b}(i);
       end
       lists_average(i) = som;
    end
    
    function O_next = further(O, AR, i)
        O_next = contract(O, [1 2 -3], AR(i), [1 3 -1], conj(AR(i)), [2 3 -2]);
    end

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