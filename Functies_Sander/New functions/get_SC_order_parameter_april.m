function [list_bb, list_br, list_rb, list_rr] = get_SC_order_parameter_april(gs_mps, N, P, Q, max_dist, L1, L2, L3, L4)
    %[L1, L2, L3, L4] = Hubbard_operators_superconductivity_ordered(P, Q);
%    load('operators.mat');
    disp('Operators are defined');
    AC = gs_mps.AC;
    AR = gs_mps.AR;
    w = period(gs_mps);
    b = 2;

    warning('fix minus signs!!');

    vspaces = zeros(1, w);
    for j = 1:w
        space = AC(j).space.dims;
        vspaces(j) = space(1);
    end
    %{
    for b = 1:w
        if mod(b,N) == 1
            disp('TBA');
        elseif mod(b,N) == 0
            disp('TBA');
        else
    %}
    option_l = -contract(AC(loop(b,-N+1,w)), [1 2 -1], L2, [-3 3 -4 2], conj(AC(loop(b,-N+1,w))), [1 3 -2]);
    option_l = few_steps_4(option_l, AR, loop(b,-N+2,w), N-2, w);
    option_l = contract(option_l, [1 5 3 -3], AC(b), [1 2 -1], conj(AC(b)), [5 4 -2], L1, [4 3 2]);
    option_l = few_steps(option_l, AR, loop(b,1,w), 2*N-1, w);

    option_t = contract(AC(loop(b,-1,w)), [1 2 -1], L2, [-3 3 -4 2], conj(AC(loop(b,-1,w))), [1 3 -2]);
    option_t = contract(option_t, [1 5 3 -3], AC(b), [1 2 -1], conj(AC(b)), [5 4 -2], L1, [4 3 2]);
    option_t = few_steps(option_t, AR, loop(b,1,w), 2*N-1, w);

    base_br = contract(AC(b), [1 2 -1], L1, [4 -3 2], conj(AC(b)), [1 4 -2]);

    option_b = contract(base_br, [1 2 3], AR(loop(b,1,w)), [1 4 -1], L2, [3 5 -3 4], conj(AR(loop(b,1,w))), [2 5 -2]);
    option_b = few_steps(option_b, AR, loop(b,2,w), 2*N-2, w);

    option_r = few_steps(base_br, AR, loop(b,1,w), N-1, w);
    option_r = -contract(option_r, [1 2 3], AR(loop(b,N,w)), [1 4 -1], L2, [3 5 -3 4], conj(AR(loop(b,1,w))), [2 5 -2]);
    option_r = few_steps(option_r, AR, loop(b,N+1,w), N-1, w);

    left_side = option_l + option_t + option_b + option_r;


    %{
    left_side = contract(AC(b), [1 2 -1], L1, [4 -3 2], conj(AC(b)), [1 4 -2]);
    option_b = contract(left_side, [1 2 3], AR(loop(b,1,w)), [1 4 -1], L2, [3 5 -3 4], conj(AR(loop(b,1,w))), [2 5 -2]);
    option_b = further(option_b, AR, loop(b,2,w));
    option_b = further(option_b, AR, loop(b,3,w));

    option_r = further(left_side, AR, loop(b,1,w));
    option_r = -contract(option_r, [1 2 3], AR(loop(b,2,w)), [1 4 -1], L2, [3 5 -3 4], conj(AR(loop(b,2,w))), [2 5 -2]);
    option_r = further(option_r, AR, loop(b,3,w));
    %}

    for i = 1:max_dist
        fprintf('Started with iteration %d of %d \n', i, max_dist);
        
        option_l = -contract(left_side, [1  2 -4], AC(loop(b,N+N*i,w)), [1 3 -1], L4, [-3 4 3], conj(AC(loop(b,N+N*i,w))), [2 4 -2]);
        option_l = few_steps_4(option_l, AR, loop(b, N+1+N*i, w), N-1, w);
        option_l = contract(option_l, [1 6 4 3], AC(loop(b,2*N+N*i,w)), [1 2 7], L3, [3 5 4 2], conj(AC()), [6 5 7]);

        option_t = few_steps(left_side, AR, loop(b, N + N*i, w), N-2, w);
        option_t = contract(option_t, [1  2 -3], AC(2*N-2+N*i), [1 3 -1], L4, [4 -4 3], conj(AC(N+N*i)), [2 4 -2]);

        base = few_steps(left_side, AR, loop(b,N+N*i), N, w);
        base = contract(base, [1 5 3], AC(loop(b,2*N+N*i,w)), [1 2 -1], L3, [3 4 -3 2], conj(AC(loop(b,2*N+N*i,w))), [5 4 -2]);

        option_b = contract(base, [1 5 3], AC(loop(b,2*N+1+N*i,w)), [1 2 6], L4, [3 4 2], conj(AC(loop(b,2*N+1+N*i,w))), [5 4 6]);

        option_r = few_steps(base, AR, loop(b,2*N+1+N*i,w), N-2, w);
        option_r = -contract(option_r, [1 5 3], AC(loop(b,3*N+N*i,w)), [1 2 6], L4, [3 4 2], conj(AC(loop(b,3*N+N*i,w))), [5 4 6]);

        list_tot(i) = option_l + option_t + option_b + option_r;

        %{
        finalizing_b = contract(option_b, [1 2 3], AR(loop(b,2+2*i,w)), [1 4 -1], L3, [3 5 -3 4], conj(AR(loop(b,2+2*i,w))), [2 5 -2]);
        finalizing_r = contract(option_r, [1 2 3], AR(loop(b,2+2*i,w)), [1 4 -1], L3, [3 5 -3 4], conj(AR(loop(b,2+2*i,w))), [2 5 -2]);
        
        el_bb = contract(finalizing_b, [1 2 3], AR(loop(b,3+2*i,w)), [1 4 6], L4, [3 5 4], conj(AR(loop(b,3+2*i,w))), [2 5 6]);
        el_rb = contract(finalizing_r, [1 2 3], AR(loop(b,3+2*i,w)), [1 4 6], L4, [3 5 4], conj(AR(loop(b,3+2*i,w))), [2 5 6]);

        el_b = further(finalizing_b, AR, loop(b,3+2*i,w));
        el_r = further(finalizing_r, AR, loop(b,3+2*i,w));

        el_br = contract(el_b, [1 2 3], AR(loop(b,4+2*i,w)), [1 4 6], L4, [3 5 4], conj(AR(loop(b,4+2*i,w))), [2 5 6]);
        el_rr = -contract(el_r, [1 2 3], AR(loop(b,4+2*i,w)), [1 4 6], L4, [3 5 4], conj(AR(loop(b,4+2*i,w))), [2 5 6]);

        list_bb(i) = el_bb;
        list_rb(i) = el_rb;
        list_br(i) = el_br;
        list_rr(i) = el_rr;
        %}
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