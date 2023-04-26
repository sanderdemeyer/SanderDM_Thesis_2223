function [list_bb, list_br, list_rb, list_rr] = get_SC_order_parameter_april(gs_mps, P, Q, max_dist, L1, L2, L3, L4)
    %[L1, L2, L3, L4] = Hubbard_operators_superconductivity_ordered(P, Q);
%    load('operators.mat');
    disp('Operators are defined');
    AC = gs_mps.AC;
    AR = gs_mps.AR;
    w = period(gs_mps);
    b = 1;

    vspaces = zeros(1, w);
    for j = 1:w
        space = AC(j).space.dims;
        vspaces(j) = space(1);
    end
    list_bb = zeros(1, max_dist);
    list_br = zeros(1, max_dist);
    list_rb = zeros(1, max_dist);
    list_rr = zeros(1, max_dist);

    left_side = contract(AC(b), [1 2 -1], L1, [4 -3 2], conj(AC(b)), [1 4 -2]);
    option_b = contract(left_side, [1 2 3], AR(loop(b,1,w)), [1 4 -1], L2, [3 5 -3 4], conj(AR(loop(b,1,w))), [2 5 -2]);
    option_b = further(option_b, AR, loop(b,2,w));
    option_b = further(option_b, AR, loop(b,3,w));

    option_r = further(left_side, AR, loop(b,1,w));
    option_r = -contract(option_r, [1 2 3], AR(loop(b,2,w)), [1 4 -1], L2, [3 5 -3 4], conj(AR(loop(b,2,w))), [2 5 -2]);
    option_r = further(option_r, AR, loop(b,3,w));


    for i = 1:max_dist
        fprintf('Started with iteration %d of %d \n', i, max_dist);
        
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

        option_b = further(option_b, AR, loop(b,2+2*i,w));
        option_b = further(option_b, AR, loop(b,3+2*i,w));

        option_r = further(option_r, AR, loop(b,2+2*i,w));
        option_r = further(option_r, AR, loop(b,3+2*i,w));
    end
    
    function O_next = further(O, AR, i)
        O_next = contract(O, [1 2 -3], AR(i), [1 3 -1], conj(AR(i)), [2 3 -2]);
    end
end