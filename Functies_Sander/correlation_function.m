function corr_list = correlation_function(O, gs_mps, number_site, fermion, max_dist, varargin)
% Takes an 2-site operator O_{ij}, together with the properties of an
% optimized MPS (AC1, AC2, AR1, AR2), which is thus a 2-site MPS
% It returns a list, where the expectation value of H_{i,i+dist} is
% calculated, where dist is the index of the list. The list runs until
% max_dist.
    if number_site == 1
        AC = gs_mps.AC;
        AR = gs_mps.AR;
        corr_list = zeros(1,max_dist);
        x = contract(AC, [1 2 -1], conj(twist(AC,3)), [1 3 -2], O, [2 -3 3 -4]);
        for i = 1:max_dist
            if fermion
                x_final = contract(x, [1 3 2 4], AR, [1 2 5], conj(twist(AR,3)), [3 4 5]);
            else
                x_final = contract(x, [1 3 2 4], AR, [1 2 5], conj(AR), [3 4 5]);
            end
            corr_list(i) = x_final;
            if i ~= max_dist
                x = contract(x, [1 2 -3 -4], AR, [1 3 -1], conj(AR), [2 3 -2]);
            end
        end
    elseif number_site == 2
        AC1 = gs_mps.AC(1);
        AC2 = gs_mps.AC(2);
        AR1 = gs_mps.AR(1);
        AR2 = gs_mps.AR(2);

        corr_list = zeros(1,max_dist);
        x1 = contract(AC1, [1 2 -1], conj(AC1), [1 3 -2], O, [2 -3 3 -4]);
        x2 = contract(AC2, [1 2 -1], conj(AC2), [1 3 -2], O, [2 -3 3 -4]);
        for i = 1:max_dist
            if mod(i, 2) == 0
                if fermion
                    x_final1 = contract(x1, [1 3 2 4], AR1, [1 2 5], conj(twist(AR1,3)), [3 4 5]);
                    x_final2 = contract(x2, [1 3 2 4], AR2, [1 2 5], conj(twist(AR2,3)), [3 4 5]);
                else
                    x_final1 = contract(x1, [1 3 2 4], AR1, [1 2 5], conj(AR1), [3 4 5]);
                    x_final2 = contract(x2, [1 3 2 4], AR2, [1 2 5], conj(AR2), [3 4 5]);
                end
                corr_list(i) = (x_final1 + x_final2)/2;
                if i ~= max_dist
                    x1 = contract(x1, [1 2 -3 -4], AR1, [1 3 -1], conj(AR1), [2 3 -2]);
                    x2 = contract(x2, [1 2 -3 -4], AR2, [1 3 -1], conj(AR2), [2 3 -2]);
                end
            else
                if fermion
                    x_final1 = contract(x1, [1 3 2 4], AR2, [1 2 5], conj(twist(AR2,3)), [3 4 5]);
                    x_final2 = contract(x2, [1 3 2 4], AR1, [1 2 5], conj(twist(AR1,3)), [3 4 5]);
                else
                    x_final1 = contract(x1, [1 3 2 4], AR2, [1 2 5], conj(AR2), [3 4 5]);
                    x_final2 = contract(x2, [1 3 2 4], AR1, [1 2 5], conj(AR1), [3 4 5]);
                end
                corr_list(i) = (x_final1 + x_final2)/2;
                if i ~= max_dist
                    x1 = contract(x1, [1 2 -3 -4], AR2, [1 3 -1], conj(AR2), [2 3 -2]);
                    x2 = contract(x2, [1 2 -3 -4], AR1, [1 3 -1], conj(AR1), [2 3 -2]);
                end
            end    
        end
    end
end