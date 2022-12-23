function corr_list = correlation_function(O, gs_mps, max_dist, varargin)
% Takes an 2-site operator O_{ij}, together with the properties of an
% optimized MPS (AC1, AC2, AR1, AR2), which is thus a 2-site MPS
% It returns a list, where the expectation value of H_{i,i+dist} is
% calculated, where dist is the index of the list. The list runs until
% max_dist.
% 2 casess:
% separate = false. The correlation function of a (2,2) mpo is calculated
% separate = true. The correlation function of 2 (1,1) mpo's is calculated.
% The mpo's should be given in a cell.
    separate = iscell(O);
    if separate
        assert(length(O) == 2, 'Length of O should be 2')
        assert(nspaces(O{1}) == 2, 'O1 should have 2 legs')
        assert(nspaces(O{2}) == 2, 'O2 should have 2 legs')
    else
        assert(nspaces(O) == 4, 'O should have 4 legs')
    end

    if nargin == 4
        tol_check = true;
        tol = varargin{1};
    else
        tol_check = false;
    end

    AC = gs_mps.AC;
    AR = gs_mps.AR;
    %{
    % Code from before 11/12/2022, code is more general below.
    if number_site == 1
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
    %}
    if separate
        O_alpha = O{1};
        O_beta = O{2};
        w = period(gs_mps);
        corr_list = zeros(1,max_dist);
        x = num2cell(zeros(1, w));
        for b = 1:w
            x{b} = contract(AC(b), [1 2 -1], conj(AC(b)), [1 3 -2], O_alpha, [2 3]);
        end
        for i = 1:max_dist
            x_final = num2cell(zeros(1, w));
            for b = 1:w
                x_final{b} = contract(x{b}, [1 4], AR(loop(b,i,w)), [1 2 5], conj(twist(AR(loop(b,i,w)),3)), [4 3 5], O_beta, [2 3]);
            end
            new_value = mean(cell2mat(x_final));
            corr_list(i) = new_value;
            if tol_check && abs(new_value) < tol
                corr_list = corr_list(1:i);
                return
            end
            if i ~= max_dist
                for b = 1:w
                    x{b} = contract(x{b}, [1 2], AR(loop(b,i,w)), [1 3 -1], conj(AR(loop(b,i,w))), [2 3 -2]);
                end
            end
        end
    else
        w = period(gs_mps);
        corr_list = zeros(1,max_dist);
        x = num2cell(zeros(1, w));
        for b = 1:w
            x{b} = contract(AC(b), [1 2 -1], conj(AC(b)), [1 3 -2], O, [2 -3 3 -4]);
        end
        for i = 1:max_dist
            x_final = num2cell(zeros(1, w));
            for b = 1:w
                x_final{b} = contract(x{b}, [1 3 2 4], AR(loop(b,i,w)), [1 2 5], conj(twist(AR(loop(b,i,w)),3)), [3 4 5]);
            end
            new_value = mean(cell2mat(x_final));
            corr_list(i) = new_value;
            if tol_check && abs(new_value) < tol
                corr_list = corr_list(1:i);
                return
            end
            if i ~= max_dist
                for b = 1:w
                    x{b} = contract(x{b}, [1 2 -3 -4], AR(loop(b,i,w)), [1 3 -1], conj(AR(loop(b,i,w))), [2 3 -2]);
                end
            end
        end
    end
end