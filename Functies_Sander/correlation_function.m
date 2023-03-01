function corr_list = correlation_function(O, gs_mps, max_dist, operator_type, convention, varargin)
    % operator_type is either separate, joint, or twosite
    % separate means that 2 one-site operators are given, not connected with a
    % joint means that 2 one-site operators are given, connected with a virtual leg
    % twosite means a two-site operator is given
    % virtual leg.
    % In the 'twosite' case, O is an operator
    % In the other cases, O is a cell containing the left operator as first
    % element and the right operator as second element. They are contracted
    % over their second index.
    % If O is a cell, it consists of 2 elements, O_alpha and O_beta
    % O_alpha is put on the 1st site, O_beta on the n^th site.
    % Extra argument denotes a tolerance check.
    
    if strcmp(operator_type, 'separate')
        assert(length(O) == 2, 'Length of O should be 2')
        assert(nspaces(O{1}) == 2, 'O1 should have 2 legs')
        assert(nspaces(O{2}) == 2, 'O2 should have 2 legs')
    elseif strcmp(operator_type, 'joint')
        assert(length(O) == 2, 'Length of O should be 2')
        assert(nspaces(O{1}) == 3, 'O1 should have 3 legs')
        assert(nspaces(O{2}) == 3, 'O2 should have 3 legs')
    elseif strcmp(operator_type, 'twosite')
        assert(nspaces(O) == 4, 'O should have 4 legs')
    else
        error('Invalid operator_type');
    end

    assert(xor(strcmp(operator_type, 'twosite'),iscell(O)), 'Not a valid combination of operator type and O');

    %{
    if strcmp(convention, 'first')
        disp('');
    elseif strcmp(convention, 'conventional')
        assert(~separate, 'not implemented')
        O = tpermute(O, [4 3 1 2]);
    end
    %}
    if nargin == 6
        tol_check = true;
        tol = varargin{1};
    else
        tol_check = false;
    end

    AC = gs_mps.AC;
    AR = gs_mps.AR;

    if strcmp(operator_type, 'separate')
        O_alpha = O{1};
        O_beta = O{2};
        w = period(gs_mps);
        corr_list = zeros(1,max_dist);
        x = num2cell(zeros(1, w));

        if strcmp(convention, 'first')
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
        elseif strcmp(convention, 'conventional')
            for b = 1:w
                x{b} = contract(AC(b), [1 2 -1], conj(AC(b)), [1 3 -2], O_alpha, [3 2]);
            end
            for i = 1:max_dist
                x_final = num2cell(zeros(1, w));
                for b = 1:w
                    x_final{b} = contract(x{b}, [1 4], AR(loop(b,i,w)), [1 2 5], conj(twist(AR(loop(b,i,w)),3)), [4 3 5], O_beta, [3 2]);
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
        end
    elseif strcmp(operator_type, 'joint')
        O_alpha = O{1};
        O_beta = O{2};
        w = period(gs_mps);
        corr_list = zeros(1,max_dist);
        x = num2cell(zeros(1, w));

        if strcmp(convention, 'first')
            for b = 1:w
                x{b} = contract(AC(b), [1 2 -1], conj(AC(b)), [1 3 -2], O_alpha, [2 -3 3]);
            end
            for i = 1:max_dist
                x_final = num2cell(zeros(1, w));
                for b = 1:w
                    x_final{b} = contract(x{b}, [1 4 6], AR(loop(b,i,w)), [1 2 5], conj(twist(AR(loop(b,i,w)),3)), [4 3 5], O_beta, [2 6 3]);
                end
                new_value = mean(cell2mat(x_final));
                corr_list(i) = new_value;
                if tol_check && abs(new_value) < tol
                    corr_list = corr_list(1:i);
                    return
                end
                if i ~= max_dist
                    for b = 1:w
                        x{b} = contract(x{b}, [1 2 -3], AR(loop(b,i,w)), [1 3 -1], conj(AR(loop(b,i,w))), [2 3 -2]);
                    end
                end
            end
        elseif strcmp(convention, 'conventional')
            for b = 1:w
                x{b} = contract(AC(b), [1 2 -1], conj(AC(b)), [1 3 -2], O_alpha, [3 -3 2]);
            end
            for i = 1:max_dist
                x_final = num2cell(zeros(1, w));
                for b = 1:w
                    x_final{b} = contract(x{b}, [1 4 6], AR(loop(b,i,w)), [1 2 5], conj(twist(AR(loop(b,i,w)),3)), [4 3 5], O_beta, [3 6 2]);
                end
                new_value = mean(cell2mat(x_final));
                corr_list(i) = new_value;
                if tol_check && abs(new_value) < tol
                    corr_list = corr_list(1:i);
                    return
                end
                if i ~= max_dist
                    for b = 1:w
                        x{b} = contract(x{b}, [1 2 -3], AR(loop(b,i,w)), [1 3 -1], conj(AR(loop(b,i,w))), [2 3 -2]);
                    end
                end
            end
        end
    elseif strcmp(operator_type, 'twosite')
        w = period(gs_mps);
        corr_list = zeros(1,max_dist);
        x = num2cell(zeros(1, w));

        if strcmp(convention, 'first')
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
        elseif strcmp(convention, 'conventional')
            for b = 1:w
                x{b} = contract(AC(b), [1 2 -1], conj(AC(b)), [1 3 -2], O, [3 -4 -3 2]);
            end
            for i = 1:max_dist
                disp(i);
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
    else
        error('Not a valid operator_type')
    end
end