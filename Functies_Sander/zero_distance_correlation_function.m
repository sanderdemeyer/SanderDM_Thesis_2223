function corr = zero_distance_correlation_function(O, gs_mps, operator_type, convention, kwargs)
    arguments
        O
        gs_mps
        operator_type
        convention
        kwargs.swap = false
    end
    % Implements the zero-distance value of the correlation function
    % If O is a cell, it consists of 2 elements, O_alpha and O_beta
    % O_alpha is put on the 1st site, O_beta on the n^th site in the
    % general correlation function.
    % Is a is not a cell, it consists of a single operator defined as 
    % O = O_alpha O_beta. This causes O_beta to be put on top of O_alpha
    % If kwargs.swap = true, O_alpha is put on top of O_beta. In this
    % function, this is done by simply swapping the definitions of O_alpha
    % and O_beta.

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

    AC = gs_mps.AC;

    if strcmp(operator_type, 'separate')
        if kwargs.swap
            O_alpha = O{2};
            O_beta = O{1};
        else
            O_alpha = O{1};
            O_beta = O{2};
        end
        w = period(gs_mps);
        x = num2cell(zeros(1, w));
        if strcmp(convention, 'first')
            for b = 1:w
                x{b} = contract(AC(b), [1 2 5], O_beta, [2 3], O_alpha, [3 4], twist(conj(AC(b)), 3), [1 4 5]);
            end
        elseif strcmp(convention, 'conventional')
            for b = 1:w
                x{b} = contract(AC(b), [1 2 5], O_beta, [3 2], O_alpha, [4 3], twist(conj(AC(b)), 3), [1 4 5]);
            end
        end
    elseif strcmp(operator_type, 'joint')
        if kwargs.swap
            O_alpha = O{2};
            O_beta = O{1};
        else
            O_alpha = O{1};
            O_beta = O{2};
        end
        w = period(gs_mps);
        x = num2cell(zeros(1, w));
        if strcmp(convention, 'first')
            for b = 1:w
                x{b} = contract(AC(b), [1 2 6], O_beta, [2 3 4], O_alpha, [4 3 5], twist(conj(AC(b)), 3), [1 5 6]);
            end
        elseif strcmp(convention, 'conventional')
            for b = 1:w
                x{b} = contract(AC(b), [1 2 6], O_beta, [4 3 2], O_alpha, [5 3 4], twist(conj(AC(b)), 3), [1 5 6]);
            end
        end
    elseif strcmp(operator_type, 'twosite')
        w = period(gs_mps);
        x = num2cell(zeros(1, w));
        if strcmp(convention, 'first')
            if kwargs.swap
                for b = 1:w
                    x{b} = contract(AC(b), [1 2 5], O, [2 3 3 4], twist(conj(AC(b)), 3), [1 4 5]);
                end
            else
                for b = 1:w
                    x{b} = contract(AC(b), [1 2 5], O, [3 2 4 3], twist(conj(AC(b)), 3), [1 4 5]);
                end
            end
        elseif strcmp(convention, 'conventional')
                % The extra twist in the contraction is put because contract 
                % automatically takes the supertrace, while the trace
                % is wanted.
            if kwargs.swap
                for b = 1:w
                    x{b} = contract(AC(b), [1 2 5], twist(O,1), [3 4 3 2], twist(conj(AC(b)), 3), [1 4 5]);
                end
            else
                for b = 1:w
                    x{b} = contract(AC(b), [1 2 5], twist(O,2), [4 3 2 3], twist(conj(AC(b)), 3), [1 4 5]);
                end
            end
        end
    end
    corr = mean(cell2mat(x));

    %{
    if strcmp(convention, 'first')
        disp('');
    elseif strcmp(convention, 'conventional')
        assert(~separate, 'not implemented')
        O = tpermute(O, [4 3 1 2]);
    end
    AC = gs_mps.AC;

    if separate
        O_alpha = O{1};
        O_beta = O{2};
        w = period(gs_mps);
        x = num2cell(zeros(1, w));
        for b = 1:w
            x{b} = contract(AC(b), [1 2 5], O_alpha, [2 3], O_beta, [3 4], twist(conj(AC(b)), 3), [1 4 5]);
        end
        corr = mean(cell2mat(x));
    else
         corr = mean(cell2mat(x));
    end
    %}
end
