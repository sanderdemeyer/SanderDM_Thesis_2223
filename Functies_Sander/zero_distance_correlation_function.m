function corr = zero_distance_correlation_function(O, gs_mps, convention)
    separate = iscell(O);

    if separate
        assert(length(O) == 2, 'Length of O should be 2')
        assert(nspaces(O{1}) == 2, 'O1 should have 2 legs')
        assert(nspaces(O{2}) == 2, 'O2 should have 2 legs')
    else
        assert(nspaces(O) == 4, 'O should have 4 legs')
    end

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