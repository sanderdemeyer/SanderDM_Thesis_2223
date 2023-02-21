function [charge_occupancies, spin_occupancies] = get_occupancies(gs_mps, P, Q, SU2, plot, N, rungs)
    % Using convention = conventional
    spin_occupancies = 0;
    if SU2
        [pspace, ~, trivspace] = get_spaces_Hubbard_SU2(P, Q);
        onesite = Tensor(pspace, pspace);
        number_tot = fill_tensor(onesite, num2cell([0 1 2]));
    else
        [pspace, ~, trivspace] = get_spaces_Hubbard_asymmetric(P, Q);
        onesite = Tensor(pspace, pspace);
        number_tot = fill_tensor(onesite, num2cell([0 1 1 2]));
        spin_tot = fill_tensor(onesite, num2cell([0 -1 1 0]));
    end

    AC = gs_mps.AC;
    w = period(gs_mps);
    charge_occupancies = zeros(1, w);
    for b = 1:w
        charge_occupancies(b) = contract(AC(b), [1 2 3], number_tot, [4 2], twist(conj(AC(b)),3), [1 4 3]);
    end
    
    if ~SU2
        spin_occupancies = zeros(1, w);
        for b = 1:w
            spin_occupancies(b) = contract(AC(b), [1 2 3], spin_tot, [4 2], twist(conj(AC(b)),3), [1 4 3]);
        end
    end

    if plot
        X = arrayfun(@(x) repmat(x, 1, N), 1:rungs, 'UniformOutput', false);
        X = [X{:}];
        Y = repmat(1:N, 1, rungs);
        figure
        scatter(X, Y, (1 - charge_occupancies)*500, repmat(-1, 1, length(charge_occupancies)));
        if ~SU2
            figure
            X = arrayfun(@(x) repmat(x, 1, N), 1:rungs, 'UniformOutput', false);
            X = [X{:}];
            Y = repmat(1:N, 1, rungs);
            colors = zeros(length(X), 3);
            for i = 1:length(X)
                if sign(spin_occupancies(i)) == 1
                    colors(i,:) = [0 1 0];
                elseif sign(spin_occupancies(i)) == -1
                    colors(i,:) = [1 0 0];
                else
                    error('neither 1 or -1');
                end
            end
            scatter(X, Y, (abs(spin_occupancies))*200, colors);
        end
    end
end