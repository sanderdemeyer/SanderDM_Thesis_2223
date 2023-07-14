function [charge_occupancies_matrix, spin_occupancies, Cu_occupancies] = get_occupancies(gs_mps, P, Q, N, rungs, kwargs)
    arguments
        gs_mps
        P
        Q
        N
        rungs
        kwargs.symmetries = 'U1_SU2'
        kwargs.plot = true
        kwargs.scale1 = 500
        kwargs.scale2 = 200
        kwargs.Cu = false
    end
    % Using convention = conventional
    spin_occupancies = 0;
    if strcmp(kwargs.symmetries, 'U1_SU2')
        pspace = get_spaces_Hubbard_SU2(P, Q);
        onesite = Tensor(pspace, pspace);
        number_tot = fill_tensor(onesite, num2cell([0 1 2]));
    elseif strcmp(kwargs.symmetries, 'U1_U1')
        pspace = get_spaces_Hubbard_asymmetric(P, Q);
        onesite = Tensor(pspace, pspace);
        number_tot = fill_tensor(onesite, num2cell([0 1 1 2]));
        spin_tot = fill_tensor(onesite, num2cell([0 -1 1 0]));
    elseif strcmp(kwargs.symmetries, 'None_SU2')
        pspace = get_spaces_Hubbard_None_SU2();
        onesite = Tensor(pspace, pspace);
        number_tot = fill_tensor(onesite, {reshape([0 0; 0 2], 2, 2), 1});
    end

    AC = gs_mps.AC;
    w = period(gs_mps);

    charge_occupancies = zeros(1, w);
    for b = 1:w
        charge_occupancies(b) = contract(AC(b), [1 2 3], number_tot, [4 2], twist(conj(AC(b)),3), [1 4 3]);
    end

        charge_occupancies_matrix = reshape(real(charge_occupancies), N, rungs);
    
    if strcmp(kwargs.symmetries, 'U1_U1')
        spin_occupancies = zeros(1, w);
        for b = 1:w
            spin_occupancies(b) = contract(AC(b), [1 2 3], spin_tot, [4 2], twist(conj(AC(b)),3), [1 4 3]);
        end
    end

    if kwargs.plot
        X = arrayfun(@(x) repmat(x, 1, N), 1:rungs, 'UniformOutput', false);
        X = [X{:}];
        Y = repmat(1:N, 1, rungs);
        figure
        colors = zeros(length(X), 3);
        for i = 1:length(X)
            if sign(real(1-charge_occupancies(i))) == -1
                colors(i,:) = [0 1 0];
            elseif sign(real(1-charge_occupancies(i))) == 1
                colors(i,:) = [1 0 0];
            else
                error('neither 1 or -1');
            end
        end
        scatter(X, Y, abs(1-charge_occupancies)*kwargs.scale1, colors);
        title('Charge occupancies');
        xlabel('X');
        ylabel('Y');
        if strcmp(kwargs.symmetries, 'U1_U1')
            figure
            colors = zeros(length(X), 3);
            for i = 1:length(X)
                if sign(real(spin_occupancies(i))) == 1
                    colors(i,:) = [0 1 0];
                elseif sign(real(spin_occupancies(i))) == -1
                    colors(i,:) = [1 0 0];
                else
                    error('neither 1 or -1');
                end
            end
            scatter(X, Y, (abs(spin_occupancies))*kwargs.scale2, colors);
            title('Spin occupancies');
            xlabel('X');
            ylabel('Y');
        end
    end
    if kwargs.Cu
        assert(mod(N, 3) == 0, 'N should be a multiple of 3 to be interpreted as a threeband model');
        Cu_occupancies = zeros(N/3, rungs);
        for i = 0:N/3-1
            Cu_occupancies(i+1,:) = charge_occupancies_matrix(3*i + 1,:);
        end
        O_occupancies = zeros(2*N/3, rungs);
        for i = 0:N/3-1
            O_occupancies(2*i+1,:) = charge_occupancies_matrix(3*i + 2,:);
            O_occupancies(2*i+2,:) = charge_occupancies_matrix(3*i + 3,:);
        end
    end
end