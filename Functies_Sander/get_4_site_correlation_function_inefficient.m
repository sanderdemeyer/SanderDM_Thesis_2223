function corr = get_4_site_correlation_function_inefficient(O1, O2, AL, AC, AR, N, indices, max_dist, kwargs)
    % This function was written for the function get_SC_order_parameter.
    % This is consequently not the most general implementation.
    arguments
        O1
        O2
        AL
        AC
        AR
        N
        indices
        max_dist
        kwargs.ri = [2 4] % in a more general implementation, this can be used.
    end
    assert(sum(O1.var.rank) == 4, 'O1 must have 4 physical legs');
    assert(sum(O2.var.rank) == 4, 'O2 must have 4 physical legs');
    disp(O1);

    x = num2cell(zeros(1, max_dist));
    T = contract(AR, [-1 1 -3], conj(AR), [-2 1 -4]);
    for i = 1:max_dist
        x{i} = T;
        T = contract(T, [-1 -2 1 2], AR, [1 3 -3], conj(AR), [2 3 -4]);
    end
    
    for j = [1, -1, N, -N]
        for k = [1, -1, N, -N]
            for i = 0:max_dist
                Tens = contract(AC, [1 -1 -3], conj(AC), [1 -2 -4]);
                for index = 2:max([0 j i i+k]) - min([0 j i i+k])
                    Tens = contract(Tens, [])
                end
            end
        end
    end

    for j = [1, -1, N, -N]
        for k = [1, -1, N, -N]
            for i = 1:max_dist
                min_i = min([N, N + j, N + k]);
                index_list = [N - min_i 4; N + i - min_i 3; N + j - min_i 2; N + i + k - min_i 1];
                index_list = sortrows(index_list);

                %initialisation
                if index_list(2, 1) == 0
                    if index_list(3, 1) == 0
                        Tens = contract(AC, [])
                    else
                        2 times
                    end
                else
                    if index_list(1, 2) < 3
                        Tens = contract(AC, [1 2 -1], O1)
                end

                for index = 1:4
                    if index_list(index,1) == index_list(index+1,1)
                        if index_list(index, 1) == index_list(index+2,1)
                            wtf
                        else
                            nogsteedswtf
                        end
                    else
                        fkjsdm
                    end
                end
            end
        end
    end
end