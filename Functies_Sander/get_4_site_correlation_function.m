function corr = get_4_site_correlation_function(O1, O2, AL, AC, AR, N, indices, max_dist, kwargs)
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
            for i = 2:max_dist
                if j > 0
                    Tens = contract(AC(b), [1 2 -1], O1, [2 -3 3 -4], conj(AC(b)), [1 3 -2]);
                    if j == i + k
                        if j < i
                            contract(Tens, [1 2 ? ? ], T{})
                        else % j can not be i here
    
                        end
                    elseif j < i + k
                        if k < 0
                        
                        elseif j < i
                         
                        elseif j == i

                        else
                        
                        end
                    else
                        if j < i

                        elseif j == i

                        elseif k < 0

                        else
                    end
                else

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