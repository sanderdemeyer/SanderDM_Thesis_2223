function summation = get_4_site_correlation_function_new(O1, O2, mps, N, max_dist, kwargs)
    % This function was written for the function get_SC_order_parameter.
    % This is consequently not the most general implementation.
    arguments
        O1
        O2
        mps
        N
        max_dist
        kwargs.ri = [2 4] % in a more general implementation, this can be used.
        kwargs.P = 1
        kwargs.Q = 1
        kwargs.SU2
    end
    summation = 0;
    if kwargs.SU2
        [O1, O2] = Hubbard_operators('cdc_both', kwargs.P, kwargs.Q);
    else
        [O1, O2] = Hubbard_operators('cdc_both', kwargs.P, kwargs.Q);
    end
    assert(sum(O1.var.rank) == 4, 'O1 must have 4 physical legs');
    assert(sum(O2.var.rank) == 4, 'O2 must have 4 physical legs');

    AR = mps.AR;

    open_indices_tensors_2 = cell(max_dist+N, max_dist+N);
    open_indices_tensors_3 = cell(max_dist+N, max_dist+N, max_dist+N);
    open_indices_tensors_4 = cell(max_dist+N, max_dist+N, max_dist+N, max_dist+N);

    w = period(mps);
    T = cell(w, max_dist);
    for b = 1:period(mps)
        x = contract(AR(b), [-1 1 -3], conj(AR(b)), [-2 1 -4]);
        for i = 1:max_dist
            T{b, i} = x;
            x = contract(x, [-1 -2 1 2], AR(loop(b,i,w)), [1 3 -3], conj(AR(loop(b,i,w))), [2 3 -4]);
        end
    end
    for j = [1, -1, N, -N]
        for k = [1, -1, N, -N]
            factor_j = 1*(abs(j) == 1) - 2*(abs(j) == N);
            factor_k = 1*(abs(k) == 1) - 2*(abs(k) == N);

            fprintf('j == %d and k == %d', j, k);
            for i = 0:max_dist
                min_i = min([N, N + i, N + j, N + i + k]);
                index_list = [N - min_i 4; N + i - min_i 3; N + j - min_i 2; N + i + k - min_i 1] + 1;
                index_list_sorted = sortrows(index_list);
                open_indices = unique(index_list_sorted(:,1));
                len_open_indices = length(open_indices);                
                if len_open_indices == 2
                    if isempty(open_indices_tensors_2{open_indices(1), open_indices(2)})
                        disp('here, empty, 2 open');
                        Tens = get_open_indices(mps, open_indices, T);
                        open_indices_tensors_2{open_indices(1), open_indices(2)} = Tens;
                    else
                        disp('here, not empty, 2 open');
                        Tens = open_indices_tensors_2{open_indices(1), open_indices(2)};
                    end

                    indices_base = -(1:4);
                    place = find(open_indices==index_list(4,1));
                    indices_base(2*place-1) = 1;
                    Tens = contract(Tens, indices_base, O2, [-5 1 -6 -2*place+1]);

                    indices_base = -(1:6);
                    place = find(open_indices==index_list(3,1));
                    indices_base(2*place-1) = 1;
                    indices_base(5) = 1;
                    indices_base(6) = -2*place+1;
                    Tens = contract(Tens, indices_base);

                    indices_base = -(1:4);
                    place = find(open_indices==index_list(2,1));
                    indices_base(2*place-1) = 1;
                    Tens = contract(Tens, indices_base, O1, [-5 1 -6 -2*place+1]);

                    indices_base = -(1:6);
                    place = find(open_indices==index_list(1,1));
                    indices_base(2*place-1) = 1;
                    indices_base(5) = 1;
                    indices_base(6) = -2*place+1;
                    Tens = contract(Tens, indices_base);

                    Tens_final = contract(Tens, [1 1 2 2]);

                elseif len_open_indices == 3
                    if isempty(open_indices_tensors_3{open_indices(1), open_indices(2), open_indices(3)})
                        disp('here, empty, 3 open');
                        Tens = get_open_indices(mps, open_indices, T);
                        open_indices_tensors_3{open_indices(1), open_indices(2), open_indices(3)} = Tens;
                    else
                        disp('here, not empty, 3 open');
                        Tens = open_indices_tensors_3{open_indices(1), open_indices(2), open_indices(3)};
                    end

                    indices_base = -(1:6);
                    place = find(open_indices==index_list(4,1));
                    indices_base(2*place-1) = 1;
                    Tens = contract(Tens, indices_base, O2, [-7 1 -8 -2*place+1]);

                    indices_base = -(1:8);
                    place = find(open_indices==index_list(3,1));
                    indices_base(2*place-1) = 1;
                    indices_base(7) = 1;
                    indices_base(8) = -2*place+1;
                    Tens = contract(Tens, indices_base);

                    indices_base = -(1:6);
                    place = find(open_indices==index_list(2,1));
                    indices_base(2*place-1) = 1;
                    Tens = contract(Tens, indices_base, O1, [-7 1 -8 -2*place+1]);

                    indices_base = -(1:8);
                    place = find(open_indices==index_list(1,1));
                    indices_base(2*place-1) = 1;
                    indices_base(7) = 1;
                    indices_base(8) = -2*place+1;
                    Tens = contract(Tens, indices_base);

                    Tens_final = contract(Tens, [1 1 2 2 3 3]);

                elseif len_open_indices == 4
                    if isempty(open_indices_tensors_4{open_indices(1), open_indices(2), open_indices(3)})
                        disp('here, empty, 4 open');
                        Tens = get_open_indices(mps, open_indices, T);
                        open_indices_tensors_4{open_indices(1), open_indices(2), open_indices(3)} = Tens;
                    else
                        disp('here, not empty, 4 open');
                        Tens = open_indices_tensors_4{open_indices(1), open_indices(2), open_indices(3)};
                    end

                    indices_base = -(1:8);
                    place = find(open_indices==index_list(4,1));
                    indices_base(2*place-1) = 1;
                    Tens = contract(Tens, indices_base, O2, [-9 1 -10 -2*place+1]);

                    indices_base = -(1:10);
                    place = find(open_indices==index_list(3,1));
                    indices_base(2*place-1) = 1;
                    indices_base(9) = 1;
                    indices_base(10) = -2*place+1;
                    Tens = contract(Tens, indices_base);

                    indices_base = -(1:8);
                    place = find(open_indices==index_list(2,1));
                    indices_base(2*place-1) = 1;
                    Tens = contract(Tens, indices_base, O1, [-9 1 -10 -2*place+1]);

                    indices_base = -(1:10);
                    place = find(open_indices==index_list(1,1));
                    indices_base(2*place-1) = 1;
                    indices_base(9) = 1;
                    indices_base(10) = -2*place+1;
                    Tens = contract(Tens, indices_base);

                    Tens_final = contract(Tens, [1 1 2 2 3 3 4 4]);

                else
                    error('should not happen');
                end
                if i == 0
                    summation = summation - factor_j*factor_k*Tens_final;
                else
                    summation = summation - 2*  factor_j*factor_k*Tens_final;
                end
            end
        end
    end
%{
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
%}
end