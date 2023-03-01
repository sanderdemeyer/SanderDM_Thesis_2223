function [S_sing, S_trip, value_sing, value_trip] = get_spin_operators(P, Q, gs_mps, N)
    [pspace, ~, trivspace] = get_spaces_Hubbard_SU2(P, Q);
    tens = Tensor([pspace pspace], [pspace pspace]);

    data_sing = zeros(1, 20);
    data_trip = zeros(1, 20);
    data_sing(10) = 1;
    data_trip(15) = 1;

    S_sing = fill_tensor(tens, num2cell(data_sing));
    S_trip = fill_tensor(tens, num2cell(data_trip));
    disp('starting with sing');
    corr_list_sing = correlation_function(S_sing, gs_mps, N, 'twosite', 'conventional');
    disp('starting with trip');
    corr_list_trip = correlation_function(S_trip, gs_mps, N, 'twosite', 'conventional');
    disp('done');
    value_sing = corr_list_sing(N);
    value_trip = corr_list_trip(N);
end