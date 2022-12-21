function s = get_structure_function(O, gs_mps, q_list, range, dimension)
    % returns the static structure function, based on 
    % equation (41-44) in Tangent-space methods...
    % The infinite sum is truncated after a certain distance,
    % or when the correlations drop below a certain treshold.
    % O is the operator, gotten from the contraction of O_alpha
    % and O_beta.
    if dimension == 2
        N = period(gs_mps)/2;
        corr_list = correlation_function(O, gs_mps, range);
        corr_list = [1.5 corr_list];
        for ny_i = N-1:-1:0
            qy = (2*pi/N)*ny_i;
            s{ny_i+1} = zeros(1, length(q_list));
            for q_ind = 1:length(q_list)
                qx = q_list(q_ind);
                for n = 0:range
                    s{ny_i+1}(q_ind) = s{ny_i+1}(q_ind) + exp(1i * (qx * floor(n/N) + qy * mod(n,N)))*corr_list(n+1);
                end
            end
        end
    elseif dimension == 1
        corr_list = correlation_function(O, gs_mps, range);
        s = zeros(1, length(q_list));
        for q_ind = 1:length(q_list)
            q = q_list(q_ind);
            for n = 1:range
                s(q_ind) = s(q_ind) + exp(1i * q * n)*corr_list(n);
            end
        end
    else
        error('Dimension should be 1 or 2');
    end
end