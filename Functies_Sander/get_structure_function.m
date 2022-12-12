function s = get_structure_function(O, gs_mps, q_list, range)
    % returns the static structure function, based on 
    % equation (41-44) in Tangent-space methods...
    % The infinite sum is truncated after a certain distance,
    % or when the correlations drop below a certain treshold.
    % O is the operator, gotten from the contraction of O_alpha
    % and O_beta.
    corr_list = correlation_function(O, gs_mps, range);
    s = zeros(1, length(q_list));
    for q_ind = 1:length(q_list)
        q = q_list(q_ind);
        for n = 1:range
            s(q_ind) = s(q_ind) + exp(1i * q * n)*corr_list(n);
        end
    end
end