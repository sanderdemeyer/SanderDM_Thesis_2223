function s = get_structure_function(O, gs_mps, q_list, range, N, type)
    % returns the static structure function, based on 
    % equation (41-44) in Tangent-space methods...
    % The infinite sum is truncated after a certain distance,
    % or when the correlations drop below a certain treshold.
    % O is the operator, gotten from the contraction of O_alpha
    % and O_beta.
    if N == 0
        % Implements 1D system
        corr_list = correlation_function(O, gs_mps, range);
        s = zeros(1, length(q_list));
        for q_ind = 1:length(q_list)
            q = q_list(q_ind);
            for n = 1:range
                s(q_ind) = s(q_ind) + exp(1i * q * n)*corr_list(n);
            end
        end
    else
        % implements 2D system with radius N
        % This returns a cell, index i of the cell contains a list 
        % with the same size as q_lijst, this corresponds to the structure
        % factor with qy = 2pi/N * (i-1) and qx = q_list
        corr_list = correlation_function(O, gs_mps, range);
        %corr_list = [zero_dist corr_list];
        %corr_list = [1.5 corr_list];
        for ny_i = N-1:-1:0
            s{ny_i+1} = zeros(1, length(q_list));
            for q_ind = 1:length(q_list)

                qx = q_list(q_ind);
                if strcmp('Helix', type)
                    qy = (2*pi/N)*(ny_i + qx/(2*pi));
                elseif strcmp('Cylinder', type)
                    qy = (2*pi/N)*ny_i;
                end

                for n = 1:range
                    s{ny_i+1}(q_ind) = s{ny_i+1}(q_ind) + exp(1i * (qx * floor(n/N) + qy * mod(n,N)))*corr_list(n);
                end
            end
        end
    end
end