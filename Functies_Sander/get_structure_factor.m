function s = get_structure_factor(O, gs_mps, max_dist, operator_type, kwargs)
    arguments
        O
        gs_mps
        max_dist
        operator_type
        kwargs.corr = []
        kwargs.convention = 'conventional'
        kwargs.system = {'1D'}
        kwargs.qx_size = 500
        kwargs.tol_check = false
    end
    % returns the static structure factor, based on 
    % equation (41-44) in Tangent-space methods...
    % The infinite sum is truncated after a certain distance,
    % or when the correlations drop below a certain treshold.
    % O is the operator, gotten from the contraction of O_alpha
    % and O_beta.
    % operator_type is either separate, joint, or twosite
    % separate means that 2 one-site operators are given, not connected with a
    % joint means that 2 one-site operators are given, connected with a virtual leg
    % twosite means a two-site operator is given
    % virtual leg.
    qx_list = linspace(0, 2*pi, kwargs.qx_size);

    if ~isempty(kwargs.corr)
        corr_list = kwargs.corr;
        max_dist = length(corr_list);
    else
        if kwargs.tol_check
            corr_list = correlation_function(O, gs_mps, max_dist, operator_type, kwargs.convention, tol_check);
        else
            corr_list = correlation_function(O, gs_mps, max_dist, operator_type, kwargs.convention);
        end
    end
    zero_dist = zero_distance_correlation_function(O, gs_mps, operator_type, kwargs.convention, 'swap', true);
    
    if strcmp(kwargs.system{1}, '1D')
        % Implements 1D system
        s = repmat(zero_dist, 1, kwargs.qx_size);
        for q_ind = 1:kwargs.qx_size
            q = qx_list(q_ind);
            for n = 1:max_dist
                s(q_ind) = s(q_ind) + exp(1i * q * n)*corr_list(n);
            end
        end
    else
        N = kwargs.system{2};
        % implements 2D system with radius N
        % This returns a cell, index i of the cell contains a list 
        % with the same size as q_lijst, this corresponds to the structure
        % factor with qy = 2pi/N * (i-1) and qx = qx_list
        
        %corr_list = [zero_dist corr_list];
        %corr_list = [1.5 corr_list];
        for ny_i = N-1:-1:0
            s{ny_i+1} = repmat(zero_dist, 1, kwargs.qx_size);
            for q_ind = 1:kwargs.qx_size

                qx = qx_list(q_ind);
                if strcmp('Helix', kwargs.system{1}) || strcmp('helix', kwargs.system{1})
                    qy = (2*pi/N)*(ny_i + qx/(2*pi));
                elseif strcmp('Cylinder', kwargs.system{1}) || strcmp('cylinder', kwargs.system{1})
                    qy = (2*pi/N)*ny_i;
                else
                    error('kwargs.system{1} must be either Helix, Cylinder or 1D')
                end
                for n = 1:max_dist
                    s{ny_i+1}(q_ind) = s{ny_i+1}(q_ind) + exp(1i * (qx * floor(n/N) + qy * mod(n,N)))*corr_list(n);
                end
            end
        end
    end
end