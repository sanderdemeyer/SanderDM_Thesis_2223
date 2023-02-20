function [s, rel_error, rel_error_transpose] = get_2D_structure_factor(O, gs_mps, N, qx_list, qy_list, max_dist, kwargs)
    arguments
        O
        gs_mps
        N
        qx_list
        qy_list
        max_dist
        kwargs.corr = []
        kwargs.colorplot = false
        kwargs.threedimplot = false
        kwargs.check_XY_symmetry = false
        kwargs.check_XY_symmetry_transpose = false
    end
    warning('Not sure about physical meaning');
    if ~isempty(kwargs.corr)
        corr_list = kwargs.corr;
        max_dist = length(corr_list);
    else
        corr_list = correlation_function(O, gs_mps, max_dist);
    end
    zero_dist = zero_distance_correlation_function(O, gs_mps);

    qx_max = length(qx_list);
    qy_max = length(qy_list);
    s = repmat(zero_dist, qx_max, qy_max);
    for i = 1:qx_max
        for j = 1:qy_max
            qx = qx_list(i);
            qy = qy_list(j);
            for n = 1:max_dist
                s(i, j) = s(i, j) + exp(1i * (qx * floor(n/N) + qy * mod(n,N)))*corr_list(n);
            end
        end
    end
    
    if kwargs.check_XY_symmetry
        assert(max_dist >= N^2)
        errors = zeros(1, N);
        for i = 1:N
            errors(i) = abs(corr_list(i) - corr_list(i*N));
        end
        rel_error = errors;
    end

    if kwargs.check_XY_symmetry_transpose
        if qx_max == qy_max
            rel_error_transpose = norm(s - transpose(s), 'fro')/(2*norm(s, 'fro')*qx_max);

        else
            error('TBA');
        end
    end

    if kwargs.colorplot
        figure
        pcolor(qx_list, qy_list, real(s));
    end
    
    if kwargs.threedimplot
        figure
        X = cell2mat(arrayfun(@(x) repmat(x, 1, qy_max), qx_list, 'UniformOutput', false));
        Y = repmat(qy_list, 1, qx_max);
        Z = reshape(s, 1, qx_max*qy_max);
        scatter3(X, Y, Z);
    end

end