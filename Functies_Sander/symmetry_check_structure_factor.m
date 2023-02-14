function [error, errors, error_new] = symmetry_check_structure_factor(s)
    N = length(s);
    len = length(s{1});
    errors = zeros(N, N);
    
    for i = 0:N-1
        for j = 0:N-1
           qxi = (j + i/N)/(1-1/N^2)*(2*pi/N); % qx of i and qy of j
           qxj = (2*pi*i)/N + qxi/N; % qx of j and qy of i
           qxi_index = 1 + floor((len-1)*qxi/(2*pi));
           qxj_index = 1 + floor((len-1)*qxj/(2*pi));
           A1 = s{i+1}(qxi_index);
           A2 = s{j+1}(qxj_index);
           error_ij = abs((A1 - A2)/sqrt(A1*A2));
           errors(i+1, j+1) = error_ij;
        end
    end
    error = norm(errors)/N;
    error_new = 0;
    for i = 1:N-1
        error_new = error_new + errors5(i,i+1)^2;
    end
    error_new = sqrt(error_new);
end