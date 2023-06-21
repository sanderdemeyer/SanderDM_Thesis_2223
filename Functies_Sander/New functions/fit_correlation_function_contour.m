function [exponent, maximum_x, maximum_corr] = fit_correlation_function_contour(corr_list, kwargs)
    arguments
        corr_list
        kwargs.connected = true
    end
    % fit the correlation function
    if kwargs.connected
        corr = log(abs(corr_list));
    else
        corr = log(abs(corr_list-corr_list(end)));
    end
    log_x = log(1:length(corr_list));

    len = length(corr_list);
    fraction = 2/3;
    start = round(fraction*len);
    %corr = corr(start:end);
    %log_x = log_x(start:end);

    bin_width = 10;
    starting_point = 5;
    end_point = 7;
    k = 1;
    lijst = [];
    
    maximum_x = zeros(1, floor((end_point-starting_point)/bin_width)-1);
    maximum_corr = zeros(1, floor((end_point-starting_point)/bin_width)-1);
    
    for i = 1:length(corr_list)
        if i > 30
            if i < 30 + bin_width*k
                lijst = [lijst corr(i)];
            else
                maximum = max(lijst);
                maximum_corr(k) = maximum;
                maximum_x(k) = find(corr == maximum);
                lijst = [];
                k = k + 1;
            end
        end
    end

    scatter(log(maximum_x), maximum_corr);

    derivatives_3d = zeros(1, length(maximum_corr)-4);
    derivatives_2d = zeros(1, length(maximum_corr)-4);
    for i = 3:length(maximum_x)-2
        hi = log(maximum_x(i+1)) - log(maximum_x(i));
        him1 = log(maximum_x(i)) - log(maximum_x(i-1));
        derivatives_3d(i-2) = (maximum_corr(i+1)-(hi/him1)^2*maximum_corr(i-1)-(1-(hi/him1)^2)*maximum_corr(i))/(hi*(1+hi/him1));
        derivatives_2d(i-2) = (maximum_corr(i+1)-maximum_corr(i-1))/(log(maximum_x(i+1))-log(maximum_x(i-1)));
    end

    %figure; scatter(1:90, derivatives_3d); hold on; scatter(1:90, derivatives_2d); legend('2', '3'); hold off;
    exponent = max(derivatives_3d);
    return
    %maximum_x = maximum_x(1:length(maximum_corr));
    zero_list = find(maximum_corr == 0);
    zero_first = zero_list(1);
    maximum_x = maximum_x(1:zero_first-1);
    maximum_corr = maximum_corr(1:zero_first-1);

    p = polyfit(maximum_x(end-3:end), maximum_corr(end-3:end), 1);
    exponent = p(1);
    if length(maximum_corr) < 2
        exponent = 0;
    end
end