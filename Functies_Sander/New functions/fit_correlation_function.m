function [exponent, maximum_x, maximum_corr] = fit_correlation_function(corr_list, kwargs)
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

    bin_width = 0.1;
    starting_point = 5;
    end_point = 7;
    k = 1;
    lijst = [];
    
    maximum_x = zeros(1, floor((end_point-starting_point)/bin_width)-1);
    maximum_corr = zeros(1, floor((end_point-starting_point)/bin_width)-1);
    
    for i = 1:length(corr_list)
        if log_x(i) > starting_point && log_x(i) < end_point
            if log_x(i) < starting_point + bin_width*k
                lijst = [lijst corr(i)];
            else
                maximum = max(lijst);
                maximum_corr(k) = maximum;
                maximum_x(k) = log(find(corr == maximum));
                lijst = [];
                k = k + 1;
            end
        end
    end
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