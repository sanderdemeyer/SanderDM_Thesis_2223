function [exponent, maximum_x, maximum_corr] = fit_correlation_function(corr_list, kwargs)
    arguments
        corr_list
        kwargs.connected = true
        kwargs.edges = [5.4 5.95]

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

    bin_width = 0.05;
    starting_point = 5;
    end_point = 5.7;
    starting_point = 5.4;
    end_point = 5.95;
    starting_point = kwargs.edges(1);
    ending_point = kwargs.edges(2);
    %starting_point = 3.5;
    %end_point = 4.5;
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
                max_x_index = find(corr == maximum);
                if length(max_x_index) > 1
                    maximum_x(k) = log(max_x_index(1));
                else
                    maximum_x(k) = log(max_x_index);
                end
                lijst = [];
                k = k + 1;
            end
        end
    end
    %maximum_x = maximum_x(1:length(maximum_corr));
    zero_list = find(maximum_corr == 0);
    if ~isempty(zero_list)
        zero_first = zero_list(1);
        maximum_x = maximum_x(1:zero_first-1);
        maximum_corr = maximum_corr(1:zero_first-1);
    end

    p = polyfit(maximum_x(end-3:end), maximum_corr(end-3:end), 1);
    exponent = p(1);
    if length(maximum_corr) < 2
        exponent = 0;
    end
end