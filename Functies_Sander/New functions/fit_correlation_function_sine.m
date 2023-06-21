function fit_correlation_function_sine(corr_list, kwargs)
    arguments
        corr_list
        kwargs.connected = true
        kwargs.starting_point = 400
        kwargs.end_point = 600 %length(corr_list)
    end

    if ~kwargs.connected
        corr_list = corr_list-corr_list(end);
    end

    x = (kwargs.starting_point:kwargs.end_point)';
    y = corr_list(kwargs.starting_point:kwargs.end_point)';
    y = y(end) + 1000*sin((2*pi/8)*x+1.1).*(x.^(-2.2));

    %model = @(n0, A, lambda, phi, exponent, x) n0 + A*sin((2*pi/lambda)*x+phi).*(x.^exponent);
    model = @(A, phi, lambda, exponent, x) y(end) + A*sin((2*pi/lambda)*x+phi).*(x.^exponent);
    fitfun = fittype(model);
    fittedmodel = fit(x, y, fitfun, 'StartPoint', [1000 1.1 8 -2.2], 'MaxIter', 100000, 'TolFun', 10^(-10));
    bounds = confint(fittedmodel);
    disp(fittedmodel);
    
    x_fit = linspace(kwargs.starting_point, kwargs.end_point, 1000);
    %y_fit = model(fittedmodel.n0, fittedmodel.A, fittedmodel.lambda, fittedmodel.phi, fittedmodel.exponent, x_fit);
    y_fit = model(fittedmodel.A, fittedmodel.phi, fittedmodel.lambda, fittedmodel.exponent, x_fit);
    figure;
    %scatter(1:length(corr_list), log(abs(corr_list)));
    hold on;
    scatter(x_fit, log(abs(y_fit)));
    scatter(x, log(abs(y)));
    legend('correlation function', 'fit', 'spec');
    hold off;
end