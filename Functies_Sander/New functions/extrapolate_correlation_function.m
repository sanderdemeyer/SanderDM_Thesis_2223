true_corr_pair = zeros(1, max_dist);
start_point_extrapolation = 2;

for i = 1:max_dist
    corr_values = zeros(1, length(phi_yy_list_list));
    for k = 1:length(phi_yy_list_list)
        %corr_values(k) = phi_yy_list_list{k}(i);
        corr_values(k) = log(phi_yy_list_list{k}(i));
    end
  
    %p = polyfit(cell2mat(deltas_42(start_point_extrapolation:end)), corr_values(start_point_extrapolation:end),1);
    p = polyfit(cell2mat(deltas_42(start_point_extrapolation:end)), log(corr_values(start_point_extrapolation:end)),1);
    x = abs(log(cell2mat(epsilon1s(start_point_extrapolation:end))))';
    y = abs(log(abs(corr_values(start_point_extrapolation:end))))';
    p = polyfit(x,y,1);


    true_corr_pair(i) = exp(p(2));
    if i == 200
        model =  @(f,a,b,x) f + a*x.^b;
        fitfun = fittype(model);
        fittedmodel = fit(x, y, fitfun);
        disp('ok');
    end
end


true_corr_dens = zeros(1, max_dist);

for i = 1:max_dist*2
    corr_values = zeros(1, length(corr_list_density_list));
    for k = 1:length(phi_yy_list_list)
        %corr_values(k) = corr_list_density_list{k}(i);
        corr_values(k) = log(corr_list_density_list{k}(i)-corr_list_density_list{k}(end));
    end
  
    p = polyfit(cell2mat(deltas_42(start_point_extrapolation:end)), corr_values(start_point_extrapolation:end),1);
    true_corr_dens(i) = p(2);
    if i == 250
        disp('ok');
    end
end

figure;
scatter(log(1:2*max_dist), true_corr_dens);
hold on;
scatter(log(1:max_dist), true_corr_pair);
legend('density', 'pair-pair');
title('extrapolated correlation functions');
hold off;

x = cell2mat(deltas_42(3:end));
y = cell2mat(exp_SC_list(3:end));
[p, S] = polyfit(x, y, 1);
[y_fit, delta] = polyval(p,x,S);
RMSE = sqrt(mean(y-y_fit).^2);
disp(RMSE);

x = cell2mat(deltas_42(3:end));
y = cell2mat(log(exp_SC_list(3:end)));
[p, S] = polyfit(x, y, 1);
[y_fit, delta] = polyval(p,x,S);
RMSE = sqrt(mean(y-y_fit).^2);
disp(RMSE);

x = cell2mat(deltas_42(3:end));
y = cell2mat(exp_SC_list(3:end));
[p, S] = polyfit(x, y, 1);
[y_fit, delta] = polyval(p,x,S);
RMSE = sqrt(mean(y-y_fit).^2);
disp(RMSE);

%%

figure;
scatter(log(1:2*max_dist), log(true_corr_dens-true_corr_dens(end)));
hold on;
scatter(log(1:max_dist), log(true_corr_pair));
legend('density', 'pair-pair');
title('extrapolated correlation functions');
hold off;


%%
legend_bonddim = cell(1,6);

figure;
for k = 6:length(phi_yy_list_list)
    scatter(log(1:max_dist), log(phi_yy_list_list{k}-phi_yy_list_list{k}(end)));
    hold on;
    legend_bonddim{k} = num2str(tot_bonddim_list{k});
end
legend(legend_bonddim);

%% Using only the exponents
figure;
hold on;
scatter(cell2mat(deltas_42), cell2mat(exp_SC_disconnected_list));
scatter(cell2mat(deltas_42), cell2mat(exp_density_list));
legend('SC', 'density');
xlabel('$$ \delta_{42} $$', 'interpreter', 'latex');
ylabel('Exponent')
title('Comparison between charge and pair-pair correlation function exponent for $$ mu = -10.1 $$ and $$ f = 0.8732 $$', 'interpreter', 'latex');
hold off;


%% Only using the exponents - SC

x = cell2mat(deltas_42(3:end))';
y = cell2mat(exp_SC_list(3:end))';

model =  @(f,a,b,x) f + a*x.^b;
fitfun = fittype(model);
fittedmodel = fit(x, y, fitfun, 'StartPoint', [0 1 0]);
bounds = confint(fittedmodel);
disp(fittedmodel);

x_fit = linspace(0, max(x), 1000);
y_fit = model(fittedmodel.f, fittedmodel.a, fittedmodel.b, x_fit);
y_fit_min = model(bounds(1,1), bounds(1,2), bounds(2,3), x_fit);
y_fit_plus = model(bounds(2,1), bounds(2,2), bounds(1,3), x_fit);

figure; 
hold on;
scatter(x, y);
plot(x_fit, y_fit);
legend('Data', 'Fit');
hold off;


%% Only using the exponents - density

x = cell2mat(deltas_42(3:end))';
y = cell2mat(exp_density_list(3:end))';

model =  @(f,a,b,x) f + a*x.^b;
fitfun = fittype(model);
fittedmodel = fit(x, y, fitfun, 'StartPoint', [0 0 1], Lower = [-10^5 -10^5 1]);
bounds = confint(fittedmodel);
disp(fittedmodel);

x_fit = linspace(0, max(x), 1000);
y_fit = model(fittedmodel.f, fittedmodel.a, fittedmodel.b, x_fit);
y_fit_min = model(bounds(1,1), bounds(1,2), bounds(2,3), x_fit);
y_fit_plus = model(bounds(2,1), bounds(2,2), bounds(1,3), x_fit);

figure; 
hold on;
scatter(x, y);
plot(x_fit, y_fit);
legend('Data', 'Fit');
hold off;

