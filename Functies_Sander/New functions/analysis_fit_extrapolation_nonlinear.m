i = 6;
x = cell2mat(deltas(i:end));
y = cell2mat(epsilon1s(i:end));
x(5) = deltas{4};
y(5) = epsilon1s{4};
x = x';
y = y';

model =  @(f,a,b,x) f + a*x.^b;
fitfun = fittype(model);
fittedmodel = fit(x, y, fitfun, 'StartPoint', [0 0 1]);
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
