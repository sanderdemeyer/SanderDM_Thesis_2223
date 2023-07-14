x = cell2mat(deltas(end-3:end))';
y = (amplitudes(end-3:end))';

x = [deltas{6} deltas{7} deltas{8} deltas{9}]';
y = [amplitudes(6) amplitudes(7) amplitudes(8) amplitudes(9)]';



model =  @(f,a,x) f + a*x;
fitfun = fittype(model);
fittedmodel = fit(x, y, fitfun);
bounds = confint(fittedmodel);
disp(fittedmodel);

x_fit = linspace(0, max(x), 1000);
y_fit = model(fittedmodel.f, fittedmodel.a, x_fit);
%y_fit_min = model(bounds(1,1), bounds(1,2), bounds(2,3), x_fit);
%y_fit_plus = model(bounds(2,1), bounds(2,2), bounds(1,3), x_fit);

figure; 
hold on;
scatter(x, y);
plot(x_fit, y_fit);
legend('Data', 'Fit');
hold off;

%%
mean(occ_O)
mean(occ_Cu)
(max(occ_O)-min(occ_O))/mean(occ_O)
(max(occ_Cu)-min(occ_Cu))/mean(occ_Cu)


