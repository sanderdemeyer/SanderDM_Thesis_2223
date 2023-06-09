mu_list = [-9.2 -8.9 -10.12 -6.9 -10.1 -15 -6.5 -9.5 -10.2 -8 -7.5 -8.5 -7 -9.8753];

exponents_SC = cell(1, 14);
exponents_dens = cell(1, 14);
phases = cell(1, 14);
for i = 1:length(mu_list)
    mu = mu_list(i);
    fprintf('Started with mu = %d\n', mu);
    [exp_SC_list, exp_density_list, phase_list] = analysis_SC_vs_density(mu);
    close;
    exponents_SC{i} = exp_SC_list;
    exponents_dens{i} = exp_density_list;
    phases{i} = phase_list;
end
save('analysis_everything');