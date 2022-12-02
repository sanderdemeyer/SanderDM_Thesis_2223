l = length(corr);
p1 = polyfit(log(1:l), log(abs(corr)), 1);
a = p1(2);
b = p1(1);
%%

plot(1:l, abs(corr), 'Color', [0.8500 0.3250 0.0980]);
hold on

%%

scatter(1:l, abs(corr), "blue", "x");
%set(gca,'yscale', 'log')
set(gca, 'xscale','log', 'yscale', 'log')
%%
scatter(1:l, abs(tot_numb_corr_new), "blue", "x");
set(gca,'xscale','log', 'yscale', 'log')

hold on

N = 1:20;
%plot(N, 10.^(a*log(N)+b));
set(gca,'xscale','log', 'yscale', 'log')
xlabel('Distance between site i and site j');
ylabel('$Correlation: \langle n_{i,up}n_{j,down} \rangle -  \langle n_{i,up} \rangle \langle n_{j,down} \rangle$', 'interpreter', 'latex');
title('Correlation function of the 1D Hubbard Model, not redefined')
disp('hey')