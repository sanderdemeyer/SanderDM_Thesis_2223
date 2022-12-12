lijstje = [10 50 100 500];
q_lijst = linspace(-20, 20, 100);
for n = 1 : length(lijstje)
    disp(n)
    leng = lijstje(n);
    disp(leng);
    s_list{n} = get_structure_function(H_A, gs_mps, q_lijst, leng, false);
end

%%

for i = 4:4
    s = s2; % s_list{i};
    plot(q_lijst/pi, (s));
    hold on
end
xlabel('$q/\pi$', 'interpreter', 'latex')
ylabel('Static structure function $s^{\alpha \beta}(q)$', 'interpreter', 'latex')
title('Static structure function for the Hamiltonian of the 1D Heisenberg model with $\Delta = 1$', 'interpreter', 'latex')
hold off
%legend('10', '30', '50', '100', '200', '500')

%%

plot(q_lijst, s4);
hold on
plot(q_lijst, s6);
hold on
plot(q_lijst, s8);

legend('N = 4', 'N = 6', 'N = 8')
xlabel('$q$', 'interpreter', 'latex')
ylabel('Structure function $s^{\alpha \beta}(q)$', 'interpreter', 'latex')
title('Static structure function for the Hamiltonian of the 1D Heisenberg model with $\Delta = 1$', 'interpreter', 'latex')
hold off
