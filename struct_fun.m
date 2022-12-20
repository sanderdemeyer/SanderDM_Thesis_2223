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

plot(q_lijst/pi, s4);
hold on
plot(q_lijst/pi, s6);
hold on
plot(q_lijst/pi, s8);
hold on

legend('N = 4', 'N = 6', 'N = 8');
xlabel('$q/ \pi$', 'interpreter', 'latex')
ylabel('Structure function $s^{\alpha \beta}(q)$', 'interpreter', 'latex')
title('Static structure function of the Heisenberg model with $\Delta = 1$ on a cylinder', 'interpreter', 'latex')

%%

tens = Tensor([pspace' pspace'], [pspace' pspace']);
%swap = num2cell([0 1 0 0 1 0]); actual swap
swap = num2cell([1 0 -1 -1 0 1]);
swap_tens = fill_tensor(tens, swap);

q_lijst = linspace(0, 2*pi, 500);

s4_swapped = get_structure_function(swap_tens, gs_mps, q_lijst, 500, 2, 8);

%%
for i = 1:8
    plot(q_lijst, s4_swapped{i})
    hold on
end
legend('0', '1', '2', '3', '4', '5', '6', '7');
xlabel('$q/ \pi$', 'interpreter', 'latex')
ylabel('Structure function $s^{\alpha \beta}(q)$', 'interpreter', 'latex')
title('Static structure function of the Heisenberg model with $\Delta = 1$ on a cylinder of N = 8, for $S_z^{i} S_z^{i+1}$', 'interpreter', 'latex')

hold off