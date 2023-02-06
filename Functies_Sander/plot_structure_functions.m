figure()
for i = 1:6
    plot(q_lijst, s{i});
    hold on
end
legend('1', '2', '3', '4', '5', '6');
hold off

figure()
for i = 1:6
    plot(q_lijst, s_cyl{i});
    hold on
end
legend('1', '2', '3', '4', '5', '6');
hold off