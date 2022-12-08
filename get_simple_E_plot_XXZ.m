deltas = linspace(4, -2, 50);
energies = zeros(1, 50);
corr_lengths = zeros(1, 50);

for i = 1:50
    disp(i);
    delta = deltas(i);
    [gs_mps, energy] = XXZ_Cylinder(0, delta, 70, 8, 100, 0,0,0);
    energies(i) = energy;
    fprintf('for delta = %s, energy is %s\n', delta, energy);
    save('energies', 'energies')

    [V, D] = transfereigs(gs_mps, gs_mps, 2);
    epsilons = zeros(1,2);
    for j = 1:2
        epsilons(j) = -log(norm(D(j,j)));
    end
    disp(epsilons);
    corr_lengths(i) = 1/epsilons(2);
    save('corr_lengths', 'corr_lengths');
end