% to check which bitstring is correct
gs_energies = zeros(1, 256);
for i = 0:255
    fprintf('Started with i = %d \n', i);
    bitstring = get_bitstring(i);
    [gs_mps, gs_energy] = Hubbard_1D(1, 8, 1, 1, 3, [10], 7, -2, 0, 0, 'symmetries', 'None_SU2', 'bitstring', bitstring);
    fprintf('total bonddim is %d', get_tot_bonddim(gs_mps));
    if -0.73 < gs_energy && gs_energy < -0.57
        fprintf('For i = %d, the energy is \n', i);
        disp(gs_energy);
    end
    gs_energies(i+1) = gs_energy;
end