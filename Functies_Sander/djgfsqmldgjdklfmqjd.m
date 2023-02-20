%[mps_1, energy_1] = Hubbard_1D(1, 5, 1, 1, 20, [5 15], 7, -1, 0, 0, 'convention', 'first', 'symmetries', 'U1_U1');

[mps_2, energy_2] = Hubbard_1D(1, 5, 1, 1, 10, [4 15], 7, -1, 0, 0, 'convention', 'conventional', 'symmetries', 'U1_SU2');

[mps_3, energy_3] = Hubbard_1D(1, 5, 1, 1, 10, [4 15], 7, -1, 0, 0, 'convention', 'conventional', 'symmetries', 'None_U1');