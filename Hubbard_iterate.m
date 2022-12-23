for U = linspace(10, 2, 9)
    disp(U);
    Hubbard_1D_external_field_half_filling(1, U, 500, [75 5 1], 6, -1, 0, 0, 0, 0)
end