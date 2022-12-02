values = [0.1 0.05 0.025];
for stag_h_field = values
    disp(stag_h_field);
    XXZ_Cylinder(N, 70, 3.7, 200, stag_h_field, 0, 1)
end
