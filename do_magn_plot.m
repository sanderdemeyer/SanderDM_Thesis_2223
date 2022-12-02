function do_magn_plot(N, eps)
    doPath
    for i = 1:125
        stag_h_field = eps + 0.001*i;
        disp(stag_h_field);
        XXZ_Cylinder(N, 70, 3.7, 25, stag_h_field, 0, 1)
    end
end