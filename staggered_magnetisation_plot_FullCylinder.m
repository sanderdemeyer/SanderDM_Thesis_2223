N = 4;
delta = 1;
trunc = 500;

[gs_mps, gs_energy] = FullCylinder(4, trunc, [40 8 5], -1, 0, 0, 0);
for i = 1:10
    stag_h_field = 0.001*i;
    fprintf('Iteration %d just started, with h_stag = %d \n', i, stag_h_field)
    starting_name = 'XXZ_FullCylinder_vumps_' + string(N) + '_delta_' + string(delta) + '_trunctotdim_' + string(trunc) + '_stagh_' + string((i-1)*0.001);
    [gs_mps, gs_energy] = FullCylinder(4, trunc, [40 8 5], -1, stag_h_field, starting_name, 2);
end


