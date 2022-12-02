for cut = 3:8
    disp(cut);
    load('XXZ_1D_2_idmrg_' + string(cut) + '_final');
    disp(gs_mps.AL);
    plot_entanglementspectrum(gs_mps);
end