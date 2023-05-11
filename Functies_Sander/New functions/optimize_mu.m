function optimize_mu(starting_mu, second_mu, target_filling, error)
    model = 'Hg_oneband1';
    N = 2;
    P = 0;
    Q = 0;
    rungs = 1;

    gs_mps = Hubbard_cylinder_oneband(2, 'Hg_oneband1', 0, 0, 1, 4, [20 20 20 50], 7, -2, 0, 0, 'symmetries', 'None_SU2', 'mu', starting_mu);
    mus{1} = starting_mu;
    fillings{1} = get_filling(gs_mps);
    prev_name = 'Hubbard_FullCylinder_oneband_' + string(model) + '_N_' + string(N) + '_P_' + string(P) + '_Q_' + string(Q) + '_rungs_' + string(rungs) + '_mu_' + string(starting_mu);

    gs_mps = Hubbard_cylinder_oneband(2, 'Hg_oneband1', 0, 0, 1, 4, [20 20 20 50], 7, -2, prev_name, 2, 'symmetries', 'None_SU2', 'mu', second_mu);
    prev_name = 'Hubbard_FullCylinder_oneband_' + string(model) + '_N_' + string(N) + '_P_' + string(P) + '_Q_' + string(Q) + '_rungs_' + string(rungs) + '_mu_' + string(second_mu);
    mus{2} = second_mu;
    fillings{2} = get_filling(gs_mps);

    k = 3;
    while abs(filling - target_filling) > error
        p = polyfit([mus{end-1} mus{end}], [fillings{end-1} fillings{end}], 1);
        mu = (target_filling - p(2))/p(1);
        gs_mps = Hubbard_cylinder_oneband(2, 'Hg_oneband1', 0, 0, 1, 4, [20 20 20 50], 7, -2, prev_name, 2, 'symmetries', 'None_SU2', 'mu', mu);

        mus{k} = mu;
        fillings{k} = get_filling(gs_mps);
        fprintf('For mu = %d, filling = %d \n', mu, fillings{k});
        k = k+1;
        prev_name = 'Hubbard_FullCylinder_oneband_' + string(model) + '_N_' + string(N) + '_P_' + string(P) + '_Q_' + string(Q) + '_rungs_' + string(rungs) + '_mu_' + string(mu);
    end
end