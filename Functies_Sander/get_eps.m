function [epsilon1, delta] = get_eps(gs_mps)
    disp(gs_mps.AL)
    plot_entanglementspectrum(gs_mps);
    % Returns epsilon1 and delta = epsilon2 - epsilon1
    [~, D] = transfereigs(gs_mps, gs_mps, 3);

    epsilons = zeros(1,3);
    for i = 1:3
        epsilons(i) = -log(norm(D(i,i)));
    end
    epsilon1 = epsilons(2);
    delta = epsilons(3)-epsilons(2);
    fprintf('epsilon1 is %s \n', epsilon1);
    fprintf('delta is %s \n', delta);
end