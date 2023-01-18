function [gs_mps, gs_energy, eta] = doVumps(H1, mps, vumps_way, maxiter, trunc, tol, name_base)
    trunc_tot = ~iscell(trunc);
    if trunc_tot
        trunc_way = {'TruncTotalDim', trunc};
        name = name_base + '_trunctotdim_' + string(trunc);
    else
        trunc_way = {'TruncBelow', 10^(-trunc{2}), 'TruncDim', trunc{1}};
        name = name_base + '_truncbond_' + string(trunc{1}) + '_cut_' + string(trunc{2});
    end

    if vumps_way == 1
        alg1 = Vumps('which', 'smallestreal', 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'doSave', true, 'name', strcat(name, '.mat'), 'tol', 10^(-tol), 'doplot', true);
        [gs_mps, gs_energy, ~, ~, eta] = fixedpoint(alg1, H1, mps);
    elseif vumps_way == 2
        alg2 = Vumps2('which', 'smallestreal', 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'doSave', true, 'trunc', trunc_way, 'name', strcat(name, '.mat'), 'tol', 10^(-tol), 'doplot', true);
        [gs_mps, gs_energy, ~, ~, eta] = fixedpoint(alg2, H1, mps);
    elseif vumps_way == 3
         alg = IDmrg2('dynamical_tols', true, 'which', 'smallestreal', 'trunc', trunc_way, 'tol', 10^(-tol), 'maxiter', maxiter, 'verbosity', Verbosity.iter, 'name', strcat(name, '.mat'), 'doSave', true, 'saveIterations', 1, 'doplot', true);
        [gs_mps, gs_energy] = fixedpoint(alg, H1, mps);
    else
        iterations = length(maxiter);
        gs_mps = mps;
        for i = 1:iterations
            if mod(i,2) == 1
                fprintf('Big iteration %d of %d \n', i, iterations);
                alg2 = Vumps2('which', 'smallestreal', 'miniter', 1, 'maxiter', maxiter(i), 'verbosity', Verbosity.iter, 'doSave', true, 'trunc', trunc_way, 'name', strcat(name, '.mat'), 'tol', 10^(-tol), 'doplot', true);
                [gs_mps, gs_energy, ~, ~, eta] = fixedpoint(alg2, H1, gs_mps);
            else
                fprintf('Big iteration %d of %d \n', i, iterations);
                alg1 = Vumps('which', 'smallestreal', 'miniter', 1, 'maxiter', maxiter(i), 'verbosity', Verbosity.iter, 'doSave', true, 'name', strcat(name, '.mat'), 'tol', 10^(-tol), 'doplot', true);
                [gs_mps, gs_energy, ~, ~, eta] = fixedpoint(alg1, H1, gs_mps);
            end
            if eta < 10^(-tol)
                return
            end
        end
    end
    save(strcat(name, '_final.mat'));
    disp('Done, Hooray!');
end