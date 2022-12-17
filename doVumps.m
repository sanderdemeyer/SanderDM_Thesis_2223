function [gs_mps, gs_energy, eta] = doVumps(mpo, mps, naam, maxiter, tol, varargin)
    maxiter1 = maxiter(1);
    maxiter2 = maxiter(2);
    iterations = maxiter(3);
    trunc = varargin{1};
    alg1 = Vumps('which', 'smallestreal', 'miniter', 1, 'maxiter', maxiter1, 'verbosity', Verbosity.iter, 'doSave', true, 'name', strcat(naam, '.mat'), 'tol', 10^(-tol), 'doplot', true);
    alg2 = Vumps2('which', 'smallestreal', 'miniter', 1, 'maxiter', maxiter2, 'verbosity', Verbosity.iter, 'doSave', true, 'trunc', trunc, 'name', strcat(naam, '.mat'), 'tol', 10^(-tol), 'doplot', true);
    gs_mps = mps;
    for i = 1:iterations
        fprintf('Big iteration %s \n', iterations)
        [gs_mps, gs_energy, ~, ~, eta] = fixedpoint(alg2, mpo, gs_mps);
        [gs_mps, gs_energy, ~, ~, eta] = fixedpoint(alg1, mpo, gs_mps);
        if eta < 10^(-tol)
            return
        end
    end
end
