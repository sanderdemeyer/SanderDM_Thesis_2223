disp('Code started running');
load('Hubbard_1D_half_filling_idmrg_5.mat');
mps = canonicalize(mps, 'Order', 'rl');
[pspace, vspace, trivspace, fusion_trees] = get_spaces('Hubbard', false, 1, 1, 12, 3);

%%

% For atomic separation = 2å of 1D hydrogen chain: 
% t = 1.5 en abs(U/t) = 6
t = 1.5;
U = 9;
mu = 0;
h_field = 0;

H = get_hamiltonian('Hubbard_external_field', fusion_trees, pspace, t, U, mu, h_field);
H = tpermute(H, [3 4 1 2], [2 2]);
mpo = get_mpo(H, 0);

H1 = InfJMpo({mpo mpo});

%mps = UniformMps.randnc(pspace, vspace);

disp('initialization correct');
%%

%alg = Vumps('which', 'smallestreal', 'miniter', 1, 'maxiter', 10, 'verbosity', Verbosity.iter);
%[gs_mps, gs_energy] = fixedpoint(alg, H1, mps);
cut = 4;
naam = 'Hubbard_1D_half_filling_VUMPS_5';
%alg = IDmrg2('dynamical_tols', true, 'which', 'smallestreal', 'tol_min', 10^(-7), 'tol_max', 10^(-5), 'trunc', {'TruncBelow', 10^(-cut)}, 'tol', 10^(-6), 'maxiter', 100, 'verbosity', Verbosity.iter, 'name', naam, 'doSave', true, 'saveIterations', 1);
alg = Vumps('which', 'smallestreal', 'maxiter', 200, 'verbosity', Verbosity.iter, 'dynamical_tols', false, 'doSave', true, 'name', naam, 'tol', 10^(-6));

[gs_mps, gs_energy] = fixedpoint(alg, H1, mps);

save(naam + '_final');
disp('Done, Hooray!');

