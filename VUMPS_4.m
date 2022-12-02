disp('Code started running');
load('Hubbard_1D_half_filling_VUMPS_4.mat');
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
naam = 'Hubbard_1D_half_filling_VUMPS_4';
naam = 'ignore';
%alg = IDmrg2('dynamical_tols', true, 'which', 'smallestreal', 'tol_min', 10^(-7), 'tol_max', 10^(-5), 'trunc', {'TruncBelow', 10^(-cut)}, 'tol', 10^(-6), 'maxiter', 100, 'verbosity', Verbosity.iter, 'name', naam, 'doSave', true, 'saveIterations', 1);
%alg = Vumps('which', 'smallestreal', 'maxiter', 200, 'verbosity', Verbosity.iter, 'dynamical_tols', false, 'tol', 10^(-6)); %'doSave', true, 'name', naam);

%[gs_mps, gs_energy] = fixedpoint(alg, H1, mps);

save(naam + '_finaljklmjklmj');
disp('Done, Hooray!');

%%

AC1 = mps.AC(1);
AC2 = mps.AC(2);
AL1 = mps.AL(1);
AL2 = mps.AL(2);


E1 = contract(AL1, [1 2 3], AC2, [3 4 5], conj(AL1), [1 6 7], conj(twist(AC2,3)), [7 8 5], H, [2 4 6 8]);
E2 = contract(AL2, [1 2 3], AC1, [3 4 5], conj(AL2), [1 6 7], conj(twist(AC1,3)), [7 8 5], H, [2 4 6 8]);
E = (E1+E2)/2;
disp(E);

%%
[GL, GR] = environments(H1, mps, mps);
%%
E1 = contract(AL1, [1 2 3], AC2, [3 4 5], conj(AL1), [1 6 7], conj(twist(AC2,3)), [7 8 5], H, [2 4 6 8]);
E2 = contract(AL2, [1 2 3], AC1, [3 4 5], conj(AL2), [1 6 7], conj(twist(AC1,3)), [7 8 5], H, [2 4 6 8]);
E = (E1+E2)/2;

E1_e = contract(GL{1}, [1 3 5], conj((AC1)), [1 2 9], mpo, [3 2 8 6], AC1, [5 6 7], GR{2}, [7 8 9]);
E2_e = contract(GL{2}, [1 3 5], conj((AC2)), [1 2 9], mpo, [3 2 8 6], AC2, [5 6 7], GR{1}, [7 8 9]);
E_e = (E1_e + E2_e)/6;

H_squared1 = contract(GL{1}, [1 7 8], conj(AL1), [1 2 3], conj(AC2), [3 4 15], H, [5 6 2 4], mpo, [7 5 11 9], AL1, [8 9 10], AC2, [10 12 13], mpo, [11 6 14 12], GR{1}, [13 14 15]);
H_squared2 = contract(GL{2}, [1 7 8], conj(AL2), [1 2 3], conj(AC1), [3 4 15], H, [5 6 2 4], mpo, [7 5 11 9], AL2, [8 9 10], AC1, [10 12 13], mpo, [11 6 14 12], GR{2}, [13 14 15]);
H_squared = (H_squared1 + H_squared2)/2;



