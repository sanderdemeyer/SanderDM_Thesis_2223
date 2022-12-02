N = 6;
delta = 1;
stag_h_field = 0;
shear = 5;
charges = U1([1 -1]);
fusion_trees = FusionTree.new([2 2], charges, false, charges, false, charges, false, charges, false);
pspace = GradedSpace.new(charges, [1 1], false);
trivspace = GradedSpace.new(U1(0), 1, false);

%%
H = get_hamiltonian('XXZ', delta, pspace);
H = tpermute(H, [3 4 1 2], [2 2]);
H_one_site = get_hamiltonian('one_site_XXZ', pspace, trivspace, stag_h_field, shear);
H_one_site_A = H_one_site{1};
H_one_site_B = H_one_site{2};

mpo = get_mpo(H, N);
mpo_joint = {mpo mpo};

%mpo_joint{1}(1, 1, N+3, 1) = tpermute(H_one_site_A, [4 3 2 1], [2 2]);
%mpo_joint{2}(1, 1, N+3, 1) = tpermute(H_one_site_B, [4 3 2 1], [2 2]);

%mpo_joint{1}(1, 1, N+3, 1) = tpermute(H_one_site_A, [4 1 2 3], [2 2]);
%mpo_joint{2}(1, 1, N+3, 1) = tpermute(H_one_site_B, [4 1 2 3], [2 2]);

mpo_joint{1}(1, 1, N+3, 1) = tpermute(H_one_site_A, [2 1 4 3], [2 2]);
mpo_joint{2}(1, 1, N+3, 1) = tpermute(H_one_site_B, [2 1 4 3], [2 2]);

H1 = InfJMpo(mpo_joint);

%%
D = 12;
vspace1 = GradedSpace.new(U1([-1 1]), [D D], false);
vspace2 = GradedSpace.new(U1([2, 0, -2]), [D/2 D/2 D/2], false);
mps = UniformMps.randnc([pspace pspace], [vspace1 vspace2]);

%%
alg = IDmrg2('dynamical_tols', true, 'which', 'smallestreal', 'tol', 10^(-3), 'tol_max', 10^(-3), 'trunc', {'TruncBelow', 10^(-5)}, 'maxiter', 200, 'verbosity', Verbosity.iter);
%alg = Vumps('which', 'smallestreal', 'maxiter', 2, 'verbosity', Verbosity.iter);
[gs_mps, gs_energy] = fixedpoint(alg, H1, mps);

%%

%%

magn = get_magnetisation(gs_mps, pspace, trivspace, 2, true);

[prob_up, prob_down] = get_probabilities(gs_mps, pspace, trivspace, 2);
N = 6;
delta = 1;
stag_h_field = 0.5;
shear = 0;
charges = U1([1 -1]);
fusion_trees = FusionTree.new([2 2], charges, false, charges, false, charges, false, charges, false);
pspace = GradedSpace.new(charges, [1 1], false);
trivspace = GradedSpace.new(U1(0), 1, false);

%%
H = get_hamiltonian('XXZ', delta, pspace);
H = tpermute(H, [3 4 1 2], [2 2]);
H_one_site = get_hamiltonian('one_site_XXZ', pspace, trivspace, stag_h_field, shear);
H_one_site_A = H_one_site{1};
H_one_site_B = H_one_site{2};

mpo = get_mpo(H, N);
mpo_joint = {mpo mpo};

%mpo_joint{1}(1, 1, N+3, 1) = tpermute(H_one_site_A, [4 3 2 1], [2 2]);
%mpo_joint{2}(1, 1, N+3, 1) = tpermute(H_one_site_B, [4 3 2 1], [2 2]);

%mpo_joint{1}(1, 1, N+3, 1) = tpermute(H_one_site_A, [4 1 2 3], [2 2]);
%mpo_joint{2}(1, 1, N+3, 1) = tpermute(H_one_site_B, [4 1 2 3], [2 2]);

mpo_joint{1}(1, 1, N+3, 1) = tpermute(H_one_site_A, [2 1 4 3], [2 2]);
mpo_joint{2}(1, 1, N+3, 1) = tpermute(H_one_site_B, [2 1 4 3], [2 2]);

H1 = InfJMpo(mpo_joint);

%%
D = 12;
vspace1 = GradedSpace.new(U1([-1 1]), [D D], false);
vspace2 = GradedSpace.new(U1([2, 0, -2]), [D/2 D/2 D/2], false);
mps = UniformMps.randnc([pspace pspace], [vspace1 vspace2]);

%%
%alg = IDmrg2('dynamical_tols', true, 'which', 'smallestreal', 'tol', 10^(-3), 'tol_max', 10^(-3), 'trunc', {'TruncBelow', 10^(-5)}, 'maxiter', 200, 'verbosity', Verbosity.iter);
alg = Vumps('which', 'smallestreal', 'maxiter', 5, 'verbosity', Verbosity.iter);
[gs_mps, gs_energy] = fixedpoint(alg, H1, mps);

%%

%%

magn = get_magnetisation(gs_mps, pspace, trivspace, 2, true);

[prob_up, prob_down] = get_probabilities(gs_mps, pspace, trivspace, 2);
