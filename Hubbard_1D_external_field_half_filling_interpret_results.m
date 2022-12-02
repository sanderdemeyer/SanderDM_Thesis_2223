disp('Code started running');

%load('Hubbard_1D_external_field_0_half_filling.mat');
%mps = canonicalize(mps, 'Order', 'rl');

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

E2_A_L = contract(AL1, [1 2 3], AC2, [3 4 5], conj(AL1), [1 6 8], conj(twist(AC2,3)), [8 7 5], H, [2 4 6 7]);
E2_B_L = contract(AL2, [1 2 3], AC1, [3 4 5], conj(AL2), [1 6 8], conj(twist(AC1,3)), [8 7 5], H, [2 4 6 7]);


E2_A_R = contract(AC1, [1 2 3], AR2, [3 4 5], conj(AC1), [1 6 8], conj(twist(AR2,3)), [8 7 5], H, [2 4 6 7]);
E2_B_R = contract(AC2, [1 2 3], AR1, [3 4 5], conj(AC2), [1 6 8], conj(twist(AR1,3)), [8 7 5], H, [2 4 6 7]);


E1_1 = contract(AC1, [1 2 3], conj(twist(AC1,3)), [1 4 3], H_one_site, [-1 4 -2 2]);
E1_1 = E1_1.var.var;
E1_2 = contract(AC2, [1 2 3], conj(twist(AC2,3)), [1 4 3], H_one_site, [-1 4 -2 2]);
E1_2 = E1_2.var.var;

%%
%alg = Vumps('which', 'smallestreal', 'miniter', 1, 'maxiter', 10, 'verbosity', Verbosity.iter);
%[gs_mps, gs_energy] = fixedpoint(alg, H1, mps);

%alg = IDmrg2('dynamical_tols', true, 'which', 'smallestreal', 'tol_min', 10^(-6), 'tol_max', 10^(-4), 'trunc', {'TruncBelow', 10^(-6)}, 'tol', 10^(-5), 'maxiter', 200, 'verbosity', Verbosity.iter, 'name', 'Hubbard_Cylinder_Twosite_half_filling', 'doSave', true, 'saveIterations', 1);
[gs_mps, gs_energy] = fixedpoint(alg, H1, mps);


%%
deltas = zeros(1,3);
epsilon1s = zeros(1,3);
finals = [0 1 1];
maxdims = [70 100 130];
for i = 1:3
    disp(i);
    name = 'Hubbard_1D_half_filling_idmrg2_t_1.5_U_9_maxdim_' + string(maxdims(i)) + '_cut_15_redef_0';
    if finals(i)
        name = strcat(name, '_final.mat');
    else
        name = strcat(name, '.mat');
    end
    load(name)
    if ~finals(i)
        gs_mps = canonicalize(mps, 'Order', 'rl');
    end

    [V, D] = transfereigs(gs_mps, gs_mps, 5);
    epsilons = zeros(1,5);
    for j = 1:5
        epsilons(j) = -log(norm(D(j,j)));
    end
    disp(epsilons);
    epsilon1s(i) = epsilons(2);
    deltas(i) = epsilons(3) - epsilons(2);
end



%%
values = zeros(1, 20);
for i = 1:20
    values(i) = real(D(i,i));
end
disp(values);

diff_values = zeros(1, 19);
for i = 1:19
    diff_values(i) = values(i) - values(i+1);
end
disp(diff_values);
