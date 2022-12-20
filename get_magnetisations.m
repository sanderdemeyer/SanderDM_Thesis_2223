%naam = 'XXZ_Cylinder_4_5.mat';
%naam = 'XXZ_Cylinder_3_5_0_final.mat';
%naam = 'XXZ_Cylinder_5_cut_5_stagh_0.mat';
%load('XXZ_Cylinder_3_5_0_final.mat');
%load('XXZ_Cylinder_3_5_stagh_m0p1_final.mat');
%load('XXZ_Cylinder_4_maxdim_130_stagh_0_final.mat');
%load('XXZ_Cylinder_idmrg2_4_maxdim_70_stagh_1_final.mat');
%load('XXZ_Cylinder_idmrg2_4_maxdim_70_cut_3.7_stagh_10_final.mat')

%load(naam);
%disp(naam);
%gs_mps = canonicalize(mps, 'Order', 'rl');

disp(gs_mps.AL);
N = 0;
delta = 1;
stag_h_field = 0;
h_field = 0;

charges = U1([1 -1]);
pspace = GradedSpace.new(charges, [1 1], false);
trivspace = GradedSpace.new(U1(0), 1, false);

%H = get_hamiltonian('XXZ', delta, pspace);
%H = tpermute(H, [3 4 1 2], [2 2]);

%%
%gs_mps = canonicalize(mps, 'Order', 'rl');
AC1 = gs_mps.AC(1);
AC2 = gs_mps.AC(2);
AL1 = gs_mps.AL(1);
AL2 = gs_mps.AL(2);
%{
AC1 = mps.AC(1);
AC2 = mps.AC(2);
AL1 = mps.AL(1);
AL2 = mps.AL(2);
%}

%{
E1 = contract(AL1, [1 2 3], AC2, [3 4 5], conj(AL1), [1 6 7], conj(twist(AC2,3)), [7 8 5], H_A, [2 4 6 8]);
E2 = contract(AL2, [1 2 3], AC1, [3 4 5], conj(AL2), [1 6 7], conj(twist(AC1,3)), [7 8 5], H_B, [2 4 6 8]);
E = (E1+E2)/2;
fprintf('Energy is %s \n', E);


[V, D] = transfereigs(gs_mps, gs_mps, 5);

epsilons = zeros(1,5);
for i = 1:5
    epsilons(i) = -log(norm(D(i,i)));
end
disp(epsilons);
fprintf('Epsilon1 is %s \n', epsilons(2));
fprintf('This gives a correlation length of %s \n', 1/epsilons(2));
fprintf('Delta is %s \n', epsilons(4)-epsilons(2));
%}

magn = get_magnetisation('XXZ', gs_mps, pspace, trivspace, 1.5, false, false);
fprintf('Magnetisation is %s \n', magn)

stag_magn = get_magnetisation('XXZ', gs_mps, pspace, trivspace, 1.5, false, true);
fprintf('Staggered magnetisation is %s \n', stag_magn)
% above converges to 1

stag_magn_new = get_magnetisation('XXZ', gs_mps, pspace, trivspace, 2, false, true);
fprintf('New staggered magnetisation is %s \n', stag_magn_new)
% above converges to 2

stag_magn = get_magnetisation('XXZ', gs_mps, pspace, trivspace, 2.1, false, true);
fprintf('Staggered magnetisation is %s \n', stag_magn)
