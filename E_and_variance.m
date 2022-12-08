%load('XXZ_Cylinder_3_5_0_final.mat');
%load('XXZ_Cylinder_idmrg2_3_maxdim_70_cut_3.7_stagh_0.05_final.mat');
%load('XXZ_Cylinder_idmrg2_4_maxdim_70_cut_3.7_stagh_0_VUMPS_final.mat')
%load('XXZ_Cylinder_VUMPS_0_maxdim_70_cut_5_stagh_0_VUMPS_final.mat')

AC1 = gs_mps.AC(1);
AC2 = gs_mps.AC(2);
AL1 = gs_mps.AL(1);
AL2 = gs_mps.AL(2);
AR1 = gs_mps.AR(1);
AR2 = gs_mps.AR(2);

E1_A = contract(AC1, [1 2 3], conj(twist(AC1,3)), [1 4 3], H_one_site_A, [4 -1 2 -2]);
E1_A = E1_A.var.var;
E1_B = contract(AC2, [1 2 3], conj(twist(AC2,3)), [1 4 3], H_one_site_B, [4 -1 2 -2]);
E1_B = E1_B.var.var;

E1 = (E1_A + E1_B)/2;
%%
E2_A_L = contract(AL1, [1 2 3], AC2, [3 4 5], conj(AL1), [1 6 8], conj(twist(AC2,3)), [8 7 5], H_A, [2 4 6 7]);
E2_B_L = contract(AL2, [1 2 3], AC1, [3 4 5], conj(AL2), [1 6 8], conj(twist(AC1,3)), [8 7 5], H_B, [2 4 6 7]);

E2_A_R = contract(AC1, [1 2 3], AR2, [3 4 5], conj(AC1), [1 6 8], conj(twist(AR2,3)), [8 7 5], H_A, [2 4 6 7]);
E2_B_R = contract(AC2, [1 2 3], AR1, [3 4 5], conj(AC2), [1 6 8], conj(twist(AR1,3)), [8 7 5], H_B, [2 4 6 7]);

E2 = (E2_A_L + E2_B_L + E2_A_R + E2_B_R)/4;

E = E1 + E2;
%%
fprintf('Energy, as calculated with expectation values, is %s x2 = %s \n', E, 2*E);
exp = expectation_value(gs_mps, H1, gs_mps);
fprintf('Correct energy, as calculated with expation_value, is ¨%s \n', exp);

%{
if isa(O, 'InfJMpo')
    [GL, GR] = environments(O, mps1, mps2);
    H = AC_hamiltonian(O, mps1, GL, GR);
    i = length(H);
    N = size(H{i}.R.var, 2);
    H{i}.R.var = H{i}.R.var(1, N, 1);
    H{i}.O{1} = H{i}.O{1}(:, :, N, :);
    AC_ = apply(H{i}, mps1.AC(i));
    E = dot(AC_, mps2.AC(i));
%}


W1 = mpo_joint{1};
W2 = mpo_joint{2};
assert(all(size(W1) == size(W2)), 'W1 and W2 dont have the same size');
[P, ~, Q] = size(W1);

[GL, GR] = environments(H1, gs_mps, gs_mps);

GL1 = GL{1};
GL2 = GL{2};
%%
% E_new1_R correct! => enkel GR{end} gebruiken
E_new1 = contract(AC1, [3 2 1], W1(1,:,:,:), [-1 5 4 2], conj(AC1), [3 5 6], GR{2}, [1 4 6]);
E_new2 = contract(AC2, [3 2 1], W2(1,:,:,:), [-1 5 4 2], conj(AC2), [3 5 6], GR{1}, [1 4 6]);
E_new1_R = E_new1.var.var;
E_new2_R = E_new2.var.var;
disp(E_new1_R);
disp(E_new2_R);

% Doesn't work
%{
E_new1 = contract(AC1, [3 2 1], W1(1,:,:,:), [-1 5 4 2], conj(AC1), [3 5 6], GR{1}, [1 4 6]);
E_new2 = contract(AC2, [3 2 1], W2(1,:,:,:), [-1 5 4 2], conj(AC2), [3 5 6], GR{2}, [1 4 6]);
E_new1_R = E_new1.var.var;
E_new2_R = E_new2.var.var;
disp(E_new1_R);
disp(E_new2_R);
%}

%%
% GL
% E_new2_L correct! => Enkel GL{end} gebruiken
E_new1 = contract(GL{1}, [6 4 1], AC1, [1 2 3], W1(:,:,Q,:), [4 5 -1 2], conj(twist(AC1,3)), [6 5 3]);
E_new2 = contract(GL{2}, [6 4 1], AC2, [1 2 3], W2(:,:,Q,:), [4 5 -1 2], conj(twist(AC2,3)), [6 5 3]);
E_new1_L = E_new1.var.var;
E_new2_L = E_new2.var.var;

disp(E_new1_L);
disp(E_new2_L);


E_new1 = contract(GL{1}, [6 4 1], AC1, [1 2 3], W2(:,:,Q,:), [4 5 -1 2], conj(twist(AC1,3)), [6 5 3]);
E_new2 = contract(GL{2}, [6 4 1], AC2, [1 2 3], W1(:,:,Q,:), [4 5 -1 2], conj(twist(AC2,3)), [6 5 3]);
E_new1_L = E_new1.var.var;
E_new2_L = E_new2.var.var;

disp(E_new1_L);
disp(E_new2_L);
%{
%E_new1 = contract(AC1, [1 2 3], W1(1,:,:,:), [-1 -2 -3 2], conj(AC1), [1 3 -4]);
E_nefsdfqdsw2 = contract(AC2, [3 2 1], W2(1,:,:,:), [-1 5 4 2], conj(AC2), [3 5 6], GR{1}, [1 4 6]);

E_new2 = contract(AC2, [1 2 3], W2(1,:,:,:), [-1 4 5 2], conj(AC2), [1 4 6], GR{2}, [3 4 6]);
E_new1 = E_new1.var.var;
E_new2 = E_new2.var.var;
%}
fprintf('Energy with W is %s \n', (E_new1_L+E_new2_L+E_new1_R+E_new2_R)/4);

%%

test = contract(AR1, [-1 3 1], W1, [-2 4 2 3], conj(AR1), [-3 4 5], GR{2}, [1 2 5]);
testvar = test.var;
testvarr = testvar.var;

GLvar = GL{2}.var;
GLvarr = GLvar.var;


%%
%{
E_new1 = contract(GL{1}, [6 4 1], AC2, [1 2 3], W1(:,:,Q,:), [4 5 -1 2], conj(twist(AC2,3)), [6 5 3]);
E_new1_L = E_new1.var.var;
disp(E_new1_L);
%}



%%
% Variance

%top_str = contract(AL1, [-1 1 2], AC2, [2 3 -4], H_A, [1 3 -2 -3]);
top_str = contract(AL1, [-1 -2 2], AC2, [2 -3 -4]);
hoera = contract(top_str, [1 3 5 -5], GL{1}, [-1 2 1], W1, [2 -2 4 3], W1(:,:,Q,:), [4 -3 -4 5]);
var_A = contract(hoera, [1 2 3 -1 4], conj(AL1), [1 2 5], conj(AC2), [5 3 4]);
disp(var_A.var.var);

%top_str = contract(AC1, [-1 1 2], AR2, [2 3 -4], H_A, [1 3 -2 -3]);
top_str = contract(AC1, [-1 -2 2], AR2, [2 -3 -4]);
hoera = contract(top_str, [1 3 5 -5], GL{1}, [-1 2 1], W1, [2 -2 4 3], W1(:,:,Q,:), [4 -3 -4 5]);
var_A = contract(hoera, [1 2 3 -1 4], conj(AC1), [1 2 5], conj(AR2), [5 3 4]);
disp(var_A.var.var);

% INCORRECT
%{
%top_str = contract(AL1, [-1 1 2], AC2, [2 3 -4], H_A, [1 3 -2 -3]);
top_str = contract(AL1, [-1 -2 2], AC2, [2 -3 -4]);
hoera = contract(top_str, [1 3 5 -5], GL{1}, [-1 2 1], W1, [2 -2 4 3], W1, [4 -3 -4 5]);
var_A = contract(hoera, [1 2 3 4 5], conj(AL1), [1 2 6], conj(AC2), [6 3 7], GR{1}, [5 4 7]);
disp(var_A);

% top_str = contract(AC1, [-1 1 2], AR2, [2 3 -4], H_A, [1 3 -2 -3]);
top_str = contract(AC1, [-1 -2 2], AR2, [2 -3 -4]);
hoera = contract(top_str, [1 3 5 -5], GL{1}, [-1 2 1], W1, [2 -2 4 3], W1, [4 -3 -4 5]);
var_A = contract(hoera, [1 2 3 4 5], conj(AC1), [1 2 6], conj(AR2), [6 3 7], GR{1}, [5 4 7]);
disp(var_A);
%}
%var_A_new = contract(AL1, [1 2 3], AC2, [3 4 5], H_A, [2 4 7 9], GL{1}, [10 6 1], W1, [6 11 8 7], W1, [8 13 14 9], conj(AL1), [10 11 12], conj(AC2), [12 13 15], GR{1}, [5 14 15]);

%error('fdjssqkl');
%NIEUW: alles met GR en GL (end)

%top_str = contract(AL2, [-1 1 2], AC1, [2 3 -4], H_B, [1 3 -2 -3]);
top_str = contract(AL2, [-1 -2 2], AC1, [2 -3 -4]);
hoera = contract(top_str, [1 3 5 -5], GL{2}, [-1 2 1], W1, [2 -2 4 3], W1, [4 -3 -4 5]);
var_B = contract(hoera, [1 2 3 4 5], conj(AL2), [1 2 6], conj(AC1), [6 3 7], GR{2}, [5 4 7]);
disp(var_B);

%top_str = contract(AL2, [-1 1 2], AC1, [2 3 -4], H_A, [1 3 -2 -3]);
top_str = contract(AL2, [-1 -2 2], AC1, [2 -3 -4]);

hoera = contract(top_str, [1 3 5 -5], GL{2}, [-1 2 1], W1, [2 -2 4 3], W1, [4 -3 -4 5]);
var_A = contract(hoera, [1 2 3 4 5], conj(AL2), [1 2 6], conj(AC1), [6 3 7], GR{2}, [5 4 7]);
disp(var_A);

%top_str = contract(AL2, [-1 1 2], AC1, [2 3 -4], H_A, [1 3 -2 -3]);
top_str = contract(AL2, [-1 -2 2], AC1, [2 -3 -4]);

hoera = contract(top_str, [1 3 5 -5], GL{2}, [-1 2 1], W1, [2 -2 4 3], W1(:,:,Q,:), [4 -3 -4 5]);
var_A = contract(hoera, [1 2 3 -1 4], conj(AL2), [1 2 5], conj(AC1), [5 3 4]);
disp(var_A.var.var);

%top_str = contract(AL2, [-1 1 2], AC1, [2 3 -4], H_A, [1 3 -2 -3]);
top_str = contract(AL2, [-1 -2 2], AC1, [2 -3 -4]);
hoera = contract(top_str, [1 3 5 -5], GL{2}, [-1 2 1], W1, [2 -2 4 3], W1(:,:,Q,:), [4 -3 -4 5]);
var_A = contract(hoera, [1 2 3 -1 4], conj(AL2), [1 2 5], conj(AC1), [5 3 4]);
disp(var_A.var.var);

%top_str = contract(AC2, [-1 1 2], AR1, [2 3 -4], H_A, [1 3 -2 -3]);
top_str = contract(AC2, [-1 -2 2], AR1, [2 -3 -4]);
hoera = contract(top_str, [1 3 5 -5], GL{2}, [-1 2 1], W1, [2 -2 4 3], W1(:,:,Q,:), [4 -3 -4 5]);
var_A = contract(hoera, [1 2 3 -1 4], conj(AC2), [1 2 5], conj(AR1), [5 3 4]);
disp(var_A.var.var);

ENDD
%%

E_new1 = contract(AC1, [3 2 1], W1(1,:,:,:), [-1 5 4 2], conj(AC1), [3 5 6], GR{2}, [1 4 6]);
E_new1 = contract(AC1, [4 2 1], W1(1,:,:,:), [-1 5 3 2], conj(AC1), [4 5 6], GR{2}, [1 3 6]);
disp(E_new1.var.var);

%ttest = contract(GR{1}, [1 2 3], AR2, [-1 4 1], W1(1,:,:,:), [-1 5 2 4], conj(AR2), [-2 5 3]);

%E_new1_again = contract(AC1, [1 2 3], W1(1,:,:,:), [-1 5 4 2], conj(AC1), [1 5 6], AR2, [3 7 9], conj(AR2), [6 8 11], W2, [4 8 10 7], GR{2}, [9 10 11]);
%E_new1_again = contract(AC1, [5 2 3], W1(1,:,:,:), [-1 6 1 2], conj(AC1), [5 6 7], AR2, [3 4 10], conj(AR2), [7 8 9], W2, [1 8 11 4], GR{2}, [10 11 9]);
%disp(E_new1_again.var.var);
%E_new2 = contract(GL{2}, [6 4 1], AC2, [1 2 3], W2(:,:,Q,:), [4 5 -1 2], conj(twist(AC2,3)), [6 5 3]);

%%

% CORRECT VARIANCE

top_str = contract(AL1, [-1 1 2], AC2, [2 3 -4], H_A, [1 3 -2 -3]);
%top_str = contract(AL1, [-1 -2 2], AC2, [2 -3 -4]);
hoera = contract(top_str, [1 3 5 -5], GL{1}, [-1 2 1], W1, [2 -2 4 3], W2(:,:,Q,:), [4 -3 -4 5]);
var_A = contract(hoera, [1 2 3 -1 4], conj(AL1), [1 2 5], conj(AC2), [5 3 4]);
disp(var_A.var.var);
