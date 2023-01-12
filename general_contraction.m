%{
%N = 3;
delta = 1;
%stag_h_field = 0;

charges = U1([1 -1]);
fusion_trees = FusionTree.new([2 2], charges, false, charges, false, charges, false, charges, false);
pspace = GradedSpace.new(charges, [1 1], false);
trivspace = GradedSpace.new(U1(0), 1, false);

D = 12;
vspace1 = GradedSpace.new(U1([-1 1]), [D D], false);
vspace2 = GradedSpace.new(U1([2, 0, -2]), [D/2 D/2 D/2], false);

%%
H = get_hamiltonian('XXZ', 1, pspace);
H_one_site = get_hamiltonian('one_site_XXZ', pspace, trivspace, stag_h_field, 0);
H_one_site_A = H_one_site{1};
H_one_site_B = H_one_site{2};

mpo_joint = get_mpo(H, N, 'FullCylinder');
H1 = InfJMpo(mpo_joint);
%%

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

%a = contraction_general([AL1, [1 2 3], AC2, [3 4 5], conj(AL1), [1 6 8], conj(twist(AC2,3)), [8 7 5], H_A, [2 4 6 7]]);
a = contraction_general({AL1 [1 2 3] AC2 [3 4 5] conj(AL1) [1 6 8] conj(twist(AC2,3)) [8 7 5] H_A [2 4 6 7]});

E2_A_R = contract(AC1, [1 2 3], AR2, [3 4 5], conj(AC1), [1 6 8], conj(twist(AR2,3)), [8 7 5], H_A, [2 4 6 7]);
E2_B_R = contract(AC2, [1 2 3], AR1, [3 4 5], conj(AC2), [1 6 8], conj(twist(AR1,3)), [8 7 5], H_B, [2 4 6 7]);

E2 = (E2_A_L + E2_B_L + E2_A_R + E2_B_R)/4;

E = E1 + E2;
%}
%%

AL = gs_mps.AL;
AC = gs_mps.AC;

p = [1 2 3 4];
c = get_energy(p);
p = [2 3 4 1];
c = get_energy(p);
p = [3 4 1 2];
c = get_energy(p);
p = [4 1 2 3];
c = get_energy(p);


function c = get_energy(p)
    load('XXX_FullCylinder_idmrg2_2_maxdim_50_cut_6_stagh_0.002_final.mat')
    AL = gs_mps.AL;
    AC = gs_mps.AC;

    N = 4;
    contr_list = zeros(1,1);
    contr_list = num2cell(contr_list);
    % General Energy calculation
    n_tens = 2*N+1;
    contr_list{n_tens*2} = 0;
    contr_list{1} = AL(p(1));
    contr_list{2} = [1 2 5];
    contr_list{3} = conj(AL(p(1)));
    contr_list{4} = [1 3 4];
    for i = 1 : N-2
        contr_list{4*i+1} = AL(p(i+1));
        contr_list{4*i+2} = [4*i+1 4*i+2 4*i+5];
        contr_list{4*i+3} = conj(AL(p(i+1)));
        contr_list{4*i+4} = [4*i 4*i+3 4*i+4];
    end
    
    contr_list{4*N-3} = AC(p(N));
    contr_list{4*N-2} = [4*N-3 4*N-2 4*N]; % Normally the last one would be 4*N+1, but contracted with conj
    contr_list{4*N-1} = conj(twist(AC(p(N)),3));
    contr_list{4*N} = [4*N-4 4*N-1 4*N];
    
    
    contr_list{4*N+1} = H;
    contr_list{4*N+2} = [2 6 3 7];
    contr_list{4*N+3} = H;
    contr_list{4*N+4} = [10 14 11 15];
    
    %{
    contr_list{4*N+5} = H;
    contr_list{4*N+6} = [18 22 19 23];
    contr_list{4*N+7} = H;
    contr_list{4*N+8} = [26 30 27 31];
    %}    
    c = contraction_general(contr_list);
    disp(c);
end
