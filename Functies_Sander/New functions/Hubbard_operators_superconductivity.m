function [L1, L2, L3, L4] = Hubbard_operators_superconductivity(P, Q)
    % assumes order c^dagger, c, c^dagger, c
    [pspace] = get_spaces_Hubbard_SU2(P, Q);

    pos_charge = ProductCharge(U1(Q), SU2(2), fZ2(1));
    pos_charge_space = GradedSpace.new(pos_charge, 1, false);

    zero_charge_singlet = ProductCharge(U1(0), SU2(1), fZ2(0));
    zero_charge_singlet_space = GradedSpace.new(zero_charge_singlet, 1, false);
    zero_charge_triplet = ProductCharge(U1(0), SU2(3), fZ2(0));
    zero_charge_triplet_space = GradedSpace.new(zero_charge_triplet, 1, false);

    Tens1 = Tensor(pspace, [pspace pos_charge_space]);
    Tens2_s = Tensor([pos_charge_space pspace], [pspace zero_charge_singlet_space]);
    Tens2_t = Tensor([pos_charge_space pspace], [pspace zero_charge_triplet_space]);
    Tens3_s = Tensor([zero_charge_singlet_space pspace], [pspace pos_charge_space]);
    Tens3_t = Tensor([zero_charge_triplet_space pspace], [pspace pos_charge_space]);
    Tens4 = Tensor([pos_charge_space pspace], pspace);
    
    Tens1 = fill_tensor(Tens1, {1, sqrt(2)});
    Tens2_s = fill_tensor(Tens2_s, {1, 1});
    Tens2_t = fill_tensor(Tens2_t, {1, 1});
    Tens3_s = fill_tensor(Tens3_s, {1, 1});
    Tens3_t = fill_tensor(Tens3_t, {1, 1});
    Tens4 = fill_tensor(Tens4, {1, 1});

    %{
    Delta_s = contract(Tens1, [-1 1 -4], Tens2_s, [1 -2 -5 -3]);
    Delta_s = contract(Delta_s, [-1 -2 -5 -6 1], Tens3_s, [1 -3 -7 -4]);
    Delta_s = contract(Delta_s, [-1 -2 -3 -6 -7 -8 1], Tens4, [1 -4 -5]);    
    %}

    Delta_s = contract(Tens1, [-1 1 -8], Tens2_s, [1 -2 2 -7], Tens3_s, [2 -3 3 -6], Tens4, [3 -4 -5]);
    Delta_t = contract(Tens1, [-1 1 -8], Tens2_t, [1 -2 2 -7], Tens3_t, [2 -3 3 -6], Tens4, [3 -4 -5]);

    Delta = Delta_s - Delta_t;

    Delta_tot = Delta + tpermute(conj(Delta), [8 7 6 5 4 3 2 1]);
   
    [L1, L2, L3, L4] = tsvd_foursite_conventional(Delta_tot);
    
    %{
    Delta_new = contract(L1, [-1 1 -8], L2, [1 -2 2 -7], L3, [2 -3 3 -6], L4, [3 -4 -5]);

    verschil = reshape(double(Delta), 256, 256) - reshape(double(Delta_new), 256, 256);
    disp(norm(verschil));
    disp('done');
    %}
end