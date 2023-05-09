function [L1, L2, L3, L4] = Hubbard_operators_superconductivity_ordered(P, Q)
    % assumes order c^dagger, c^dagger, c, c
    [pspace] = get_spaces_Hubbard_SU2(P, Q);

    pos_charge = ProductCharge(U1(Q), SU2(2), fZ2(1));
    pos_charge_space = GradedSpace.new(pos_charge, 1, false);

    double_charge_singlet = ProductCharge(U1(2*Q), SU2(1), fZ2(0));
    double_charge_singlet_space = GradedSpace.new(double_charge_singlet, 1, false);
    double_charge_triplet = ProductCharge(U1(2*Q), SU2(3), fZ2(0));
    double_charge_triplet_space = GradedSpace.new(double_charge_triplet, 1, false);

    Tens1 = Tensor(pspace, [pspace pos_charge_space]);
    Tens2_s = Tensor([pos_charge_space pspace], [pspace double_charge_singlet_space]);
    Tens2_t = Tensor([pos_charge_space pspace], [pspace double_charge_triplet_space]);
    Tens3_s = Tensor([double_charge_singlet_space pspace], [pspace pos_charge_space]);
    Tens3_t = Tensor([double_charge_triplet_space pspace], [pspace pos_charge_space]);
    Tens4 = Tensor([pos_charge_space pspace], pspace);
    
    % Tens1(:,1,:) is c^dagger_down
    % Tens1(:,2,:) is c^dagger up
    % Tens2_s(1,:,1,:) is c^dagger_up
    % Tens_s(2,:,1,:) is -c^dagger_down
    % Tens3_s(1,:,1,:) is c_up
    % Tens3_s(1,:,2,:) is -c_down
    % Tens4(1,:,:) is c_down
    % Tens4(2,:,:) is c_up

    Tens1 = fill_tensor(Tens1, {1, sqrt(2)});
    Tens2_s = fill_tensor(Tens2_s, {sqrt(2), 1});
    %Tens2_t = fill_tensor(Tens2_t, {sqrt(2), sqrt(3)});
    Tens3_s = fill_tensor(Tens3_s, {-sqrt(2), 1});
    %Tens3_t = fill_tensor(Tens3_t, {sqrt(2), -sqrt(3)});
    Tens4 = fill_tensor(Tens4, {1, -sqrt(2)});

    %{
    Delta_s = contract(Tens1, [-1 1 -4], Tens2_s, [1 -2 -5 -3]);
    Delta_s = contract(Delta_s, [-1 -2 -5 -6 1], Tens3_s, [1 -3 -7 -4]);
    Delta_s = contract(Delta_s, [-1 -2 -3 -6 -7 -8 1], Tens4, [1 -4 -5]);    
    %}

    Delta_s = contract(Tens1, [-1 1 -8], Tens2_s, [1 -2 2 -7], Tens3_s, [2 -3 3 -6], Tens4, [3 -4 -5]);
    %Delta_t = contract(Tens1, [-1 1 -8], Tens2_t, [1 -2 2 -7], Tens3_t, [2 -3 3 -6], Tens4, [3 -4 -5]);

    %Delta = Delta_s - Delta_t;

    Delta_tot = Delta_s + tpermute(conj(Delta_s), [8 7 6 5 4 3 2 1]);
   
    [L1, L2, L3, L4] = tsvd_foursite_conventional(Delta_tot);
   end