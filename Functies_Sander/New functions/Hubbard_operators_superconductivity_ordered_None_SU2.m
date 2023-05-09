function [L1, L2, L3, L4] = Hubbard_operators_superconductivity_ordered_None_SU2()
    % assumes order c^dagger, c^dagger, c, c
    [pspace] = get_spaces_Hubbard_None_SU2();

    pos_charge = ProductCharge(SU2(2), fZ2(1));
    pos_charge_space = GradedSpace.new(pos_charge, 1, false);

    double_charge_singlet = ProductCharge(SU2(1), fZ2(0));
    double_charge_singlet_space = GradedSpace.new(double_charge_singlet, 1, false);

    Tens1 = Tensor(pspace, [pspace pos_charge_space]);
    Tens2_s = Tensor([pos_charge_space pspace], [pspace double_charge_singlet_space]);
    Tens3_s = Tensor([double_charge_singlet_space pspace], [pspace pos_charge_space]);
    Tens4 = Tensor([pos_charge_space pspace], pspace);
    
    data_1 = {reshape([0 sqrt(2)], 2, 1) reshape([1 0], 1, 1, 2)};
    data_2 = {reshape([sqrt(2) 0], 1, 1, 1, 2) [0 1]};
    data_3 = {[-sqrt(2) 0] reshape([0 -1], 1, 1, 1, 2)};
    data_4 = {reshape([0 sqrt(2)], 1, 1, 2) [1 0]};

    Tens1 = fill_tensor(Tens1, data_1);
    Tens2_s = fill_tensor(Tens2_s, data_2);
    Tens3_s = fill_tensor(Tens3_s, data_3);
    Tens4 = fill_tensor(Tens4, data_4);

    %{
    dTens1 = double(Tens1);
    dTens2_s = double(Tens2_s);
    dTens3_s = double(Tens3_s);
    dTens4 = double(Tens4);
    %}

    % equivalence: 1 = down. 2 = up.
    % Tens1(:,1,:) is c^dagger_1
    % Tens1(:,2,:) is c^dagger_2
    % Tens2_s(1,:,1,:) is c^dagger_2
    % Tens_s(2,:,1,:) is -c^dagger_1
    % Tens3_s(1,:,1,:) is c_2
    % Tens3_s(1,:,2,:) is -c_1
    % Tens4(1,:,:) is c_1
    % Tens4(2,:,:) is c_2

    Delta_s = contract(Tens1, [-1 1 -8], Tens2_s, [1 -2 2 -7], Tens3_s, [2 -3 3 -6], Tens4, [3 -4 -5]);

    Delta_tot = Delta_s + tpermute(conj(Delta_s), [8 7 6 5 4 3 2 1]);
   
    [L1, L2, L3, L4] = tsvd_foursite_conventional(Delta_tot);
   end