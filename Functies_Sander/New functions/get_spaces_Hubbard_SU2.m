function [pspace, vspace, trivspace, prodspace, fusion_trees] = get_spaces_Hubbard_SU2(P, Q, kwargs)
    arguments
        P
        Q
        kwargs.D = 1
    end

    if gcd(P,Q) ~= 1
        warning('P = %d and Q = %d should be relative prime, working with P = %d and Q = %d', P, Q, P/gcd(P,Q), Q/gcd(P,Q))
    end
    % Implementation of the physical charge
    a = ProductCharge(U1(-P), SU2(1), fZ2(0));
    b = ProductCharge(U1(Q-P), SU2(2), fZ2(1));
    d = ProductCharge(U1(2*Q-P), SU2(1), fZ2(0));
    charges1 = [a b d];

    % Implementation of the charges of the productspace
    e = ProductCharge(U1(-2*P), SU2(2), fZ2(0));
    f = ProductCharge(U1(Q-2*P), SU2(2), fZ2(1));
    h = ProductCharge(U1(2*Q-2*P), SU2(2), fZ2(0));

    i = ProductCharge(U1(2*Q-2*P), SU2(2), fZ2(0)); % same as h
    j = ProductCharge(U1(2*Q-2*P), SU2(2), fZ2(0)); % same as h
    k = ProductCharge(U1(3*Q-2*P), SU2(2), fZ2(1));
    l = ProductCharge(U1(2*Q-2*P), SU2(2), fZ2(0)); % same as h
    m = ProductCharge(U1(3*Q-2*P), SU2(2), fZ2(1));
    n = ProductCharge(U1(4*Q-2*P), SU2(2), fZ2(0));
    charges2 = [e f h i k l m n];

    trivcharge = ProductCharge(U1(0), SU2(1), fZ2(0));

    fusion_trees = FusionTree.new([2 2], charges1, false, charges1, false, charges1, false, charges1, false);
    pspace = GradedSpace.new(charges1, [1 1 1], false);
    trivspace = GradedSpace.new(trivcharge, 1, false);
    vspace = [pspace pspace*pspace];
    prodspace = 0;
    %prodspace = GradedSpace.new(charges2, [1 1 1 1 1 1 1 1 1], false);
    %{
    if kwargs.half_filling
        vspace = [pspace pspace*pspace];
    else
        error('TBA: now only half filling')
    end
    %}
    U1_values = @(x) [mod(x-1, 2*Q) + 1 mod(x-1, 2*Q) + 1 - 2*Q];
    if mod(P, 2) == 1
        length_vspaces = 2*Q;
    else
        length_vspaces = Q;
    end

    vspace = repmat(trivspace, 1, length_vspaces);
    for ind = 1:length_vspaces
        vspace(ind) = GradedSpace.new(get_asymm(U1_values(2*Q - ind*P), U1_values(Q - ind*P)), repmat(kwargs.D, 1, 6), false);
    end

    function list = get_asymm(list1, list2)
        ind_i = list1(1);
        ind_j = list1(2);
        ind_k = list2(1);
        ind_l = list2(2);
        list = [ProductCharge(U1(ind_i), SU2(1), fZ2(0)), ...
        ProductCharge(U1(ind_i), SU2(3), fZ2(0)), ...
        ProductCharge(U1(ind_j), SU2(1), fZ2(0)), ...
        ProductCharge(U1(ind_j), SU2(3), fZ2(0)),...
        ProductCharge(U1(ind_k), SU2(2), fZ2(1)), ...
        ProductCharge(U1(ind_l), SU2(2), fZ2(1))];
    end
end