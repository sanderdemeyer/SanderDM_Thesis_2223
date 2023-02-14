function [pspace, vspace, trivspace, prodspace, fusion_trees] = get_spaces_Hubbard_without_U1(kwargs)
    arguments
        kwargs.D1 = 1
        kwargs.D2 = 1
    end

    % Implementation of the physical charge
    a = ProductCharge(U1(0), fZ2(0));
    b = ProductCharge(U1(1), fZ2(1));
    c = ProductCharge(U1(-1), fZ2(1));
    charges1 = [b c a];

    % Implementation of the charges of the productspace
    d = ProductCharge(U1(0), fZ2(0));
    e = ProductCharge(U1(1), fZ2(1));
    f = ProductCharge(U1(-1), fZ2(1));
    g = ProductCharge(U1(2), fZ2(0));
    h = ProductCharge(U1(-2), fZ2(0));
    charges2 = [d e f g h];

    trivcharge = ProductCharge(U1(0), fZ2(0));

    fusion_trees = FusionTree.new([2 2], charges1, false, charges1, false, charges1, false, charges1, false);
    pspace = GradedSpace.new(charges1, [1 1 2], false);
    trivspace = GradedSpace.new(trivcharge, 1, false);
    prodspace = GradedSpace.new(charges2, [1 1 1 1 1], false);

    vspace1 = GradedSpace.new(charges1, repmat(kwargs.D1, 1, 3), false);
    vspace2 = GradedSpace.new(charges2, repmat(kwargs.D2, 1, 5), false);
    vspace = [vspace1 vspace2];
end