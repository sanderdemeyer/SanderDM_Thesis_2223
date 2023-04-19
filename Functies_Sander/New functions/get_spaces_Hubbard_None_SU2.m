function [pspace, vspace, trivspace, prodspace, fusion_trees] = get_spaces_Hubbard_None_SU2(kwargs)
    arguments
        kwargs.D1 = 1
        kwargs.D2 = 1
    end

    % Implementation of the physical charge
    a = ProductCharge(SU2(1), fZ2(0));
    b = ProductCharge(SU2(2), fZ2(1));
    charges1 = [a b];

    trivcharge = ProductCharge(SU2(1), fZ2(0));

    c = ProductCharge(SU2(1), fZ2(1));
    d = ProductCharge(SU2(2), fZ2(0));

    e = ProductCharge(SU2(3), fZ2(0));
    f = ProductCharge(SU2(3), fZ2(1));

    g = ProductCharge(SU2(4), fZ2(0));
    h = ProductCharge(SU2(4), fZ2(1));

    fusion_trees = FusionTree.new([2 2], charges1, false, charges1, false, charges1, false, charges1, false);
    pspace = GradedSpace.new(charges1, [2 1], false);
    trivspace = GradedSpace.new(trivcharge, 1, false);
    prodspace = GradedSpace.new([a b c d], [1 1 1 1], false);

    vspace1 = GradedSpace.new(charges1, [2*kwargs.D1 kwargs.D1], false);
    vspace2 = GradedSpace.new(charges1, [2*kwargs.D2 kwargs.D2], false);
    vspace = [vspace1 vspace2];
end