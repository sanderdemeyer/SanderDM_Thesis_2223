function [pspace, vspace, trivspace, prodspace, fusion_trees] = get_spaces_Heisenberg_XXZ(D1, D2)
    charges = U1([1 -1]);

    fusion_trees = FusionTree.new([2 2], charges, false, charges, false, charges, false, charges, false);
    pspace = GradedSpace.new(charges, [1 1], false);
    trivspace = GradedSpace.new(U1(0), 1, false);
    prodspace = 0;
    vspace1 = GradedSpace.new(U1([-1 1]), [D1 D1], false);
    vspace2 = GradedSpace.new(U1([2, 0, -2]), [D2 D2 D2], false);
    vspace = [vspace1 vspace2];
end