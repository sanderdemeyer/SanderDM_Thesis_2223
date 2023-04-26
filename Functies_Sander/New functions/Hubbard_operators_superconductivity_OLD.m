function Hubbard_operators_superconductivity_OLD(P, Q)
    [pspace, ~, trivspace] = get_spaces_Hubbard_SU2(P, Q);

    charge = ProductCharge(U1(2*P), SU2(1), fZ2(0));
    e2_space = GradedSpace.new(charge, 1, false);

    %{
    c_up_c_down = Tensor([pspace pspace e2_space], [pspace pspace]);
    data = {0 -1 0 0 1 1 0 0 -1 0};
    c_up_c_down = fill_tensor(c_up_c_down, data);

    test = Tensor([e2_space e2_space e2_space e2_space], [e2_space e2_space e2_space e2_space]);
    %}

    Delta = Tensor([pspace pspace pspace pspace], [pspace pspace pspace pspace]);
    fusion_trees = fusiontrees(Delta);
    for i = [543 548 916 937 1168 1184 1499 1539]    %  i = 1:length(fusion_trees)
        charges_trees = GetMD5_helper(fusion_trees(i).charges);
        U1_charges = charges_trees{1};
        charges_cod = [U1_charges(1) U1_charges(2) U1_charges(4) U1_charges(6)];
        charges_dom = [U1_charges(8) U1_charges(10) U1_charges(12) U1_charges(13)];
        charges = [charges_cod, charges_dom];
        %{
        if charges(1) == charges(8) + 1 && ...
            charges(2) == charges(7) + 1 && ...
            charges(3) == charges(4) - 1 && ...
            charges(4) == charges(5) - 1
            SU2_charges = charges_trees{2};
            fuses_dom = [SU2_charges(9) SU2_charges(11)];
            fuses_cod = [SU2_charges(3) SU2_charges(5)];
            %disp(i);
            %disp(charges);
            %disp(fuses_cod);
            %disp(fuses_dom);
            %disp(fusion_trees(i));
        end
        %}
    end
    [tc, tb] = tensorblocks(Delta);
    tens = fill_tensor(Delta, tc);

    for i = 1:length(fusion_trees)
        tc{i} = 0;
    end
    tc{543} = 3;
    disp(tc{543});
    disp(tb(543));
    i = 543;
    charges_trees = GetMD5_helper(fusion_trees(i).charges);
    U1_charges = charges_trees{1};
    charges_cod = [U1_charges(1) U1_charges(2) U1_charges(4) U1_charges(6)];
    charges_dom = [U1_charges(8) U1_charges(10) U1_charges(12) U1_charges(13)];
    charges = [charges_cod, charges_dom];
    disp(charges);
    tc{548} = 0;
    tc{916} = 0;
    tc{937} = 0;
    tc{1168} = 0;
    tc{1184} = 0;
    tc{1499} = 0;
    tc{1539} = 0;
    tens = fill_tensor(Delta, tc);
    
    matr = reshape(double(tens), 256, 256);
    for i = 1:256
        for j = 1:256
            if matr(i, j) ~= 0
                fprintf('(i,j) = (%d,%d) \n', i, j);
                % i is output, j is input
                fprintf('bitstring is (%d %d %d %d , %d %d %d %d) \n', get_bitstring(i,4,'length', 4), get_bitstring(j,4,'length', 4)); 
                fprintf('bitstring is (%s %s %s %s | %s %s %s %s) \n', get_bitstring(i,4,'length', 4, 'arrows', true), get_bitstring(j,4,'length', 4, 'arrows', true)); 
                disp(matr(i,j)); 
            end
        end
    end
end