function mpo = get_mpo_1D_ntn(Hopping_t1, Hopping_t2, H_onesite, kwargs)
    arguments
        Hopping_t1
        Hopping_t2
        H_onesite
        kwargs.len = 2;
    end

    [L1, R1, ~, ~] = tsvd_twosite(Hopping_t1);
    [L2, R2, pspace, vspace] = tsvd_twosite(Hopping_t2);

    cod = SumSpace([one(vspace) vspace one(vspace) vspace one(vspace)], pspace);
    dom = SumSpace(pspace, [one(vspace), vspace, vspace, one(vspace), one(vspace)]);
    O = MpoTensor.zeros(cod, dom);

    sz = 5;
    O(sz, 1, sz, 1) = 1;
    O(1, 1, 1, 1) = 1;
    O(1, 1, 2, 1) = L1;    
    O(2, 1, sz, 1) = R1;
    
    O(1, 1, 3, 1) = L2;
    O(3, 1, 4, 1) = 1;
    O(4, 1, sz, 1) = R2;

    O(1, 1, sz, 1) = H_onesite;
    
    mpo = cell(1, kwargs.len);
    for i = 1:kwargs.len
        mpo{i} = O;
    end
end