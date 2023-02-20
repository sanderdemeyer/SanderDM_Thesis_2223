function mpo = get_mpo_1D(H2, H1, kwargs)
    arguments
        H2
        H1
        kwargs.len = 2;
        kwargs.convention = 'conventional' % first was used prior to 16 february 2023. Conventional is what is should be
    end
    if strcmp(kwargs.convention, 'first')
        [L, R, pspace, vspace] = tsvd_twosite(H2);
    elseif strcmp(kwargs.convention, 'conventional')
        [L, R, pspace, vspace] = tsvd_twosite_conventional(H2);
    else
        error('convention does not exist.')
    end

    cod = SumSpace([one(vspace) vspace one(vspace)], pspace);
    dom = SumSpace(pspace, [one(vspace), vspace, one(vspace)]);
    O = MpoTensor.zeros(cod, dom);

    O(3, 1, 3, 1) = 1;
    O(1, 1, 1, 1) = 1;
    O(1, 1, 2, 1) = L;
    O(2, 1, 3, 1) = R;

    O(1, 1, 3, 1) = H1;
    
    mpo = cell(1, kwargs.len);
    for i = 1:kwargs.len
        mpo{i} = O;
    end
end