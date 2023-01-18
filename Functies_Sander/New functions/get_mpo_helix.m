function mpo = get_mpo_helix(H2, H1, N, kwargs)
    arguments
        H2
        H1
        N
        kwargs.len = 2;
    end

    [L, R, pspace, vspace] = tsvd_twosite(H2);

    assert(floor(N) == N, 'Radius (N) should be an integer.')
    sz = N+3;

    cod_spaces = repmat(one(vspace), 1, sz);
    dom_spaces = repmat(one(vspace), 1, sz);
    cod_spaces(N+1) = vspace;
    cod_spaces(N+2) = vspace;
    dom_spaces(2) = vspace;
    dom_spaces(sz-1) = vspace;
    O = MpoTensor.zeros(SumSpace(cod_spaces, pspace), SumSpace(pspace, dom_spaces));

    O(sz,1,sz,1) = 1;
    O(1,1,1,1) = 1;
    O(1, 1, 2, 1) = L;
    O(1, 1, sz-1, 1) = L;
    
    for i = 2:N
        O(i,1,i+1,1) = 1;
    end
    O(N+1,1, sz,1) = R;
    O(N+2,1, sz,1) = R;

    O(1, 1, sz, 1) = H1;
    
    mpo = cell(1, kwargs.len);
    for i = 1:kwargs.len
        mpo{i} = O;
    end
end