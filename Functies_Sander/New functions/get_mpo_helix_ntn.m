function mpo = get_mpo_helix_ntn(Hopping_t1, Hopping_t2, H_onesite, N, kwargs)
    arguments
        Hopping_t1
        Hopping_t2
        H_onesite
        N
        kwargs.len = 2;
    end

    [L1, R1, ~, ~] = tsvd_twosite(Hopping_t1);
    [L2, R2, pspace, vspace] = tsvd_twosite(Hopping_t2);

    assert(floor(N) == N, 'Radius (N) should be an integer.')
    sz = 3*N+5;

    % dom_spaces based on the 'L'
    % cod_spaces based on the 'R'
    cod_spaces = repmat(one(vspace), 1, sz);
    dom_spaces = repmat(one(vspace), 1, sz);
    cod_spaces(2) = vspace;
    cod_spaces(4) = vspace;
    cod_spaces(N+4) = vspace;
    cod_spaces(3*N+4) = vspace;
    dom_spaces(2) = vspace;
    dom_spaces(3) = vspace;
    dom_spaces(5) = vspace;
    dom_spaces(N+5) = vspace;
    O = MpoTensor.zeros(SumSpace(cod_spaces, pspace), SumSpace(pspace, dom_spaces));

    O(sz,1,sz,1) = 1;
    O(1,1,1,1) = 1;
    O(1, 1, 2, 1) = L1;
    O(1, 1, 3, 1) = L2;
    O(1, 1, 5, 1) = L1;
    O(1, 1, N+5, 1) = L2;
    
    O(3, 1, 4, 1) = 1;
    for i = 5:N+3
        O(i,1,i+1,1) = 1;
    end
    for i = N+5:3*N+3
        O(i,1,i+1,1) = 1;
    end
    O(2, 1, sz, 1) = R1;
    O(4, 1, sz, 1) = R2;
    O(N+4, 1, sz, 1) = R1;
    O(3*N+4, 1, sz, 1) = R2;

    O(1, 1, sz, 1) = H_onesite;
    
    mpo = cell(1, kwargs.len);
    for i = 1:kwargs.len
        mpo{i} = O;
    end
end