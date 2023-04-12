function mpo = get_mpo_cylinder_oneband_nntn(H1, H2, H2_ntn, H2_nntn, N, rungs, kwargs)
    arguments
        H1
        H2
        H2_ntn
        H2_nntn
        N
        rungs
        kwargs.convention = 'conventional' % first was used prior to 16 february 2023. Conventional is what is should be
    end

    assert(N > 3, 'N should be bigger than 3. N = 3 is not yet implemented');

    if strcmp(kwargs.convention, 'first')
        [L1, R1, pspace, vspace] = tsvd_twosite(H2);
        [L2, R2, ~, ~] = tsvd_twosite(H2_ntn);
        [L3, R3, ~, ~] = tsvd_twosite(H2_nntn);
    elseif strcmp(kwargs.convention, 'conventional')
        [L1, R1, pspace, vspace] = tsvd_twosite_conventional(H2);
        [L2, R2, ~, ~] = tsvd_twosite_conventional(H2_ntn);
        [L3, R3, ~, ~] = tsvd_twosite_conventional(H2_nntn);
    else
        error('convention does not exist.')
    end
    
    sz = 6*N+8;
        cod_spacesA = repmat(vspace, 1, sz);
        dom_spacesA = repmat(vspace, 1, sz);
        cod_spacesB = repmat(vspace, 1, sz);
        dom_spacesB = repmat(vspace, 1, sz);
        cod_spacesD = repmat(vspace, 1, sz);
        dom_spacesD = repmat(vspace, 1, sz);

        cod_spacesA(1) = one(vspace);
        cod_spacesB(1) = one(vspace);
        cod_spacesD(1) = one(vspace);
        cod_spacesA(sz) = one(vspace);
        cod_spacesB(sz) = one(vspace);
        cod_spacesD(sz) = one(vspace);
        dom_spacesA(1) = one(vspace);
        dom_spacesB(1) = one(vspace);
        dom_spacesD(1) = one(vspace);
        dom_spacesA(sz) = one(vspace);
        dom_spacesB(sz) = one(vspace);
        dom_spacesD(sz) = one(vspace);

        mpoA = MpoTensor.zeros(SumSpace(cod_spacesA, pspace), SumSpace(pspace, dom_spacesA));
        mpoB = MpoTensor.zeros(SumSpace(cod_spacesB, pspace), SumSpace(pspace, dom_spacesB));
        mpoD = MpoTensor.zeros(SumSpace(cod_spacesD, pspace), SumSpace(pspace, dom_spacesD));

        mpoA(1, 1, 1, 1) = 1;
        mpoA(sz, 1, sz, 1) = 1;
        mpoB(1, 1, 1, 1) = 1;
        mpoB(sz, 1, sz, 1) = 1;
        mpoD(1, 1, 1, 1) = 1;
        mpoD(sz, 1, sz, 1) = 1;

        % implements 1-site interaction
        mpoA(1, 1, 2, 1) = L1;

        mpoB(1, 1, 2, 1) = L1;
        mpoB(2, 1, sz, 1) = R1;

        mpoD(2, 1, sz, 1) = R1;

        % implements n N-site interaction
        mpoA(1, 1, 3, 1) = L1;
        mpoB(1, 1, 3, 1) = L1;
        mpoD(1, 1, 3, 1) = L1;
        mpoA(N+2, 1, sz, 1) = R1;
        mpoB(N+2, 1, sz, 1) = R1;
        mpoD(N+2, 1, sz, 1) = R1;
        for j = 3 : N+1
            mpoA(j, 1, j+1, 1) = 1;
            mpoB(j, 1, j+1, 1) = 1;
            mpoD(j, 1, j+1, 1) = 1;
        end


        % implements n' N+1 interaction
        mpoA(1, 1, N+3, 1) = L2;
        mpoB(1, 1, N+3, 1) = L2;

        for j = N+3 : 2*N+2
            mpoA(j, 1, j+1, 1) = 1;
            mpoB(j, 1, j+1, 1) = 1;
            mpoD(j, 1, j+1, 1) = 1;
        end

        mpoB(2*N+3, 1, sz, 1) = R2;
        mpoD(2*N+3, 1, sz, 1) = R2;

        % implements n' N-1 interaction
        
        mpoB(1, 1, 2*N+4, 1) = L2;
        mpoD(1, 1, 2*N+4, 1) = L2;

        for j = 2*N+4 : 3*N+1
            mpoA(j, 1, j+1, 1) = 1;
            mpoB(j, 1, j+1, 1) = 1;
            mpoD(j, 1, j+1, 1) = 1;
        end

        mpoA(3*N+2, 1, sz, 1) = R2;
        mpoB(3*N+2, 1, sz, 1) = R2;

        % implements n N-1 interaction

        mpoA(1, 1, 3*N+3, 1) = L1;
        mpoB(3*N+3, 1, 3*N+3, 1) = 1;
        mpoD(3*N+3, 1, sz, 1) = R1;

        % implements n' 2*N+1 interaction

        mpoA(1, 1, 3*N+4, 1) = L2;
        mpoB(3*N+4, 1, 3*N+4, 1) = 1;
        mpoB(3*N+5, 1, 3*N+5, 1) = 1;
        mpoD(3*N+4, 1, 3*N+5, 1) = 1;
        mpoD(3*N+5, 1, sz, 1) = R2;
        
        % implements n'' 2*N interaction

        mpoA(1, 1, 3*N+6, 1) = L3;
        mpoB(1, 1, 3*N+6, 1) = L3;
        mpoD(1, 1, 3*N+6, 1) = L3;
        mpoA(5*N+5, 1, sz, 1) = R3;
        mpoB(5*N+5, 1, sz, 1) = R3;
        mpoD(5*N+5, 1, sz, 1) = R3;
        for j = 3*N+6 : 5*N+4
            mpoA(j, 1, j+1, 1) = 1;
            mpoB(j, 1, j+1, 1) = 1;
            mpoD(j, 1, j+1, 1) = 1;
        end

        % implements n'' interaction in same rung

        mpoA(1, 1, 5*N+6, 1) = L3;
        mpoB(1, 1, 5*N+6, 1) = L3;

        mpoA(5*N+6, 1, 5*N+7, 1) = 1;
        mpoB(5*N+6, 1, 5*N+7, 1) = 1;
        mpoD(5*N+6, 1, 5*N+7, 1) = 1;

        mpoB(5*N+7, 1, sz, 1) = R3;
        mpoD(5*N+7, 1, sz, 1) = R3;

        % implements n'' interaction to N-2. Starting at A.

        mpoA(1, 1, 5*N+8, 1) = L3;
        mpoB(6*N+6, 1, sz, 1) = R3;
        for j = 5*N+8:6*N+5
            mpoB(j, 1, j+1, 1) = 1;
        end

        % implements n'' interaction to N-2. Starting at B.

        mpoB(1, 1, 6*N+7, 1) = L3;
        mpoB(6*N+7, 1, 6*N+7, 1) = 1;
        mpoD(6*N+7, 1, sz, 1) = R3;

        % implements onsite interaction
        mpoA(1, 1, sz, 1) = H1;
        mpoB(1, 1, sz, 1) = H1;
        mpoD(1, 1, sz, 1) = H1;
        
        for j = rungs:-1:1
            mpo{j*N} = mpoD;
            mpo{(j-1)*N+1} = mpoA;
            for i = 2:N-1
                mpo{(j-1)*N+i} = mpoB;
            end
        end
end