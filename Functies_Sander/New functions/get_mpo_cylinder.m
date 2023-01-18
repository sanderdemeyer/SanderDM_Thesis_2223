function mpo = get_mpo_cylinder(H2, H1, N, rungs)
    % double is boolean. If true, 2 rungs are included.
    [L, R, pspace, vspace] = tsvd_twosite(H2);
    
    sz = N+4;
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

        % Underlying code implements N-range interaction
        mpoA(1, 1, 4, 1) = L;
        mpoB(1, 1, 4, 1) = L;
        mpoD(1, 1, 4, 1) = L;
        mpoA(N+3, 1, sz, 1) = R;
        mpoB(N+3, 1, sz, 1) = R;
        mpoD(N+3, 1, sz, 1) = R;
        for j = 4 : N+2
            mpoA(j, 1, j+1, 1) = 1;
            mpoB(j, 1, j+1, 1) = 1;
            mpoD(j, 1, j+1, 1) = 1;
        end

        % Underlying code implements 'A' type mpo
        % for which X: L, R, 1. Y: L. B: L.
        mpoA(1, 1, 2, 1) = L;
        %mpoA(1, 1, 3, 1) = L;
        mpoA(1, 1, 3, 1) = tpermute(conj(R), [3 4 1 2]);
        
        % Underlying code implements 'B' type mpo
        % for which X: L, R, 1. Y: L, R. B: 1.
        mpoB(1, 1, 2, 1) = L;
        mpoB(2, 1, sz, 1) = R;
        mpoB(3, 1, 3, 1) = 1;

        % Underlying code implements 'D' type mpo
        % for which X: L, R, 1. Y: R. B: R.
        mpoD(2, 1, sz, 1) = R;
       % mpoD(3, 1, sz, 1) = twist(R,1);
        mpoD(3, 1, sz, 1) = twist(tpermute(conj(L), [3 4 1 2]),1);
                
        % Underlying code implements one-site interaction
        mpoA(1, 1, sz, 1) = H1;
        mpoB(1, 1, sz, 1) = H1;
        mpoD(1, 1, sz, 1) = H1;

        if rungs == 2
            mpo{2*N} = mpoD;
            mpo{N} = mpoD;
            mpo{1} = mpoA;
            mpo{N+1} = mpoA;
            for i = 2:N-1
                mpo{i} = mpoB;
                mpo{N+i} = mpoB;
            end
        elseif rungs == 1
            mpo{N} = mpoD;
            mpo{1} = mpoA;
            for i = 2:N-1
                mpo{i} = mpoB;
            end
        else
            error('TBA')
        end
end