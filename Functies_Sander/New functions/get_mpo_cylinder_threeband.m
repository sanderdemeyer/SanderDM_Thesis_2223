function mpo = get_mpo_cylinder_threeband(H1_Cu, H1_O, H2_CuCu, H2_CuO, H2_OO, N, rungs, kwargs)
    arguments
        H1_Cu
        H1_O
        H2_CuCu
        H2_CuO
        H2_OO
        N
        rungs
        kwargs.convention = 'conventional' % first was used prior to 16 february 2023. Conventional is what is should be
    end
    if strcmp(kwargs.convention, 'first')
        [L_CuCu, R_CuCu, pspace, vspace] = tsvd_twosite(H2_CuCu);
        [L_CuO, R_CuO, ~, ~] = tsvd_twosite(H2_CuO);
        [L_OO, R_OO, ~, ~] = tsvd_twosite(H2_OO);
    elseif strcmp(kwargs.convention, 'conventional')
        [L_CuCu, R_CuCu, pspace, vspace] = tsvd_twosite_conventional(H2_CuCu);
        [L_CuO, R_CuO, ~, ~] = tsvd_twosite_conventional(H2_CuO);
        [L_OO, R_OO, ~, ~] = tsvd_twosite_conventional(H2_OO);
    else
        error('convention does not exist.')
    end
        sz = 4*N+13;

        cod_spaces = repmat(vspace, 1, sz);
        dom_spaces = repmat(vspace, 1, sz);

        cod_spaces(1) = one(vspace);
        cod_spaces(sz) = one(vspace);
        dom_spaces(1) = one(vspace);
        dom_spaces(sz) = one(vspace);

        %{
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
        %}

        mpoCuA = MpoTensor.zeros(SumSpace(cod_spaces, pspace), SumSpace(pspace, dom_spaces));
        mpoCuB = MpoTensor.zeros(SumSpace(cod_spaces, pspace), SumSpace(pspace, dom_spaces));
        mpoCuD = MpoTensor.zeros(SumSpace(cod_spaces, pspace), SumSpace(pspace, dom_spaces));
        mpoOxA = MpoTensor.zeros(SumSpace(cod_spaces, pspace), SumSpace(pspace, dom_spaces));
        mpoOx = MpoTensor.zeros(SumSpace(cod_spaces, pspace), SumSpace(pspace, dom_spaces));
        mpoOy = MpoTensor.zeros(SumSpace(cod_spaces, pspace), SumSpace(pspace, dom_spaces));
        mpoOyD = MpoTensor.zeros(SumSpace(cod_spaces, pspace), SumSpace(pspace, dom_spaces));



        mpoCuA(1, 1, 1, 1) = 1;
        mpoCuA(sz, 1, sz, 1) = 1;
        mpoCuB(1, 1, 1, 1) = 1;
        mpoCuB(sz, 1, sz, 1) = 1;
        mpoCuD(1, 1, 1, 1) = 1;
        mpoCuD(sz, 1, sz, 1) = 1;
        mpoOxA(1, 1, 1, 1) = 1;
        mpoOxA(sz, 1, sz, 1) = 1;
        mpoOx(1, 1, 1, 1) = 1;
        mpoOx(sz, 1, sz, 1) = 1;
        mpoOy(1, 1, 1, 1) = 1;
        mpoOy(sz, 1, sz, 1) = 1;
        mpoOyD(1, 1, 1, 1) = 1;
        mpoOyD(sz, 1, sz, 1) = 1;

        % implements Cu-Cu interaction in y direction
        mpoCuA(1, 1, 2, 1) = L_CuCu;

        mpoCuB(1, 1, 2, 1) = L_CuCu;
        mpoCuB(2, 1, sz, 1) = R_CuCu;

        mpoCuD(2, 1, sz, 1) = R_CuCu;

        mpoOxA(2, 1, 2, 1) = 1;
        mpoOx(2, 1, 2, 1) = 1;
        mpoOy(2, 1, 2, 1) = 1;
        mpoOyD(2, 1, 2, 1) = 1;

        % implements Cu-Cu interaction in x direction
        mpoCuA(1, 1, 3, 1) = L_CuCu;
        mpoCuB(1, 1, 3, 1) = L_CuCu;
        mpoCuD(1, 1, 3, 1) = L_CuCu;

        mpoOxA(3, 1, 3, 1) = 1;
        mpoOx(3, 1, 3, 1) = 1;
        mpoOy(3, 1, 3, 1) = 1;
        mpoOyD(3, 1, 3, 1) = 1;
        for i = 3:N+1
            mpoCuA(i, 1, i+1, 1) = 1;
            mpoCuB(i, 1, i+1, 1) = 1;
            mpoCuD(i, 1, i+1, 1) = 1;

            mpoOxA(i+1, 1, i+1, 1) = 1;
            mpoOx(i+1, 1, i+1, 1) = 1;
            mpoOy(i+1, 1, i+1, 1) = 1;
            mpoOyD(i+1, 1, i+1, 1) = 1;
        end

        mpoCuA(N+2, 1, sz, 1) = R_CuCu;
        mpoCuB(N+2, 1, sz, 1) = R_CuCu;
        mpoCuD(N+2, 1, sz, 1) = R_CuCu;

        % implements Cu-Ox-Oy interaction (y direction, within unit cell)
        % (positive)
        mpoCuA(1, 1, N+3, 1) = L_CuO;
        mpoCuB(1, 1, N+3, 1) = L_CuO;
        mpoCuD(1, 1, N+3, 1) = L_CuO;

        mpoOxA(N+3, 1, sz, 1) = R_CuO;
        mpoOxA(N+3, 1, N+3, 1) = 1;
        mpoOx(N+3, 1, sz, 1) = R_CuO;
        mpoOx(N+3, 1, N+3, 1) = 1;
        
        mpoOy(N+3, 1, sz, 1) = R_CuO;
        mpoOyD(N+3, 1, sz, 1) = R_CuO;

        % implements Ox-Cu to other rung (negative)
        mpoOxA(1, 1, N+4, 1) = -L_CuO;
        mpoOx(1, 1, N+4, 1) = -L_CuO;

        mpoOxA(N+4, 1, N+4, 1) = 1;
        mpoOx(N+4, 1, N+4, 1) = 1;
        mpoOy(N+4, 1, N+4, 1) = 1;
        mpoOyD(N+4, 1, N+4, 1) = 1;
        for i = N+4:2*N+2
            mpoCuA(i, 1, i+1, 1) = 1;
            mpoCuB(i, 1, i+1, 1) = 1;
            mpoCuD(i, 1, i+1, 1) = 1;

            mpoOxA(i+1, 1, i+1, 1) = 1;
            mpoOx(i+1, 1, i+1, 1) = 1;
            mpoOy(i+1, 1, i+1, 1) = 1;
            mpoOyD(i+1, 1, i+1, 1) = 1;
        end

        mpoCuA(2*N+3, 1, sz, 1) = R_CuO;
        mpoCuB(2*N+3, 1, sz, 1) = R_CuO;
        mpoCuD(2*N+3, 1, sz, 1) = R_CuO;

        % implements Oy-Cu interaction to other unit cell (negative)
        mpoOy(1, 1, 2*N+4, 1) = -L_CuO;
        mpoCuB(2*N+4, 1, sz, 1) = R_CuO;
        mpoCuD(2*N+4, 1, sz, 1) = R_CuO;

        % implements Oy-Cu interaction that loops around the cylinder
        % (negative)

        mpoCuA(1, 1, 2*N+5, 1) = -L_CuO;
        mpoOyD(2*N+5, 1, sz, 1) = R_CuO;

        mpoCuB(2*N+5, 1, 2*N+5, 1) = 1;
        mpoCuD(2*N+5, 1, 2*N+5, 1) = 1;
        mpoOxA(2*N+5, 1, 2*N+5, 1) = 1;
        mpoOx(2*N+5, 1, 2*N+5, 1) = 1;
        mpoOy(2*N+5, 1, 2*N+5, 1) = 1;

        % implements Cu-Cu interaction that loops around the cylinder

        mpoCuA(1, 1, 2*N+6, 1) = L_CuO;
        mpoCuD(2*N+6, 1, sz, 1) = R_CuO;

        mpoCuB(2*N+6, 1, 2*N+6, 1) = 1;
        mpoOxA(2*N+6, 1, 2*N+6, 1) = 1;
        mpoOx(2*N+6, 1, 2*N+6, 1) = 1;
        mpoOy(2*N+6, 1, 2*N+6, 1) = 1;
        mpoOyD(2*N+6, 1, 2*N+6, 1) = 1;

        % implements Ox-Oy interaction within same unit cell (negative)

        mpoOxA(1, 1, 2*N+7, 1) = -L_OO;
        mpoOx(1, 1, 2*N+7, 1) = -L_OO;
        mpoOy(2*N+7, 1, sz, 1) = R_OO;
        mpoOyD(2*N+7, 1, sz, 1) = R_OO;

        % implements Oy-Ox interaction to further unit cell (positive)

        mpoOy(1, 1, 2*N+8, 1) = L_OO;
        mpoOx(2*N+8, 1, sz, 1) = R_OO;

        mpoCuB(2*N+8, 1, 2*N+8, 1) = 1;
        mpoCuD(2*N+8, 1, 2*N+8, 1) = 1;

        % implements Ox-Oy interaction that loops around the cylinder
        % (positive)
        mpoOxA(1, 1, 2*N+9, 1) = L_OO;
        mpoOyD(2*N+9, 1, sz, 1) = R_OO;

        mpoCuB(2*N+9, 1, 2*N+9, 1) = 1;
        mpoCuD(2*N+9, 1, 2*N+9, 1) = 1;
        mpoOx(2*N+9, 1, 2*N+9, 1) = 1;
        mpoOy(2*N+9, 1, 2*N+9, 1) = 1;

        % implements Ox-Oy interaction to further rung on same height.
        % (positive)

        mpoOxA(1, 1, 2*N+10, 1) = L_OO;
        mpoOx(1, 1, 2*N+10, 1) = L_OO;

        mpoOxA(2*N+10, 1, 2*N+10, 1) = 1;
        mpoOx(2*N+10, 1, 2*N+10, 1) = 1;
        mpoCuA(2*N+10, 1, 2*N+10, 1) = 1;
        mpoCuB(2*N+10, 1, 2*N+10, 1) = 1;
        mpoCuD(2*N+10, 1, 2*N+10, 1) = 1;
        for i = 2*N+10:3*N+9
            mpoOy(i, 1, i+1, 1) = 1;
            mpoOyD(i, 1, i+1, 1) = 1;

            mpoOxA(i+1, 1, i+1, 1) = 1;
            mpoOx(i+1, 1, i+1, 1) = 1;
            mpoCuA(i+1, 1, i+1, 1) = 1;
            mpoCuB(i+1, 1, i+1, 1) = 1;
            mpoCuD(i+1, 1, i+1, 1) = 1;
        end

        mpoOy(3*N+10, 1, sz, 1) = R_OO;
        mpoOyD(3*N+10, 1, sz, 1) = R_OO;

        % implements Ox-Oy interaction to further rung 1 unit cell above.
        % (negative)

        mpoOx(1, 1, 3*N+11, 1) = -L_OO;

        mpoOxA(3*N+11, 1, 3*N+11, 1) = 1;
        mpoOx(3*N+11, 1, 3*N+11, 1) = 1;
        mpoCuA(3*N+11, 1, 3*N+11, 1) = 1;
        mpoCuB(3*N+11, 1, 3*N+11, 1) = 1;
        mpoCuD(3*N+11, 1, 3*N+11, 1) = 1;
        for i = 3*N+11:4*N+9
            mpoOy(i, 1, i+1, 1) = 1;
            mpoOyD(i, 1, i+1, 1) = 1;

            mpoOxA(i+1, 1, i+1, 1) = 1;
            mpoOx(i+1, 1, i+1, 1) = 1;
            mpoCuA(i+1, 1, i+1, 1) = 1;
            mpoCuB(i+1, 1, i+1, 1) = 1;
            mpoCuD(i+1, 1, i+1, 1) = 1;
        end

        mpoOy(4*N+10, 1, sz, 1) = R_OO;

        % implements the above Ox-Oy interaction to further rung that 
        % loops around the cylinder (negative)

        mpoOxA(1, 1, 4*N+11, 1) = -L_OO;
        mpoOyD(4*N+11, 1, 4*N+12, 1) = 1;
        mpoOyD(4*N+12, 1, sz, 1) = R_OO;

        mpoOxA(4*N+11, 1, 4*N+11, 1) = 1;
        mpoOx(4*N+11, 1, 4*N+11, 1) = 1;
        mpoOy(4*N+11, 1, 4*N+11, 1) = 1;
        mpoOyD(4*N+11, 1, 4*N+11, 1) = 1;
        mpoCuA(4*N+11, 1, 4*N+11, 1) = 1;
        mpoCuB(4*N+11, 1, 4*N+11, 1) = 1;
        mpoCuD(4*N+11, 1, 4*N+11, 1) = 1;

        mpoOxA(4*N+12, 1, 4*N+12, 1) = 1;
        mpoOx(4*N+12, 1, 4*N+12, 1) = 1;
        mpoOy(4*N+12, 1, 4*N+12, 1) = 1;
        mpoOyD(4*N+12, 1, 4*N+12, 1) = 1;
        mpoCuA(4*N+12, 1, 4*N+12, 1) = 1;
        mpoCuB(4*N+12, 1, 4*N+12, 1) = 1;
        mpoCuD(4*N+12, 1, 4*N+12, 1) = 1;

        % implements one-site interactions

        mpoCuA(1, 1, sz, 1) = H1_Cu;
        mpoCuB(1, 1, sz, 1) = H1_Cu;
        mpoCuD(1, 1, sz, 1) = H1_Cu;

        mpoOxA(1, 1, sz, 1) = H1_O;
        mpoOx(1, 1, sz, 1) = H1_O;
        mpoOy(1, 1, sz, 1) = H1_O;
        mpoOyD(1, 1, sz, 1) = H1_O;

        for j = rungs:-1:1
            mpo{j*N*3} = mpoOyD;
            mpo{j*N*3-1} = mpoOx;
            mpo{j*N*3-2} = mpoCuD;
            mpo{(j-1)*N*3+1} = mpoCuA;
            mpo{(j-1)*N*3+2} = mpoOxA;
            mpo{(j-1)*N*3+3} = mpoOy;
        
            for i = 1:N-2
                mpo{(j-1)*N*3+i*3+1} = mpoCuB;
                mpo{(j-1)*N*3+i*3+2} = mpoOx;
                mpo{(j-1)*N*3+i*3+3} = mpoOy;
            end
        end
end