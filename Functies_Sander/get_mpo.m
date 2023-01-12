function O = get_mpo(H, N, type)
% Get the InfJMpo of a given two-site operator
% N is the radius of the cylinder.
% If this is 0, the model is just 1D
% If this is > 0, it should be an integer
    [U, S, V] = tsvd(H, [1 3], [2 4], 'TruncBelow', 1e-12);
    U = insert_onespace(U, 4, false);
    L = contract(U, [-1 -3 1 -4], S, [1 -2], 'Rank', [2 2]);
    R = tpermute(insert_onespace(V, 3, true), [2 3 4 1]);
    
    % switch tensor order:
    L = tpermute(L, [4 3 2 1], [2 2]);
    R = tpermute(R, [4 3 2 1], [2 2]);

    pspace = L.domain(1);
    vspace = L.domain(2);
    %{
    cod = SumSpace([one(vspace) vspace one(vspace)], pspace);
    dom = SumSpace(pspace, [one(vspace), vspace, one(vspace)]);
    L.domain = dom;
    L.codomain = cod;
    R.domain = dom;
    R.codomain = cod;
    %}
    if strcmp('Helix', type)
        if N == 0
            pspace = L.domain(1);
            vspace = L.domain(2);
%            trivspace = L.codomain(1);
            cod = SumSpace([one(vspace) vspace one(vspace)], pspace);
            dom = SumSpace(pspace, [one(vspace), vspace, one(vspace)]);
            O = MpoTensor.zeros(cod, dom);
            %O = MpoTensor.zeros([pspace' pspace'], [pspace' pspace']);
            O(3, 1, 3, 1) = 1;
            O(1, 1, 1, 1) = 1;
            O(1, 1, 2, 1) = L;
            O(2, 1, 3, 1) = R;
            return
        else
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
            return
        end
    elseif strcmp('Helix_3D', type)
        P = N(1);
        Q = N(2);
        if P == 0 || Q == 0
            error('One of the dimensions is zero.')
        else
            assert(floor(P) == P, 'Radius (P) should be an integer.')
            assert(floor(Q) == Q, 'Radius (Q) should be an integer.')
            sz = P+Q+3;
            mpo = MpoTensor.zeros(sz, 1, sz, 1);
            mpo(sz,1,sz,1) = MpoTensor(1);
            mpo(1,1,1,1) = MpoTensor(1);

            mpo(1, 1, 2, 1) = MpoTensor(L);
            mpo(2, 1, sz, 1) = MpoTensor(R);

            mpo(1, 1, 3, 1) = MpoTensor(L);
            for i = 3:P+1
                mpo(i,1,i+1,1) = MpoTensor(1);
            end
            mpo(P+2, 1, sz, 1) = MpoTensor(R);

            mpo(1, 1, P+3, 1) = MpoTensor(L);
            for j = P+3:P+Q+1
                mpo(j,1,j+1,1) = MpoTensor(1);
            end
            mpo(P+Q+2, 1, sz, 1) = MpoTensor(R);
            mpo_full = mpo;
        end
    elseif strcmp('FullCylinder_inefficient', type)
        % This has the same (I think) functionality as FullCylinder
        if N == 0
            warning('Not implemented for FullCylinder, using Helix instead')
            mpo_full = get_mpo(H, 0, 'Helix');
            return
        %{
        elseif N == 2
            sz = 5;
            mpo = MpoTensor.zeros(sz, 1, sz, 1);
            mpo(1, 1, 1, 1) = MpoTensor(1);
            mpo(sz, 1, sz, 1) = MpoTensor(1);
            mpo(1, 1, 2, 1) = MpoTensor(L);
            mpo(2, 1, sz, 1) = MpoTensor(R);
            mpo(1, 1, 3, 1) = MpoTensor(L);
            mpo(3, 1, 4, 1) = MpoTensor(1);
            mpo(4, 1, sz, 1) = MpoTensor(R);
            mpo_full = {mpo mpo mpo mpo};
            return
        %}
        else
            sz = 2*N+2;
            %{
            cod_spacesA = repmat(one(vspace), 1, sz);
            dom_spacesA = repmat(one(vspace), 1, sz);
            cod_spacesB = repmat(one(vspace), 1, sz);
            dom_spacesB = repmat(one(vspace), 1, sz);
            cod_spacesD = repmat(one(vspace), 1, sz);
            dom_spacesD = repmat(one(vspace), 1, sz);
            
            dom_spacesA(N+2) = vspace;
            dom_spacesB(N+2) = vspace;
            dom_spacesD(N+2) = vspace;
            cod_spacesA(2*N+1) = vspace;
            cod_spacesB(2*N+1) = vspace;
            cod_spacesD(2*N+1) = vspace;
            
            dom_spacesA(2) = vspace;
            dom_spacesA(3) = vspace;

            dom_spacesB(2) = vspace;
            cod_spacesB(2) = vspace;

            cod_spacesD(2) =  vspace;
            cod_spacesD(N+1) =  vspace;
            %}
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
            mpoA(1, 1, N+2, 1) = L;
            mpoB(1, 1, N+2, 1) = L;
            mpoD(1, 1, N+2, 1) = L;
            mpoA(2*N+1, 1, sz, 1) = R;
            mpoB(2*N+1, 1, sz, 1) = R;
            mpoD(2*N+1, 1, sz, 1) = R;
            for j = N+2 : 2*N
                mpoA(j, 1, j+1, 1) = 1;
                mpoB(j, 1, j+1, 1) = 1;
                mpoD(j, 1, j+1, 1) = 1;
            end

            % Underlying code implements 'A' type mpo
            mpoA(1, 1, 2, 1) = L;
            mpoA(1, 1, 3, 1) = L;

            % Underlying code implements 'B' type mpo
            mpoB(1, 1, 2, 1) = L;
            mpoB(2, 1, sz, 1) = R;
            for i = 3:N
                mpoB(i, 1, i+1, 1) = 1;
            end

            % Underlying code implements 'D' type mpo
            mpoD(2, 1, sz, 1) = R;
            mpoD(N+1, 1, sz, 1) = R;
            
            mpo_full{2*N} = mpoD;
            mpo_full{N} = mpoD;
            mpo_full{1} = mpoA;
            mpo_full{N+1} = mpoA;
            for i = 2:N-1
                mpo_full{i} = mpoB;
                mpo_full{N+i} = mpoB;
            end
            O = mpo_full;
        end
    elseif strcmp('FullCylinder_New', type)
        if N == 0
            warning('Not implemented for FullCylinder, using Helix instead')
            O = get_mpo(H, 0, 'Helix');
            return
        else
            sz = 2*N+2;
            mpo = MpoTensor.zeros(sz, 1, sz, 1);
            mpo(1, 1, 1, 1) = MpoTensor(1);
            mpo(sz, 1, sz, 1) = MpoTensor(1);

            % Underlying code implements N-range interaction
            mpo(1, 1, N+2, 1) = MpoTensor(L);
            mpo(2*N+1, 1, sz, 1) = MpoTensor(R);
            for j = N+2 : 2*N
                mpo(j, 1, j+1, 1) = MpoTensor(1);
            end
            mpo_base = mpo;

            % Underlying code implements 'A' type mpo
            % for which X: L, R, 1. Y: L. B: L.
            mpo(1, 1, 2, 1) = MpoTensor(L);
            mpo(1, 1, 3, 1) = MpoTensor(L);
            mpo_A = mpo;
            mpo = mpo_base;
            
            % Underlying code implements 'B' type mpo
            % for which X: L, R, 1. Y: L, R. B: 1.
            mpo(1, 1, 2, 1) = MpoTensor(L);
            mpo(2, 1, sz, 1) = MpoTensor(R);
            for i = 3:N
                mpo(i, 1, i+1, 1) = MpoTensor(1);
            end
            mpo_B = mpo;
            mpo = mpo_base;

            % Underlying code implements 'D' type mpo
            % for which X: L, R, 1. Y: R. B: R.
            mpo(2, 1, sz, 1) = MpoTensor(R);
            mpo(N+1, 1, sz, 1) = MpoTensor(R);
            mpo_D = mpo;
                        
            mpo_full{2*N} = mpo_D;
            mpo_full{N} = mpo_D;
            mpo_full{1} = mpo_A;
            mpo_full{N+1} = mpo_A;
            for i = 2:N-1
                mpo_full{i} = mpo_B;
                mpo_full{N+i} = mpo_B;
            end
        end
    elseif strcmp('FullCylinder', type)
        if N == 0
            warning('Not implemented for FullCylinder, using Helix instead')
            mpo = get_mpo(H, 0, 'Helix');
            O = {mpo mpo};
            return
        else
            sz = N+4;
            %{
            cod_spacesA = repmat(one(vspace), 1, sz);
            dom_spacesA = repmat(one(vspace), 1, sz);
            cod_spacesB = repmat(one(vspace), 1, sz);
            dom_spacesB = repmat(one(vspace), 1, sz);
            cod_spacesD = repmat(one(vspace), 1, sz);
            dom_spacesD = repmat(one(vspace), 1, sz);

            dom_spacesA(4) = vspace;
            dom_spacesB(4) = vspace;
            dom_spacesD(4) = vspace;
            cod_spacesA(N+3) = vspace;
            cod_spacesB(N+3) = vspace;
            cod_spacesD(N+3) = vspace;

            dom_spacesA(2) = vspace;
            dom_spacesA(3) = vspace;

            dom_spacesB(2) = vspace;
            cod_spacesB(2) = vspace;

            cod_spacesD(2) =  vspace;
            cod_spacesD(3) =  vspace;
            %}
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
            %{
            cod_spacesA(2) = one(vspace);
            cod_spacesA(3) = one(vspace);
            dom_spacesD(2) = one(vspace);
            dom_spacesD(3) = one(vspace);
            %}
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
            mpoA(1, 1, 3, 1) = L;
            
            % Underlying code implements 'B' type mpo
            % for which X: L, R, 1. Y: L, R. B: 1.
            mpoB(1, 1, 2, 1) = L;
            mpoB(2, 1, sz, 1) = R;
            mpoB(3, 1, 3, 1) = 1;

            % Underlying code implements 'D' type mpo
            % for which X: L, R, 1. Y: R. B: R.
            mpoD(2, 1, sz, 1) = R;
            mpoD(3, 1, sz, 1) = R;
                        
            mpo_full{2*N} = mpoD;
            mpo_full{N} = mpoD;
            mpo_full{1} = mpoA;
            mpo_full{N+1} = mpoA;
            for i = 2:N-1
                mpo_full{i} = mpoB;
                mpo_full{N+i} = mpoB;
            end
            O = mpo_full;
        end
    else
        error('Type not implemented, check the definitions')
    end
end