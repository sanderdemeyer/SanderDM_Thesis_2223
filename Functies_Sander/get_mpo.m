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

    cod = SumSpace([one(vspace) vspace one(vspace)], pspace);
    dom = SumSpace(pspace, [one(vspace), vspace, one(vspace)]);
    O = MpoTensor.zeros(cod, dom);
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
        else
            sz = 2*N+2;
            mpo = MpoTensor.zeros(sz, 1, sz, 1);
            mpo(1, 1, 1, 1) = MpoTensor(1);
            mpo(sz, 1, sz, 1) = MpoTensor(1);
            mpo(1, 1, N+2, 1) = MpoTensor(L);
            mpo(2*N+1, 1, sz, 1) = MpoTensor(R);
            for j = N+2 : 2*N
                mpo(j, 1, j+1, 1) = MpoTensor(1);
            end
            mpo_base = mpo;

            mpo(1, 1, 2, 1) = MpoTensor(L);
            mpo(1, 1, 3, 1) = MpoTensor(L);
            mpo_LL = mpo;
            mpo = mpo_base;

            mpo(1, 1, 2, 1) = MpoTensor(L);
            mpo(2, 1, sz, 1) = MpoTensor(R);
            for i = 3:N
                mpo(i, 1, i+1, 1) = MpoTensor(1);
            end
            mpo_LR = mpo;
            mpo = mpo_base;

            mpo(2, 1, sz, 1) = MpoTensor(R);
            mpo(N+1, 1, sz, 1) = MpoTensor(R);
            mpo_RR = mpo;
            
            mpo_full{2*N} = mpo_RR;
            mpo_full{N} = mpo_RR;
            mpo_full{1} = mpo_LL;
            mpo_full{N+1} = mpo_LL;
            for i = 2:N-1
                mpo_full{i} = mpo_LR;
                mpo_full{N+i} = mpo_LR;
            end
        end
    elseif strcmp('FullCylinder_New', type)
        if N == 0
            warning('Not implemented for FullCylinder, using Helix instead')
            mpo_full = get_mpo(H, 0, 'Helix');
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
            mpo_full = get_mpo(H, 0, 'Helix');
            return
        else
            sz = N+4;
            mpo = MpoTensor.zeros(sz, 1, sz, 1);
            mpo(1, 1, 1, 1) = MpoTensor(1);
            mpo(sz, 1, sz, 1) = MpoTensor(1);

            % Underlying code implements N-range interaction
            mpo(1, 1, 4, 1) = MpoTensor(L);
            mpo(N+3, 1, sz, 1) = MpoTensor(R);
            for j = 4 : N+2
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
            mpo(3, 1, 3, 1) = MpoTensor(1);
            mpo_B = mpo;
            mpo = mpo_base;

            % Underlying code implements 'D' type mpo
            % for which X: L, R, 1. Y: R. B: R.
            mpo(2, 1, sz, 1) = MpoTensor(R);
            mpo(3, 1, sz, 1) = MpoTensor(R);
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
    else
        error('Type not implemented, check the definitions')
    end
end