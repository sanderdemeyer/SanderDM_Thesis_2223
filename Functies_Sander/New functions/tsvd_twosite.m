function [L, R, pspace, vspace] = tsvd_twosite(H)
    % Split a two-site Hamiltonian in a contraction of L and R using svd
    
    warning('This holds for convention = first. This is not recommended.');
    [U, S, V] = tsvd(H, [1 3], [2 4], 'TruncBelow', 1e-12);
    U = insert_onespace(U, 4, false);
    L = contract(U, [-1 -3 1 -4], S, [1 -2], 'Rank', [2 2]);
    R = tpermute(insert_onespace(V, 3, true), [2 3 4 1]);
    
    % switch tensor order:
    L = tpermute(L, [4 3 2 1], [2 2]);
    R = tpermute(R, [4 3 2 1], [2 2]);

    pspace = L.domain(1);
    vspace = L.domain(2);
end