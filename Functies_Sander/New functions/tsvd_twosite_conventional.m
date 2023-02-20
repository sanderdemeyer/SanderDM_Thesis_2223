function [L, R, pspace, vspace] = tsvd_twosite_conventional(H)
% Split a two-site Hamiltonian in a contraction of L and R using svd
    [U, S, V] = tsvd(H, [4 1], [2 3], 'TruncBelow', 1e-12);
    U = insert_onespace(U, 2, false);
    L = contract(U, [-4 -1 -2 1], S, [1 -3], 'Rank', [2 2]);
    R = insert_onespace(V, 3, true);
    % make 2,2 tensor!
    pspace = L.domain(1);
    vspace = L.domain(2);
end