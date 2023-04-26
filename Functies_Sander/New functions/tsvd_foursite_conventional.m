function [L1, L2, L3, L4] = tsvd_foursite_conventional(H)
    % Split a two-site Hamiltonian in a contraction of L and R using svd
    
    [U, S, R] = tsvd(H, [1 8], [2 3 4 5 6 7], 'TruncBelow', 1e-12);
    L1 = contract(U, [-1 -3 1], S, [1 -2]);

    [U, S, R2] = tsvd(R, [1 2 7], [3 4 5 6], 'TruncBelow', 1e-12);
    L2 = contract(U, [-1 -2 -4 1], S, [1 -3]);

    [U, S, L4] = tsvd(R2, [1 2 5], [3 4], 'TruncBelow', 1e-12);
    L3 = contract(U, [-1 -2 -4 1], S, [1 -3]);
end