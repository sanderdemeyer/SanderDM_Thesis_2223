d = 3;
D = 5;
A = createMPS(D, d);

assert(isequal(A.dims, [D, d, D]), 'Generated MPS tensor has incorrect shape.')
assert(~isreal(A), 'MPS tensor should have complex values.')


A = normalizeMPS(A);
[l, r] = fixedPoints(A);

relTol = 1;
assert(isapprox(l, l', 'RelTol', relTol), 'left fixed point should be hermitian!')
assert(isapprox(r, r', 'RelTol', relTol), 'right fixed point should be hermitian!')

assert(isapprox(l, contract(A, [1, 2, -2], l, [3, 1], conj(A), [3, 2, -1]), 'RelTol', relTol), 'l should be a left fixed point!')
assert(isapprox(r, contract(A, [-1, 2, 1], r, [1, 3], conj(A), [-2, 2, 3]), 'RelTol', relTol), 'r should be a right fixed point!')
assert(abs(trace(l * r) - 1) < relTol, 'Left and right fixed points should be trace normalized!') % TODO: fix once traces are a thing

function A =  createMPS(D, d)
    % Returns a random complex MPS tensor.
    %
    % Parameters
    % ----------
    % D : int
    %     Bond dimension for MPS.
    % d : int
    %     Physical dimension for MPS.
    %
    % Returns
    % -------
    % A : :class:`Tensor` (D, d, D)
    %     MPS tensor with 3 legs,
    %     ordered left-bottom-right.
    A = Tensor.randnc([D, d, D]);
end

function E = createTransfermatrix(A)
    % Form the transfermatrix of an MPS.
    %
    % Parameters
    % ----------
    % A : :class:`Tensor` (D, d, D)
    %     MPS tensor with 3 legs,
    %     ordered left-bottom-right.

    % Returns
    % -------
    % E : :class:`Tensor` (D, D, D, D)
    %     Transfer matrix with 4 legs,
    %     ordered topLeft-bottomLeft-topRight-bottomRight.
    E = contract(A, [-1, 1, -3], conj(A), [-2, 1, -4], 'Rank', [2, 2]);
end

function Anew = normalizeMPS(A)
    % Normalize an MPS tensor.
    %
    % Parameters
    % ----------
    % A : :class:`Tensor` (D, d, D)
    %     MPS tensor with 3 legs,
    %     ordered left-bottom-right.
    %
    % Returns
    % -------
    % Anew : :class:`Tensor` (D, d, D)
    %     MPS tensor with 3 legs,
    %     ordered left-bottom-right.
    %
    % Complexity
    % ----------
    % O(D^6) algorithm,
    %     diagonalizing (D^2, D^2) matrix.

    % create transfer matrix and reshape to appropriate tensor map acting on right fixed point
    E = createTransfermatrix(A);
    E = permute(E, [1, 2, 4, 3]);

    % determine eigenvalue and rescale MPS tensor accordingly
    norm = eigsolve(E);
    Anew = A / sqrt(norm);
end

function l = leftFixedPoint(A)
    % Find left fixed point.
    %
    % Parameters
    % ----------
    % A : :class:`Tensor` (D, d, D)
    %     MPS tensor with 3 legs,
    %     ordered left-bottom-right.
    %
    % Returns
    % -------
    % l : :class:`Tensor` (D, D)
    %     left fixed point with 2 legs,
    %     ordered bottom-top.
    %
    % Complexity
    % ----------
    % O(D^6) algorithm,
    %     diagonalizing (D^2, D^2) matrix.

    % initialize transfer matrix and reshape to appropriate tensor map acting on left fixed point
    E = createTransfermatrix(A);
    E = permute(E, [4, 3, 1, 2]);

    % find fixed point
    [l, ~] = eigsolve(E);

    % interpret as matrix
    l = repartition(l, [1, 1]);

    % make left fixed point hermitian and positive semidefinite explicitly
    l = l * trace(l) / abs(trace(l)); % remove possible phase, actually forgot why I had to do this
    l = (l + l') / 2; % force hermitian
    l = l * sign(trace(l)); % force positive semidefinite
end

function r = rightFixedPoint(A)
    % Find right fixed point.
    %
    % Parameters
    % ----------
    % A : :class:`Tensor` (D, d, D)
    %     MPS tensor with 3 legs,
    %     ordered left-bottom-right.
    %
    % Returns
    % -------
    % r : :class:`Tensor` (D, D)
    %     left fixed point with 2 legs,
    %     ordered top-bottom.
    %
    % Complexity
    % ----------
    % O(D^6) algorithm,
    %     diagonalizing (D^2, D^2) matrix.

    % initialize transfer matrix and reshape to appropriate tensor map acting on right fixed point
    E = createTransfermatrix(A);
    E = permute(E, [1:E.rank(1), flip(1:E.rank(2)) + E.rank(1)]);

    % find fixed point
    [r, ~] = eigsolve(E);

    % interpret to matrix
    r = repartition(r, [1, 1]);

    % make right fixed point hermitian and positive semidefinite explicitly
    r = r * trace(r) / abs(trace(r)); % remove possible phase, actually forgot why I had to do this
    r = (r + r') / 2; % force hermitian
    r = r * sign(trace(r)); % force positive semidefinite
end

function [l, r] = fixedPoints(A)
    % Find normalized fixed points.
    %
    % Parameters
    % ----------
    % A : :class:`Tensor` (D, d, D)
    %     MPS tensor with 3 legs,
    %     ordered left-bottom-right.
    %
    % Returns
    % -------
    % l : :class:`Tensor` (D, D)
    %     left fixed point with 2 legs,
    %     ordered bottom-top.
    % r : :class:`Tensor` (D, D)
    %     right fixed point with 2 legs,
    %     ordered top-bottom.
    %
    % Complexity
    % ----------
    % O(D^6) algorithm,
    %     diagonalizing (D^2, D^2) matrix

    % find fixed points
    l = leftFixedPoint(A);
    r = rightFixedPoint(A);

    % calculate trace and normalize
    l = l / trace(l * r);
end
