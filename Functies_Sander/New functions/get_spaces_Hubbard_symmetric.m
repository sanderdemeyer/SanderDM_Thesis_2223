function [pspace, vspace, trivspace, prodspace, fusion_trees] = get_spaces_Hubbard_symmetric(P, Q, kwargs)
    arguments
        P
        Q
        kwargs.SU2_symm = false
        kwargs.D1 = 1
        kwargs.D2 = 1
    end

    if gcd(P,Q) ~= 1
        warning('P = %d and Q = %d should be relative prime, working with P = %d and Q = %d', P, Q, P/gcd(P,Q), Q/gcd(P,Q))
    end
        
% For type = 'Hubbard', which implements the fermionic hubbard
% model, the arguments are SU2, P, Q, degeneracies
% SU2 is a boolean and indicates whether spin up and down are
% equivalent. SU2 = true implements SU(2) symmetry, SU2 = false
% implements U(1) symmetry.
% filling = P/Q, with P and Q integers. Maximum is 2.
% half-filling: P = 1. Q = 1.
% degeneracies are the degeneracies of EACH charge in the virtual
% spaces, if the MPS is two-site, these are 2 integers.
% The implementation is based on the following Hilbert space: 
% 0 e- present
% 1 e- present with spin down
% 1 e- present with spin down
% 2 e- present
% In total, the hilbert space is 4-dimensional with elements:
% 0, up, down, en doubly occupated.
if kwargs.SU2_symm  
    a = ProductCharge(U1(-P), SU2(0), fZ2(0));
    b = ProductCharge(U1(Q-P), SU2(2), fZ2(1));
    c = ProductCharge(U1(2*Q-P), SU2(0), fZ2(0));
    charges1 = [a b c];

    d = ProductCharge(U1(-2*P), SU2(0), fZ2(0));
    e = ProductCharge(U1(Q-2*P), SU2(2), fZ2(1));
    f = ProductCharge(U1(2*Q-2*P), SU2(0), fZ2(0));
    g = ProductCharge(U1(3*Q-2*P), SU2(2), fZ2(1));
    h = ProductCharge(U1(4*Q-2*P), SU2(0), fZ2(0));
    charges2 = [d e f g h];

    trivcharge = ProductCharge(U1(0), SU2(0), fZ2(0));

    fusion_trees = FusionTree.new([2 2], charges1, false, charges1, false, charges1, false, charges1, false);
    pspace = GradedSpace.new(charges1, [1 1 1], false);
    trivspace = GradedSpace.new(trivcharge, 1, false);
    prodspace = GradedSpace.new(charges2, [1 1 1 1 1], false);

    vspace1 = GradedSpace.new(charges1, repmat(kwargs.D1, 1, 3), false);
    vspace2 = GradedSpace.new(charges2, repmat(kwargs.D2, 1, 5), false);
    vspace = [vspace1 vspace2];
    return
else
    a = ProductCharge(U1(-P), U1(0), fZ2(0));
    b = ProductCharge(U1(Q-P), U1(1), fZ2(1));
    c = ProductCharge(U1(Q-P), U1(-1), fZ2(1));
    d = ProductCharge(U1(2*Q-P), U1(0), fZ2(0));
    charges1 = [a b c d];

    e = ProductCharge(U1(-2*P), U1(0), fZ2(0));
    f = ProductCharge(U1(Q-2*P), U1(1), fZ2(1));
    g = ProductCharge(U1(Q-2*P), U1(-1), fZ2(1));
    h = ProductCharge(U1(2*Q-2*P), U1(0), fZ2(0));

    i = ProductCharge(U1(2*Q-2*P), U1(2), fZ2(0));
    j = ProductCharge(U1(2*Q-2*P), U1(0), fZ2(0)); % same as h
    k = ProductCharge(U1(3*Q-2*P), U1(1), fZ2(1));
    l = ProductCharge(U1(2*Q-2*P), U1(-2), fZ2(0));
    m = ProductCharge(U1(3*Q-2*P), U1(-1), fZ2(1));
    n = ProductCharge(U1(4*Q-2*P), U1(0), fZ2(0));
    charges2 = [e f g h i k l m n];

    trivcharge = ProductCharge(U1(0), U1(0), fZ2(0));

    fusion_trees = FusionTree.new([2 2], charges1, false, charges1, false, charges1, false, charges1, false);
    pspace = GradedSpace.new(charges1, [1 1 1 1], false);
    trivspace = GradedSpace.new(trivcharge, 1, false);
    prodspace = GradedSpace.new(charges2, [1 1 1 1 1 1 1 1 1], false);

    vspace1 = GradedSpace.new(charges1, repmat(kwargs.D1, 1, 4), false);
    vspace2 = GradedSpace.new(charges2, repmat(kwargs.D2, 1, 9), false);
    vspace = [vspace1 vspace2];
    return
end
