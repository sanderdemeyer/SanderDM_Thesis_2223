function [pspace, vspace, trivspace, prodspace, fusion_trees] = get_spaces(type, varargin)
    if strcmp('Hubbard', type)
        % For type = 'Hubbard', which implements the fermionic hubbard
        % model, the arguments are SU2, P, Q, degeneracies
        % SU2 is a boolean and indicates whether spin up and down are
        % equivalent. SU2 = true implements SU(2) symmetry, SU2 = false
        % implements U(1) symmetry.
        % filling = P/Q, with P and Q integers.
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

        SU2_symm = varargin{1};
        P = varargin{2};
        Q = varargin{3};
        D1 = varargin{4};
        D2 = varargin{5};
        if SU2_symm  
            a = ProductCharge(U1(-P), SU2(0), fZ2(0));
            b = ProductCharge(U1(Q-P), SU2(1), fZ2(1));
            c = ProductCharge(U1(2*Q-P), SU2(0), fZ2(0));
            charges1 = [a b c];

            d = ProductCharge(U1(-2*P), SU2(0), fZ2(0));
            e = ProductCharge(U1(Q-2*P), SU2(1), fZ2(1));
            f = ProductCharge(U1(2*Q-2*P), SU2(0), fZ2(0));
            g = ProductCharge(U1(3*Q-2*P), SU2(1), fZ2(1));
            h = ProductCharge(U1(4*Q-2*P), SU2(0), fZ2(0));
            charges2 = [d e f g h];

            trivcharge = ProductCharge(U1(0), SU2(0), fZ2(0));

            fusion_trees = FusionTree.new([2 2], charges1, false, charges1, false, charges1, false, charges1, false);
            pspace = GradedSpace.new(charges1, [1 1 1], false);
            trivspace = GradedSpace.new(trivcharge, 1, false);
            prodspace = GradedSpace.new(charges2, [1 1 1 1 1]);

            vspace1 = GradedSpace.new(charges1, [D1 D1 D1], false);
            vspace2 = GradedSpace.new(charges2, [D2 D2 D2 D2 D2], false);
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

            vspace1 = GradedSpace.new(charges1, [D1 D1 D1 D1], false);
            vspace2 = GradedSpace.new(charges2, [D2 D2 D2 D2 D2 D2 D2 D2 D2], false);
            vspace = [vspace1 vspace2];
            return
        end
    elseif strcmp('Heisenberg XXZ', type)
        if nargin == 1
            D1 = [9 9];
            D2 = [6 6 6];
        else
            D1 = varargin{1};
            D2 = varargin{2};
        end
        charges = U1([1 -1]);

        fusion_trees = FusionTree.new([2 2], charges, false, charges, false, charges, false, charges, false);
        pspace = GradedSpace.new(charges, [1 1], false);
        trivspace = GradedSpace.new(U1(0), 1, false);

        vspace1 = GradedSpace.new(U1([-1 1]), D1, false);
        vspace2 = GradedSpace.new(U1([2, 0, -2]), D2, false);
        vspace = [vspace1 vspace2];
    else
        error('Type not implemented');
    end
end