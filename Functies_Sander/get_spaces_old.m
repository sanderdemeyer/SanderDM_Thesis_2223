function [pspace, vspace, trivspace, fusion_trees] = get_spaces_old(type, varargin)
    if type == 'Hubbard'
        % For type = Hubbard, the arguments are SU2, P, Q, 
        % degeneracies
        % SU2 is a boolean and indicates whether spin up and down are
        % equivalent.
        % filling = P/Q.
        % half-filling: P = 1. Q = 1.
        % degeneracies are the degeneracies of EACH charge in the virtual
        % spaces, if the MPS is two-site, these are 2 integers.
        % Options: 
        % 0 e- present
        % 1 e- present with spin down
        % 1 e- present with spin down
        % 2 e- present
        % In total there are thus 4 options:
        % 0, up, down, en doubly occupated.

        SU2 = varargin{1};
        P = varargin{2};
        Q = varargin{3};
        D1 = varargin{4};
        D2 = varargin{5};
        if SU2  
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
            pspace = GradedSpace.new(charges1, [1 1 1 1], false);
            trivspace = GradedSpace.new(trivcharge, 1, false);
            
            vspace1 = GradedSpace.new(charges1, [D1 D1 D1 D1], false);
            vspace2 = GradedSpace.new(charges2, [D2 D2 D2 D2 D2 D2 D2 D2 D2], false);
            vspace = [vspace1 vspace2];
            %{
            if filling == 0
                D = varargin{3};

                a = ProductCharge(U1(0), U1(0), fZ2(0));
                b = ProductCharge(U1(1), U1(1), fZ2(1));
                c = ProductCharge(U1(2), U1(0), fZ2(0));
                charges = [a b c];

                fusion_trees = FusionTree.new([2 2], charges, false, charges, false, charges, false, charges, false);
                pspace = GradedSpace.new(charges, [1 1 1], false);
                trivspace = GradedSpace.new(a, 1, false);
                vspace = GradedSpace.new(charges, [D D D], false);
                return
            elseif filling == 1/2
                D1 = varargin{3};
                D2 = varargin{4};

                a = ProductCharge(U1(-1), U1(0), fZ2(0));
                b = ProductCharge(U1(0), U1(1), fZ2(1));
                c = ProductCharge(U1(1), U1(0), fZ2(0));
                charges1 = [a b c];

                e = ProductCharge(U1(-2), U1(0), fZ2(0));
                f = ProductCharge(U1(-1), U1(-1), fZ2(1));
                g = ProductCharge(U1(-1), U1(1), fZ2(1));
                h = ProductCharge(U1(0), U1(-2), fZ2(0));
                i = ProductCharge(U1(0), U1(0), fZ2(0));
                j = ProductCharge(U1(0), U1(2), fZ2(0));
                k = ProductCharge(U1(1), U1(-1), fZ2(1));
                l = ProductCharge(U1(1), U1(1), fZ2(1));
                m = ProductCharge(U1(2), U1(0), fZ2(0));
                charges2 = [e f g h i j k l m];

                fusion_trees = FusionTree.new([2 2], charges1, false, charges1, false, charges1, false, charges1, false);
                pspace = GradedSpace.new(charges1, [1 1 1 1], false);
                trivspace = GradedSpace.new(a, 1, false);
                
                vspace1 = GradedSpace.new(charges1, [D1 D1 D1 D1], false);
                vspace2 = GradedSpace.new(charges2, [D2 D2 D2 D2 D2 D2 D2 D2 D2], false);
                vspace = [vspace1 vspace2];
                return

            else
                error('filling %s is not implemented', filling);
            end
            %}
        else % if SU2 = false
            %{
            if filling == 0
                D = varargin{3};

                a = ProductCharge(U1(0), U1(0), fZ2(0));
                b = ProductCharge(U1(1), U1(1), fZ2(1));
                c = ProductCharge(U1(1), U1(-1), fZ2(1));
                d = ProductCharge(U1(2), U1(0), fZ2(0));
                charges = [a b c d];

                fusion_trees = FusionTree.new([2 2], charges, false, charges, false, charges, false, charges, false);
                pspace = GradedSpace.new(charges, [1 1 1 1], false);
                trivspace = GradedSpace.new(a, 1, false);
                vspace = GradedSpace.new(charges, [D D D D], false);
                return
            elseif filling == 1/2
                D1 = varargin{3};
                D2 = varargin{4};

                a = ProductCharge(U1(-1), U1(0), fZ2(0));
                b = ProductCharge(U1(0), U1(1), fZ2(1));
                c = ProductCharge(U1(0), U1(-1), fZ2(1));
                d = ProductCharge(U1(1), U1(0), fZ2(0));
                charges1 = [a b c d];

                e = ProductCharge(U1(-2), U1(0), fZ2(0));
                f = ProductCharge(U1(-1), U1(-1), fZ2(1));
                g = ProductCharge(U1(-1), U1(1), fZ2(1));
                h = ProductCharge(U1(0), U1(-2), fZ2(0));
                i = ProductCharge(U1(0), U1(0), fZ2(0));
                j = ProductCharge(U1(0), U1(2), fZ2(0));
                k = ProductCharge(U1(1), U1(-1), fZ2(1));
                l = ProductCharge(U1(1), U1(1), fZ2(1));
                m = ProductCharge(U1(2), U1(0), fZ2(0));
                charges2 = [e f g h i j k l m];

                fusion_trees = FusionTree.new([2 2], charges1, false, charges1, false, charges1, false, charges1, false);
                pspace = GradedSpace.new(charges1, [1 1 1 1], false);
                trivspace = GradedSpace.new(i, 1, false);
                
                vspace1 = GradedSpace.new(charges1, [D1 D1 D1 D1], false);
                vspace2 = GradedSpace.new(charges2, [D2 D2 D2 D2 D2 D2 D2 D2 D2], false);
                vspace = [vspace1 vspace2];
                return

            else
            %}
            % arbitrary filling f= P/Q

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
            
            vspace1 = GradedSpace.new(charges1, [D1 D1 D1 D1], false);
            vspace2 = GradedSpace.new(charges2, [D2 D2 D2 D2 D2 D2 D2 D2 D2], false);
            vspace = [vspace1 vspace2];
            return
        end
    elseif type == 'Heisenberg XXZ'
        error('Heisenberg not implemented');
    end
end