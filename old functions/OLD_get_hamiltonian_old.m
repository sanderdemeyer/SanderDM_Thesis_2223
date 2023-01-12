function H = get_hamiltonian(type, varargin)
% The following types are implemented:
% * Hubbard: Full Hubbard hamiltonian
% * Hubbard_external_field: Full Hubbard hamiltonian in an external field
% * Hubbard_two_site: Returns the 2-site interaction (hopping terms) of the
% Hubbard hamiltonian.
% * Hubbard_one_site: Returns the 1-site interaction (self-interaction
% terms) of the Hubbard hamiltonian.
% * Hubbard_one_site_redefined: Does the same as above, but with a
% different definition of the Hubbard hamiltonian.
% * Different Heisenberg XXZ and XXX models
    if strcmp('Hubbard', type)
        % the arguments are
        % fusion_trees, pspace, t, U, mu
        fusion_trees = varargin{1};
        pspace = varargin{2};
        t = varargin{3};
        U = varargin{4};
        mu = varargin{5};

        tens = Tensor([pspace pspace], [pspace pspace]);
        t_trees = [[0 1 1 0] [1 1 2 0] [0 2 1 1] [1 2 2 1] [0 1 1 0] [1 1 2 0] [0 2 1 1] [1 2 2 1]];
        U_trees = [[2 0 2 0] [2 1 2 1] [2 2 2 2] [0 2 0 2] [1 2 1 2]];
        len = length(fusion_trees);
        H_vars = zeros(1, len);
        for i = 1:len
            charges_trees = GetMD5_helper(fusion_trees(i).charges);
            U1_charges = charges_trees{1};
            c = [U1_charges(1) U1_charges(2) U1_charges(4) U1_charges(5)];
            if ismember(c, t_trees)
                H_vars(i) = H_vars(i) - t;
            end
            if ismember(c, U_trees)
                H_vars(i) = H_vars(i) + U;
            end
            % Underlying to include the chemical potential
            H_vars(i) = H_vars(i) - mu*(c(1) + c(2));
        end            
        tblocks = num2cell(H_vars);
        H = fill_tensor(tens, tblocks);
        return

    elseif strcmp('Hubbard_external_field', type)
        % If type is Hubbard with 0 radius, the arguments are
        % fusion_trees, pspace, t, U, mu, h_field (not staggered)
        fusion_trees = varargin{1};
        pspace = varargin{2};
        t = varargin{3};
        U = varargin{4};
        mu = varargin{5};
        h_field = varargin{6};

        len = length(fusion_trees);
        tens = Tensor([pspace pspace], [pspace pspace]);
        %tens = Tensor.new([2 2], pspace, true, pspace, true, pspace, false, pspace, false);
        t_trees = [[0 1 1 0] [1 1 2 0] [0 2 1 1] [1 2 2 1] [0 1 1 0] [1 1 2 0] [0 2 1 1] [1 2 2 1]];
        U_trees = [[2 0 2 0] [2 1 2 1] [2 2 2 2] [0 2 0 2] [1 2 1 2]];
        h_trees_pos = [[1 0 1 0] [1 0 0 1] [0 1 1 0] [0 1 0 1]];
        h_trees_double_pos = [[1 1 1 1]];
        h_trees_neg = [[-1 0 -1 0] [-1 0 0 -1] [0 -1 -1 0] [0 -1 0 -1]];
        h_trees_double_neg = [[-1 -1 -1 -1]];
        H_vars = zeros(1, len);
        for i = 1:len
            charges_trees = GetMD5_helper(fusion_trees(i).charges);
            U1_1_charges = charges_trees{1};
            U1_2_charges = charges_trees{2};
            c1 = [U1_1_charges(1) U1_1_charges(2) U1_1_charges(4) U1_1_charges(5)];
            c2 = [U1_2_charges(1) U1_2_charges(2) U1_2_charges(4) U1_2_charges(5)];
            if ismember(c1, t_trees)
                H_vars(i) = H_vars(i) - t;
            end
            if ismember(c1, U_trees)
                H_vars(i) = H_vars(i) + U;
            end
            % Underlying to include the chemical potential
            H_vars(i) = H_vars(i) - mu*(c1(1) + c1(2));
            if ismember(c2, h_trees_pos)
                H_vars(i) = H_vars(i) - h_field;
            elseif ismember(c2, h_trees_neg)
                H_vars(i) = H_vars(i) + h_field;
            elseif ismember(c2, h_trees_double_pos)
                H_vars(i) = H_vars(i) - 2*h_field;
            elseif ismember(c2, h_trees_double_neg)
                H_vars(i) = H_vars(i) + 2*h_field;
            end
        end            
        tblocks = num2cell(H_vars);
        H = fill_tensor(tens, tblocks);
        return

    elseif strcmp('Hubbard_two_site', type)
        % The arguments are fusion_trees, pspace, t, mu
        fusion_trees = varargin{1};
        pspace = varargin{2};
        t = varargin{3};
        mu = varargin{4};
        half_filling = varargin{5};
        tens = Tensor([pspace pspace], [pspace pspace]);
        if half_filling
            C2 = 1;
            C1 = 0;
            C0 = -1;
            up = 1;
            down = -1;
        else
            warning('Be careful if youre not using half-filling, correct charges might not be defined yet');
        end
        %t_trees = [[C0 C1 C1 C0] [C1 C1 C2 C0] [0 2 1 1] [1 2 2 1] [0 1 1 0] [1 1 2 0] [0 2 1 1] [1 2 2 1]];
        %t_trees = [[C0 C1 C1 C0] [C1 C0 C0 C1] [C1 C1 C2 C0] [C0 C2 C1 C1] [C1 C2 C2 C1] [C1 C1 C2 C0] [C0 C2 C1 C1] [C1 C2 C2 C1]];
        t_U1_1_trees = {[C0 C1 C1 C0] [C0 C1 C1 C0] [C0 C2 C1 C1] [C0 C2 C1 C1] [C1 C1 C2 C0] [C1 C1 C2 C0] [C1 C2 C2 C1] [C1 C2 C2 C1]};
        t_U1_2_trees = {[0 up up 0] [0 down down 0] [0 0 up down] [0 0 down up] [down up 0 0] [up down 0 0] [up 0 0 up] [down 0 0 down]};
        len = length(fusion_trees);
        H_vars = zeros(1, len);

        for i = 1:len
            charges_trees = GetMD5_helper(fusion_trees(i).charges);
            U1_1_charges = charges_trees{1};
            U1_2_charges = charges_trees{2};
            c1 = [U1_1_charges(1) U1_1_charges(2) U1_1_charges(4) U1_1_charges(5)];
            c2 = [U1_2_charges(1) U1_2_charges(2) U1_2_charges(4) U1_2_charges(5)];
            for k = 1:length(t_U1_1_trees)
                if all(t_U1_1_trees{k} == c1) && all(t_U1_2_trees{k} == c2)
                    H_vars(i) = H_vars(i) - t;
                end
            end
            %if ismember(c, t_trees)
            %    H_vars(i) = H_vars(i) - t;
            %end
            % Underlying to include the chemical potential
            %H_vars(i) = H_vars(i) - mu*(c(1) + c(2));
        end            
        tblocks = num2cell(H_vars);
        H = fill_tensor(tens, tblocks);
        return

    elseif strcmp('Hubbard_one_site', type)
        % arguments are pspace, trivspace, U
        
        pspace = varargin{1};
        trivspace = varargin{2};
        U = varargin{3};
        tens_one_site = Tensor([pspace trivspace], [trivspace pspace]);
        
        tblocks = num2cell([0 0 0 1]*U);
        
        H = fill_tensor(tens_one_site, tblocks);
        H = tpermute(H, [2 1 4 3], [2 2]);
        return

    elseif strcmp('Hubbard_one_site_redefined', type)
        % H_i = U*(n_up - 1/2)*(n_down - 1/2)
        % Other definition of the Hubbard model, with different 
        % pspace, trivspace, U
        pspace = varargin{1};
        trivspace = varargin{2};
        U = varargin{3};
        tens_one_site = Tensor([pspace trivspace], [trivspace pspace]);
        tens_one_site = tpermute(tens_one_site, [2 1 4 3], [2 2]);
        tblocks = num2cell([1/4 -1/4 -1/4 1/4]*U);
        H = fill_tensor(tens_one_site, tblocks);
        return

    elseif strcmp('XXZ', type)
        % arguments are
        % delta, pspace
        delta = varargin{1};
        pspace = varargin{2};
        tens = Tensor([pspace pspace], [pspace pspace]);

        H_vars = [delta, 2, -delta, -delta, 2, delta]/4;
        tblocks = num2cell(H_vars);
        H = fill_tensor(tens, tblocks);
        H = tpermute(H, [3 4 1 2], [2 2]);
        return

    elseif strcmp('XXX_stag', type)
        % arguments are
        % staggered pspace, trivspace, h_field, order
        pspace = varargin{1};
        trivspace = varargin{2};
        stagh = varargin{3};
        order = varargin{4};
        tens = Tensor([pspace pspace], [pspace pspace]);

        H_vars = [1, 2, -1, -1, 2, 1]/4;
        tblocks = num2cell(H_vars);
        H = fill_tensor(tens, tblocks);

        tens_one_site_l = Tensor(pspace, [trivspace pspace]);
        tens_one_site_r = Tensor(pspace, [trivspace' pspace]);

        tens_l = tpermute(tens_one_site_l, [2 3 1], [1 2]); %T_l has triv leg on the right
        tens_r = tpermute(tens_one_site_r, [2 3 1], [2 1]);
        tblocks_A = num2cell(stagh*[-1 1]);
        one_blocks = num2cell([1 1]);
        sigma_l = fill_tensor(tens_l, tblocks_A);
        sigma_r = fill_tensor(tens_r, tblocks_A);
        one_l = fill_tensor(tens_l, one_blocks);
        one_r = fill_tensor(tens_r, one_blocks);
        M_s_1 = contract(sigma_l, [-1 1 -3], one_r, [-2 1 -4]);
        M_s_2 = contract(one_l, [-1 1 -3], sigma_r, [-2 1 -4]);
        if order == 'AB'
            M_s = M_s_1 - M_s_2;
        elseif order == 'BA'
            M_s = -M_s_1 + M_s_2;
        else
            error('order not implemented, check for spelling errors')
        end
        M_s = tpermute(M_s, [1 2 3 4], [2 2]);
        H = tpermute(H, [3 4 1 2], [2 2]);
        H = plus(H, M_s);
        return
        
    elseif strcmp('one_site_XXZ', type)
        % arguments are
        % pspace, trivspace, h_field, schear
        pspace = varargin{1};
        trivspace = varargin{2};
        stag_h_field = varargin{3};
        h_field = varargin{4};

        tens_one_site = Tensor([pspace trivspace], [trivspace pspace]);
        
        tblocks_A = num2cell([-1 -1]*h_field + [-1 1]*stag_h_field);
        tblocks_B = num2cell([-1 -1]*h_field + [1 -1]*stag_h_field);
        
        H_one_site_A = fill_tensor(tens_one_site, tblocks_A);
        H_one_site_B = fill_tensor(tens_one_site, tblocks_B);
        H = {H_one_site_A H_one_site_B};
        return

    elseif strcmp('XXX_Ferromagnet', type)
    
    else
        error('type %s not implemented.\n', type)
    end        
end
