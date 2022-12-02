function spin_correlation = get_spin_correlation(gs_mps, pspace, trivspace, fusion_trees, site_number, filling, fermion, max_dist)
    % site_number is 1 or 2
    % fermion indicates whether there should be a twist on the right
    % virtual leg, as is the case with fermions.
    % max_dist is the maximal distance that will be calculated
    % filling is the filling (0 for 0 zero filling, 1/2 for half-filling,...)
    % Figure 1 of 'Ground-state properties...'
    % Calculates <n_{0,up}n_{i,down}> -  <n_{0,up}><n_{i,down}>
    % Does this by 
    if site_number > 2
        error('Only one-site and two-site implemented')
    end
    AC1 = gs_mps.AC(1);
    
    tens = Tensor([pspace pspace], [pspace pspace]);
    [tcharges, tblocks] = tensorblocks(tens);
    O_vars = zeros(1, 36);
    for i = 1:36
        charges_trees = GetMD5_helper(fusion_trees(i).charges);
        U1_1_charges = charges_trees{1};
        U1_2_charges = charges_trees{2};
        c1 = [U1_1_charges(1) U1_1_charges(2) U1_1_charges(4) U1_1_charges(5)];
        c2 = [U1_2_charges(1) U1_2_charges(2) U1_2_charges(4) U1_2_charges(5)];
        if filling == 0
            disp('Disp in get_number_correlation, erase!')
            if all([c1(1) c1(2)] == [c1(3) c1(4)]) && all([c2(1) c2(2)] == [c2(3) c2(4)])
                if c1(1) == 2 || c2(1) == 1
                    s1 = 1;
                else
                    s1 = 0;
                end
                if c1(2) == 2 || c2(2) == -1
                    s2 = 1;
                else
                    s2 = 0;
                end
                O_vars(i) = n1_up*n2_down;
            end
        elseif filling == 1/2
            %if all([c1(1) c1(2)] == [c1(3) c1(4)]) && all([c2(1) c2(2)] == [c2(3) c2(4)])
            if c1(1) == 1
                if c2(1) == 1
                    s1 = 1;
                else
                    s1 = 1; %should be -1
                end
            else
                s1 = 0;
            end

            if c1(1) == 1
                if c2(1) == 1
                    s2 = 1;
                else
                    s2 = 1; %should be -1
                end
            else
                s2 = 0;
            end
            
            O_vars(i) = s1*s2;
            %end
        else
            error('only zero and half filling implemented');
        end
    end            
    tblocks = num2cell(O_vars);
    O = fill_tensor(tens, tblocks);
    O = tpermute(O, [3 4 2 1], [2 2]);
    corr_list1 = correlation_function(O, gs_mps, 2, fermion, max_dist);
    %O = tpermute(O, [3 4 2 1], [2 2]);
    %O = tpermute(O, [1 2 4 3], [2 2]);
%     disp(O);
%     disp('This is O');
%     disp(reshape(double(O), 16, 16));
    %corr_list2 = correlation_function(O, gs_mps, 2, fermion, max_dist);
    %disp(corr_list1);
    %disp(corr_list2);
    %plot(1:max_dist, corr_list2 - corr_list1);
    tens_one_site = Tensor([pspace trivspace], [trivspace pspace]);
    
    N_up = Tensor([pspace' trivspace'], pspace');
    N_down = Tensor([pspace' trivspace], pspace');
    N_down = fill_tensor(N_down, num2cell([0 1 1 0])); % should be [0 -1 1 0]
    N_up = fill_tensor(N_up, num2cell([0 1 1 0])); % should be [0 1 -1 0]
    O = contract(N_up, [-1 1 -3], N_down, [-2 1 -4]);
    corr_list = correlation_function(O, gs_mps, 2, fermion, max_dist);
    tblocks_up = num2cell([0 0 1 1]);
    tblocks_down = num2cell([0 1 0 1]);
    
    N_up_one_site = fill_tensor(tens_one_site, tblocks_up);
    N_down_one_site = fill_tensor(tens_one_site, tblocks_down);
    
    %N_up_one_site = tpermute(N_up_one_site, [3 4 1 2], [2 2]);
    %N_down_one_site = tpermute(N_down_one_site, [3 4 1 2], [2 2]);

    N_up_one_site = tpermute(N_up_one_site, [2 1 4 3], [2 2]);
    N_down_one_site = tpermute(N_down_one_site, [2 1 4 3], [2 2]);

    %N_up = contract(AC1, [1 2 3], N_up_one_site, [2 -1 4 -2], conj(AC1), [1 4 3]);
    %N_down = contract(AC1, [1 2 3], N_down_one_site, [2 -1 4 -2], conj(AC1), [1 4 3]);
    if fermion
        N_up = contract(AC1, [1 2 3], N_up_one_site, [-1 4 -2 2], conj(twist(AC1,3)), [1 4 3]);
        N_down = contract(AC1, [1 2 3], N_down_one_site, [-1 4 -2 2], conj(twist(AC1,3)), [1 4 3]);
    else
        N_up = contract(AC1, [1 2 3], N_up_one_site, [-1 4 -2 2], conj(AC1), [1 4 3]);
        N_down = contract(AC1, [1 2 3], N_down_one_site, [-1 4 -2 2], conj(AC1), [1 4 3]);
    end
    N_up = N_up.var.var;
    N_down = N_down.var.var;
    
    if site_number == 2
        AC2 = gs_mps.AC(2);
        %N_up2 = contract(AC2, [1 2 3], N_up_one_site, [2 -1 4 -2], conj(AC2), [1 4 3]);
        %N_down2 = contract(AC2, [1 2 3], N_down_one_site, [2 -1 4 -2], conj(AC2), [1 4 3]);
        if fermion
            N_up2 = contract(AC2, [1 2 3], N_up_one_site, [-1 4 -2 2], conj(twist(AC2,3)), [1 4 3]);
            N_down2 = contract(AC2, [1 2 3], N_down_one_site, [-1 4 -2 2], conj(twist(AC2,3)), [1 4 3]);
        else
            N_up2 = contract(AC2, [1 2 3], N_up_one_site, [-1 4 -2 2], conj(AC2), [1 4 3]);
            N_down2 = contract(AC2, [1 2 3], N_down_one_site, [-1 4 -2 2], conj(AC2), [1 4 3]);
        end
        N_up2 = N_up2.var.var;
        N_down2 = N_down2.var.var;
        N_up = (N_up + N_up2)/2;
        N_down = (N_down + N_down2)/2;
    end
    number_correlation = corr_list1 - (N_up)*(N_down);
    plot(1:max_dist , number_correlation, 'Color', [0.8500 0.3250 0.0980]);
    hold on
    scatter(1:max_dist, number_correlation, "blue", "x");
    hold on
    plot(1:max_dist, zeros(1,max_dist), "black");
    hold on
    xlabel('Distance between site i and site j');
    ylabel('$Correlation: \langle n_{i,up}n_{j,down} \rangle -  \langle n_{i,up} \rangle \langle n_{j,down} \rangle$', 'interpreter', 'latex');
    title('Correlation function of the 1D Hubbard Model, not redefined')
    hold off
end