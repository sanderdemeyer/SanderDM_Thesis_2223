function [prob_up, prob_down] = get_probabilities(gs_mps, pspace, trivspace, site_number, varargin)
    if site_number == 1
        AC = gs_mps.AC;
        tens_one_site = Tensor([pspace trivspace], [trivspace pspace]);
        
        tblocks_down = num2cell([1 0]);
        prob_one_site_down = fill_tensor(tens_one_site, tblocks_down);
        prob_one_site_down = tpermute(prob_one_site_down, [2 1 4 3], [2 2]);

        prob_down = contract(AC, [1 2 3], prob_one_site_down, [-1 4 -2 2], conj(AC), [1 4 3]);
        prob_down = prob_down.var.var;
        
        tblocks_up = num2cell([0 1]);
        prob_one_site_up = fill_tensor(tens_one_site, tblocks_up);
        prob_one_site_up = tpermute(prob_one_site_up, [2 1 4 3], [2 2]);

        prob_up = contract(AC, [1 2 3], prob_one_site_up, [-1 4 -2 2], conj(AC), [1 4 3]);
        prob_up = prob_up.var.var;
        return
    elseif site_number == 2
        AC1 = gs_mps.AC(1);
        AC2 = gs_mps.AC(2);
        tens_one_site = Tensor([pspace trivspace], [trivspace pspace]);
        
        tblocks_down = num2cell([1 0]);
        prob_one_site_down = fill_tensor(tens_one_site, tblocks_down);
        prob_one_site_down = tpermute(prob_one_site_down, [2 1 4 3], [2 2]);

        prob_down_1 = contract(AC1, [1 2 3], prob_one_site_down, [-1 4 -2 2], conj(AC1), [1 4 3]);
        prob_down_2 = contract(AC2, [1 2 3], prob_one_site_down, [-1 4 -2 2], conj(AC2), [1 4 3]);
        prob_down = ((prob_down_1+prob_down_2)/2);
        prob_down = prob_down.var.var;
        
        tblocks_up = num2cell([0 1]);
        prob_one_site_up = fill_tensor(tens_one_site, tblocks_up);
        prob_one_site_up = tpermute(prob_one_site_up, [2 1 4 3], [2 2]);

        prob_up_1 = contract(AC1, [1 2 3], prob_one_site_up, [-1 4 -2 2], conj(AC1), [1 4 3]);
        prob_up_2 = contract(AC2, [1 2 3], prob_one_site_up, [-1 4 -2 2], conj(AC2), [1 4 3]);
        prob_up = ((prob_up_1+prob_up_2)/2);
        prob_up = prob_up.var.var;
    end
end