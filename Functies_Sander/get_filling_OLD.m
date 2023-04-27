function number = get_filling_OLD(gs_mps, pspace, trivspace, site_number, SU2, varargin)
    % SU2 is a boolean and indicates whether spin up and spin down are the
    % same
    if site_number == 1
        if SU2
            tens_one_site = Tensor([pspace trivspace], [trivspace pspace]);        
            tblocks = num2cell([0 1 2]);
            number_operator = fill_tensor(tens_one_site, tblocks);
            AC = gs_mps.AC;
            number = contract(AC, [1 2 3], number_operator, [2 -1 4 -2], conj(AC), [1 4 3]);
            number = number.var.var;
            return
        else
            tens_one_site = Tensor([pspace trivspace], [trivspace pspace]);        
            tblocks = num2cell([0 1 1 2]);
            number_operator = fill_tensor(tens_one_site, tblocks);
            AC = gs_mps.AC;
            number = contract(AC, [1 2 3], number_operator, [2 -1 4 -2], conj(AC), [1 4 3]);
            number = number.var.var;
            return
        end
    else
        if SU2
            tens_one_site = Tensor([pspace trivspace], [trivspace pspace]);        
            tblocks = num2cell([0 1 2]);
            number_operator = fill_tensor(tens_one_site, tblocks);
            AC1 = gs_mps.AC(1);
            AC2 = gs_mps.AC(2);
            number1 = contract(AC1, [1 2 3], number_operator, [2 -1 4 -2], conj(AC1), [1 4 3]);
            number2 = contract(AC2, [1 2 3], number_operator, [2 -1 4 -2], conj(AC2), [1 4 3]);
            number = (number1 + number2)/2;
            number = number.var.var;
            return
        else
            tens_one_site = Tensor([pspace trivspace], [trivspace pspace]);        
            tblocks = num2cell([0 1 1 2]);
            number_operator = fill_tensor(tens_one_site, tblocks);
            AC1 = gs_mps.AC(1);
            AC2 = gs_mps.AC(2);
            number1 = contract(AC1, [1 2 3], number_operator, [2 -1 4 -2], conj(AC1), [1 4 3]);
            number2 = contract(AC2, [1 2 3], number_operator, [2 -1 4 -2], conj(AC2), [1 4 3]);
            number = (number1 + number2)/2;
            number = number.var.var;
            return
        end
    end
end