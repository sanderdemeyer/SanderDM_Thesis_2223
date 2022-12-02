function magn = get_magnetisation(type, gs_mps, pspace, trivspace, site_number, fermion, staggered)
    % fermion is a boolean and indicates whether one works with fermions
    % staggered is a boolean and indicates whether or not the staggered
    % magnetisation will be calculated.
    if site_number == 1
        error('This will probably be wrong, look in function get_magnetisation')
        tens_one_site = Tensor([pspace trivspace], [trivspace pspace]);  
        if strcmp('Hubbard', type)
            tblocks = num2cell([0 -1 1 0])
        elseif strcmp('XXZ', type)
            tblocks = num2cell([-1 1]);
        else
            error('Type not implemented')
        end
            M = fill_tensor(tens_one_site, tblocks);
        AC = gs_mps.AC;
        magn = contract(AC, [1 2 3], M, [2 -1 4 -2], conj(AC), [1 4 3]);
        magn = magn.var.var;
        return
    elseif site_number == 1.5 % 2, but old
        tens_one_site = Tensor([pspace trivspace], [trivspace pspace]);
        tens_one_site = tpermute(tens_one_site, [2 1 4 3], [2 2]);
        if strcmp('Hubbard', type)
            tblocks_A = num2cell([0 -1 1 0]);
        elseif strcmp('XXZ', type)
            tblocks_A = num2cell([-1 1]);
        else
            error('Type not implemented')
        end
        M_A = fill_tensor(tens_one_site, tblocks_A);
        AC1 = gs_mps.AC(1);
        AC2 = gs_mps.AC(2);
        if fermion
            magnetisation1 = contract(AC1, [1 2 3], M_A, [-1 4 -2 2], conj(twist(AC1,3)), [1 4 3]);
        else
            magnetisation1 = contract(AC1, [1 2 3], M_A, [-1 4 -2 2], conj(AC1), [1 4 3]);
        end
        %magnetisation1 = contract(AC1, [1 2 3], M_A, [2 -1 4 -2], conj(AC1), [1 4 3]);

        if staggered
            if strcmp('Hubbard', type)
                tblocks_B = num2cell([0 1 -1 0]);
            elseif strcmp('XXZ', type)
                tblocks_B = num2cell([1 -1]);
            else
                error('Type not implemented')
            end
            M_B = fill_tensor(tens_one_site, tblocks_B);
        else
            M_B = M_A;
        end
        if fermion
            magnetisation2 = contract(AC2, [1 2 3], M_B, [-1 4 -2 2], conj(twist(AC2,3)), [1 4 3]);
        else
            magnetisation2 = contract(AC2, [1 2 3], M_B, [-1 4 -2 2], conj(AC2), [1 4 3]);
        end
        magn = (magnetisation1 + magnetisation2)/2;
        magn = magn.var.var;
        return
    elseif site_number == 2 % = 2 but new
        tens_one_site_l = Tensor(pspace, [trivspace pspace]);
        tens_one_site_r = Tensor(pspace, [trivspace' pspace]);
        tens_l = tpermute(tens_one_site_l, [2 3 1], [1 2]); %T_l has triv leg on the right
        tens_r = tpermute(tens_one_site_r, [2 3 1], [2 1]);
        if strcmp('Hubbard', type)
            tblocks_A = num2cell([0 -1 1 0]);
            one_blocks = num2cell([1 1 1 1]);
        elseif strcmp('XXZ', type)
            tblocks_A = num2cell([-1 1]);
            one_blocks = num2cell([1 1]);
        else
            error('Type not implemented')
        end
        sigma_l = fill_tensor(tens_l, tblocks_A);
        sigma_r = fill_tensor(tens_r, tblocks_A);
        one_l = fill_tensor(tens_l, one_blocks);
        one_r = fill_tensor(tens_r, one_blocks);
        M_s_1 = contract(sigma_l, [-1 1 -3], one_r, [-2 1 -4]);
        M_s_2 = contract(one_l, [-1 1 -3], sigma_r, [-2 1 -4]);
        M_s = M_s_1 - M_s_2;

        AL1 = gs_mps.AL(1);
        AL2 = gs_mps.AL(2);
        AC1 = gs_mps.AC(1);
        AC2 = gs_mps.AC(2);
        if fermion
            magnetisation1 = contract(AL1, [1 2 3], AC2, [3 4 8], conj(AL1), [1 5 6], conj(twist(AC2,3)), [6 7 8], M_s, [2 4 5 7]);
            magnetisation2 = contract(AL2, [1 2 3], AC1, [3 4 8], conj(AL2), [1 5 6], conj(twist(AC1,3)), [6 7 8], M_s, [2 4 5 7]);
        else
            magnetisation1 = contract(AL1, [1 2 3], AC2, [3 4 8], conj(AL1), [1 5 6], conj(AC2), [6 7 8], M_s, [2 4 5 7]);
            magnetisation2 = contract(AL2, [1 2 3], AC1, [3 4 8], conj(AL2), [1 5 6], conj(AC1), [6 7 8], M_s, [2 4 5 7]);
        end
        %magn = (magnetisation1 + magnetisation2)/2;
        magn = magnetisation1;
    elseif site_number == 2.1 % AB-structure, but general period of the unit cell
        w = period(gs_mps);
        assert(mod(w, 2) == 0, 'Period of the MPS should be even')

        tens_one_site = Tensor([pspace trivspace], [trivspace pspace]);
        tens_one_site = tpermute(tens_one_site, [2 1 4 3], [2 2]);
        if strcmp('Hubbard', type)
            tblocks_A = num2cell([0 -1 1 0]);
        elseif strcmp('XXZ', type)
            tblocks_A = num2cell([-1 1]);
        else
            error('Type not implemented')
        end
        M_A = fill_tensor(tens_one_site, tblocks_A);

        if staggered
            if strcmp('Hubbard', type)
                tblocks_B = num2cell([0 1 -1 0]);
            elseif strcmp('XXZ', type)
                tblocks_B = num2cell([1 -1]);
            else
                error('Type not implemented')
            end
            M_B = fill_tensor(tens_one_site, tblocks_B);
        else
            M_B = M_A;
        end

        magn = zeros(1,w);
        for i = 1:w
            if ((mod(i,2) == 0) - 1/2)*((i > w/2) - 1/2) == 1/4
                magn(i) = contract(gs_mps.AC(i), [1 2 3], M_A, [-1 4 -2 2], conj(twist(gs_mps.AC(i),3)), [1 4 3]).var.var;
            else
                magn(i) = contract(gs_mps.AC(i), [1 2 3], M_B, [-1 4 -2 2], conj(twist(gs_mps.AC(i),3)), [1 4 3]).var.var;
            end
        end
        %magnetisation1 = contract(AC1, [1 2 3], M_A, [2 -1 4 -2], conj(AC1), [1 4 3]);
        magn = mean(magn);
    end
end