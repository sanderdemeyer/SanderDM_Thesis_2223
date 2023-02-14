function Tens_final = get_open_indices(mps, indices, T)
    open = length(indices);
    assert(indices(1) == 1, 'First index must be 1')

    AC = mps.AC;
    AR = mps.AR;

    b = 1;
    w = period(mps);

    Tens = contract(AC(b), [1 -1 -3], conj(AC(b)), [1 -2 -4]);
    if indices(2) > 2
        Tens = contract(Tens, [-1 -2 1 2], T{loop(b,1,w),indices(2)-2}, [1 2 -3 -4]);
    end
    
    if open == 2
        Tens_final = contract(Tens, [-1 -2 1 2], AR(loop(b,indices(2)-1,w)), [1 -3 3], twist(conj(AR(loop(b,indices(2)-1,w))),3), [2 -4 3]);
    else
        Tens = contract(Tens, [-1 -2 1 2], AR(loop(b,indices(2)-1,w)), [1 -3 -5], conj(AR(loop(b,indices(2)-1,w))), [2 -4 -6]);
        if indices(3) - indices(2) > 1
            Tens = contract(Tens, [-1 -2 -3 -4 1 2], T{loop(b,indices(2),w),indices(3)-indices(2)-1}, [1 2 -5 -6]);
        end
        
        if open == 3
            Tens_final = contract(Tens, [-1 -2 -3 -4 1 2], AR(loop(b,indices(3)-1,w)), [1 -5 3], twist(conj(AR(loop(b,indices(3)-1,w))),3), [2 -6 3]);
        else
            Tens = contract(Tens, [-1 -2 -3 -4 1 2], AR(loop(b,indices(3)-1,w)), [1 -5 -7], conj(AR(loop(b,indices(3)-1,w))), [2 -6 -8]);
            if indices(4) - indices(3) > 1
                Tens = contract(Tens, [-1 -2 -3 -4 -5 -6 1 2], T{loop(b,indices(3),w), indices(4)-indices(3)-1}, [1 2 -7 -8]);
            end
            Tens_final = contract(Tens, [-1 -2 -3 -4 -5 -6 1 2], AR(loop(b,indices(4)-1,w)), [1 -7 3], twist(conj(AR(loop(b,indices(4)-1,w))),3), [2 -8 3]);
        end
    end
end