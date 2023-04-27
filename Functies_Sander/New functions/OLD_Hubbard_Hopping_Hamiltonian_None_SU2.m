function cdc = OLD_Hubbard_Hopping_Hamiltonian_None_SU2(t, kwargs)
    arguments
        t
        kwargs.bitstring = none
        kwargs.convention = 'conventional'
    end
        
    % based on comparison with H_double([1, 3, 4, 2], [1, 3, 4, 2])

    assert(strcmp(kwargs.convention, 'conventional'), 'Only convention = conventional is implemented.');

    pspace = get_spaces_Hubbard_None_SU2();

    cdc = Tensor([pspace pspace pspace' pspace'], []);
    %{
    data = {reshape([0 0 0 0 0 11 13 15 0 19 21 23 0 27 29 31], 2, 2, 2, 2)
        reshape([33 35 37 39], 1, 1, 2, 2)
        reshape([41 43 45 47], 1, 2, 1, 2)
        reshape([49 51 53 55], 2, 1, 1, 2)
        reshape([57 59 61 63], 1, 2, 2)
        reshape([65 67 69 71], 2, 1, 2)
        reshape([0 75 77 79], 2, 2)
        81
        83};

    
    data = {reshape([0 0 0 0 0 0 0 0 0 19 0 0 0 0 0 31], 2, 2, 2, 2)
        reshape([0 -sqrt(2) -sqrt(2) 0], 1, 1, 2, 2)
        reshape([-sqrt(2) 0 0 sqrt(2)], 1, 2, 1, 2)
        reshape([0 0 0 0], 2, 1, 1, 2)
        reshape([0 0 0 0], 1, 2, 2)
        reshape([-sqrt(2) 0 0 sqrt(2)], 2, 1, 2)
        reshape([0 sqrt(2) sqrt(2) 0], 2, 2)
        81
        83};
    %}
    %{
    data = {reshape([0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0], 2, 2, 2, 2)
        reshape([0 -sqrt(2)*t -sqrt(2)*t 0], 1, 1, 2, 2)
        reshape([-sqrt(2)*t 0 0 sqrt(2)*t], 1, 2, 1, 2)
        reshape([0 0 0 0], 2, 1, 1, 2)
        reshape([0 0 0 0], 1, 2, 2)
        reshape([-sqrt(2)*t 0 0 sqrt(2)*t], 2, 1, 2)
        reshape([0 sqrt(2)*t sqrt(2)*t 0], 2, 2)
        0
        0};
    %}
    data = {reshape([0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0], 2, 2, 2, 2)
        reshape([0 -sqrt(2)*t*(-1)^kwargs.bitstring(1) -sqrt(2)*t*(-1)^kwargs.bitstring(2) 0], 1, 1, 2, 2)
        reshape([-sqrt(2)*t*(-1)^kwargs.bitstring(3) 0 0 sqrt(2)*t*(-1)^kwargs.bitstring(4)], 1, 2, 1, 2)
        reshape([0 0 0 0], 2, 1, 1, 2)
        reshape([0 0 0 0], 1, 2, 2)
        reshape([-sqrt(2)*t*(-1)^kwargs.bitstring(5) 0 0 sqrt(2)*t*(-1)^kwargs.bitstring(6)], 2, 1, 2)
        reshape([0 sqrt(2)*t*(-1)^kwargs.bitstring(7) sqrt(2)*t*(-1)^kwargs.bitstring(8) 0], 2, 2)
        0
        0};
    %cdc = fill_tensor(cdc, num2cell(t*sqrt(2)*[0 -1 0 0 1 0 0 -1 -1 0 -1 0 1 0 0 1 0 0 1 0]));
    cdc = fill_tensor(cdc, data);
    return
end