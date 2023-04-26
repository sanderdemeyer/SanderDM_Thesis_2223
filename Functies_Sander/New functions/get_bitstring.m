function bitstring = get_bitstring(a, base, kwargs)
    arguments
        a
        base
        kwargs.length = 0
        kwargs.arrows = false
    end
    if kwargs.length == 0
        length = floor(log(a)/log(base)) + 1;
    else
        length = kwargs.length;
    end
    bitstring = zeros(1, length);
    for i = 0:length-1
        bitstring(i+1) = floor(mod(a, base^(i+1))/base^i);
    end
    bitstring = flip(bitstring);
    if kwargs.arrows == true
        bitstring_new = cell(1, length);
        for i = 1:length
            if bitstring(i) == 0
                bitstring_new{i} = '_0_';
            elseif bitstring(i) == 1
                bitstring_new{i} = '_down_';
            elseif bitstring(i) == 2
                bitstring_new{i} = '_up_';
            elseif bitstring(i) == 3
                bitstring_new{i} = '_double_';
            end
        end
        bitstring =  cell2mat(bitstring_new);
    end
end