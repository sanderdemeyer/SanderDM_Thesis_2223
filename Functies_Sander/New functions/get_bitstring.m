function bitstring = get_bitstring(a)
    bitstring = zeros(1, 8);
    for i = 0:7
        bitstring(i+1) = floor(mod(a, 2^(i+1))/2^i);
    end
end