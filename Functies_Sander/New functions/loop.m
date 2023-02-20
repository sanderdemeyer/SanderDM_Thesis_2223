function j = loop(b, step, total)
    j = mod(b+step-1, total) + 1;
end