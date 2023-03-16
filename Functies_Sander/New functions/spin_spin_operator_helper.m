[tc, charges] = tensorblocks(O);
c = GetMD5_helper(charges.charges{2});
len = 36;
data = zeros(1, 36);
for i = 1:len
    U1_values = c(i,:);
    U1_v = [U1_values(1:2) U1_values(4:5)];
    if U1_v == [-1 -1 -1 -1]
        data(i) = 1/4;
    elseif U1_v == [1 -1 -1 1]
        data(i) = 1/2;
    elseif U1_v == [-1 1 -1 1]
        data(i) = -1/4;
    elseif U1_v == [1 -1 1 -1]
        data(i) = -1/4;
    elseif U1_v == [-1 1 1 -1]
        data(i) = 1/2;
    elseif U1_v == [1 1 1 1]
        data(i) = 1/4;
    end
end

data2 = [0 0 0 0 0 0 0 0 0 1/4 0 0 0 0 0 1/2 -1/4 0 0 -1/4 1/2 0 0 0 0 0 1/4 0 0 0 0 0 0 0 0 0];

norm(data-data2)