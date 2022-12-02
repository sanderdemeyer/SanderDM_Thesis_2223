a = MpoTensor(1);
b = MpoTensor(2);
K = {a b a};
K{3} = b;
L = K;
L{2} = MpoTensor(8);

S{7} = MpoTensor(7);
N = 4;
for n = 1:2*N
    disp(n);
    if ((mod(n,2) == 0) - 1/2)*((n < N + 1) - 1/2) == 1/4
        disp('B');
    else
        disp('A');
    end
end
