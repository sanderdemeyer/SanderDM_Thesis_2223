files = dir(('Data structures/XXZ_different_deltas_interpretation'));
l = length(files);
count = 1;
for i = 3:l
    disp('ok');
    file = files(i);
    name = file.name;
    len = length(name);
    if strcmp(name(len-8:len), 'final.mat')
        delta = str2double(name(28:len-40));
        fprintf('delta is %s \n', delta);
        if 1 == 1
            load(name);
            AC1 = gs_mps.AC(1);
            AC2 = gs_mps.AC(2);
            AL1 = gs_mps.AL(1);
            AL2 = gs_mps.AL(2);
            E1 = contract(AL1, [1 2 3], AC2, [3 4 8], H_A, [2 4 5 7], conj(AL1), [1 5 6], conj(AC2), [6 7 8]);
            E2 = contract(AL2, [1 2 3], AC1, [3 4 8], H_A, [2 4 5 7], conj(AL2), [1 5 6], conj(AC1), [6 7 8]);
            E = (E1+E2)/2;
            E_list{count} = E;
            delta_list{count} = delta;

            tens = Tensor([pspace' pspace'], [pspace' pspace']);
            tens = fill_tensor(tens, num2cell([1 0 -1 -1 0 1]));
            corr1 = contract(AL1, [1 2 3], AC2, [3 4 8], tens, [2 4 5 7], conj(AL1), [1 5 6], conj(AC2), [6 7 8]);
            corr2 = contract(AL2, [1 2 3], AC1, [3 4 8], tens, [2 4 5 7], conj(AL2), [1 5 6], conj(AC1), [6 7 8]);
            corr = (corr1 + corr2)/2;
            corr_list{count} = corr;

            count = count + 1;
        end
    end
end
