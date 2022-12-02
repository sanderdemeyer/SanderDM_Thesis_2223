charges = U1([1 -1]);
pspace = GradedSpace.new(charges, [1 1], false);
trivspace = GradedSpace.new(U1(0), 1, false);

disp('Interpretation started');

files = dir(('Data structures/XXX_stagh_N_3'));
l = length(files);
count = 0;
h_list = zeros(1, l-3);
stag_magn_list = zeros(1, l-3);
E_list = zeros(1, l-3);
for i = 3:l
    disp('ok');
    file = files(i);
    name = file.name;
    len = length(name);
    if strcmp('.mat', name(len-3:len)) && strcmp('final', name(len-8:len-4))
        disp(name);
        load(name);
        if ~strcmp('final', name(len-8:len-4))
            disp('canonicalize');
            gs_mps = canonicalize(mps, 'Order', 'rl');
        else
            disp('not can');
        end
        h_field = str2double(name(47:len-10));
        disp(h_field);
        stag_magn_new = get_magnetisation('XXZ', gs_mps, pspace, trivspace, 2, false, true);
        stag_magn_list(i-2) = (stag_magn_new);
        h_list(i-2) = h_field;
        %{
        AC1 = gs_mps.AC(1);
        AC2 = gs_mps.AC(2);
        AL1 = gs_mps.AL(1);
        AL2 = gs_mps.AL(2);

        E1 = contract(AL1, [1 2 3], AC2, [3 4 5], conj(AL1), [1 6 7], conj(twist(AC2,3)), [7 8 5], H_A, [2 4 6 8]);
        E2 = contract(AL2, [1 2 3], AC1, [3 4 5], conj(AL2), [1 6 7], conj(twist(AC1,3)), [7 8 5], H_B, [2 4 6 8]);
        E = (E1+E2)/2;
        E_list(i-2) = E;
        %fprintf('New staggered magnetisation is %s \n', stag_magn_new)
        %}
        count = count + 1;
    end
end
%%
scatter(h_list, (stag_magn_list));
xlabel('staggered magnetic field');
ylabel('staggered magnetisation')
xlim([0.8 1]);
title('Staggered magnetisation for N = 2')