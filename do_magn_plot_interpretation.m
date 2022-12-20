charges = U1([1 -1]);
pspace = GradedSpace.new(charges, [1 1], false);
trivspace = GradedSpace.new(U1(0), 1, false);

disp('Interpretation started');

files = dir(('Data structures/FullCylinder_N_4_postchange'));
l = length(files);
count = 1;
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
        h_field = str2double(name(len-14:len-10));
        disp(name(len-15:len-10));
        disp(h_field);
        stag_magn = get_magnetisation('XXZ', gs_mps, pspace, trivspace, 1.5, false, true);

        stag_magn_list{count} = stag_magn;
        h_list{count} = h_field;
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
%scatter(h_list, (stag_magn_list));
h_list = [0 cell2mat(h_list)];
stag_magn_list = [0 cell2mat(stag_magn_list)];
scatter((h_list), abs((stag_magn_list)));
xlabel('staggered magnetic field');
ylabel('staggered magnetisation')
title('Staggered magnetisation for N = 4')