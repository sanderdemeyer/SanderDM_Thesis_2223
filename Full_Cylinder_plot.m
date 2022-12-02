files = dir(('Data structures/FullCylinder_N_2'));
l = length(files);
disp(l);
[pspace, vspace, trivspace, fusion_trees] = get_spaces('Heisenberg XXZ');

stag_magn = zeros(1, (l-2)/2);
h_stag = zeros(1, (l-2)/2);
j = 1;
for i = 3:l
    disp(i);
    file = files(i);
    name = file.name;
    name_length = length(name);
    disp(name);
    %load(name);
    if strcmp(name(name_length-8:name_length-4), 'final')
        load(strcat('Data structures/FullCylinder_N_2/', name))
        %gs_mps = canonicalize(mps, 'Order', 'rl');
        disp('here');
        disp(name(49:name_length-10));
        h = str2double(name(49:name_length-10));
        disp(h);
        h_stag(j) = h;
        m_stag = get_magnetisation('XXZ', gs_mps, pspace, trivspace, 2.1, false, true);
        disp(m_stag);
        stag_magn(j) = m_stag;
        j = j+1;
%         disp('no, here');
%         disp(name(49:name_length-10));
%         h = str2double(name(49:name_length-10));
%         disp(i);
%         h_stag(i-2) = h;
    end
end

%%

scatter(h_stag, abs(stag_magn));
xlabel('staggered magnetic field $h^{stag}$', 'interpreter', 'latex');
ylabel('staggered magnetisation $m^{stag}$', 'interpreter', 'latex');
title('Staggered magnetisation for cylinder of radius N = 4')


disp('done');