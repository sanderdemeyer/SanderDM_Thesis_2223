files = dir(('Data structures/XXZ_different_deltas_interpretation'));
l = length(files);
count = 1;

for i = 3:l
    disp('ok');
    file = files(i);
    name = file.name;
    disp(name);
    len = length(name);
    disp(name(28:len-40));
    try delta = str2double(name(28:len-40));
        if name(len-8:len-4) == 'final'
            disp(delta);
            deltas{count} = delta;
            %load(name);
            %corr_list = get_spin_correlation_Heisenberg(gs_mps, pspace, 500);
            %corr_list_list{count} = corr_list;
            count = count + 1;
        end
    catch
        fprintf('file was ignored: %s \n', name);
    end
end

%%
scatter(log(i:500), log(abs(corr_list_list{9})))