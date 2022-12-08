function c = contraction_general(list)
    disp('ok');
    len = length(list)/2;
    disp(len);
    if len == 1
        c = 0;
    elseif len == 4
        c = contract(list{1}, list{2}, list{3}, list{4}, list{5}, list{6}, list{7}, list{8});
    elseif len == 5
        disp(list(1));
        disp(list(2));
        c = contract(list{1}, list{2}, list{3}, list{4}, list{5}, list{6}, list{7}, list{8}, list{9}, list{10});
    elseif len == 6
        c = contract(list{1}, list{2}, list{3}, list{4}, list{5}, list{6}, list{7}, list{8}, list{9}, list{10}, list{11}, list{12});
    elseif len == 8
        c = contract(list{1}, list{2}, list{3}, list{4}, list{5}, list{6}, list{7}, list{8}, list{9}, list{10}, list{11}, list{12}, list{13}, list{14}, list{15}, list{16});
    elseif len == 10
        c = contract(list{1}, list{2}, list{3}, list{4}, list{5}, list{6}, list{7}, list{8}, list{9}, list{10}, list{11}, list{12}, list{13}, list{14}, list{15}, list{16}, list{17}, list{18}, list{19}, list{20});
    elseif len == 12
        c = contract(list{1}, list{2}, list{3}, list{4}, list{5}, list{6}, list{7}, list{8}, list{9}, list{10}, list{11}, list{12}, list{13}, list{14}, list{15}, list{16}, list{17}, list{18}, list{19}, list{20}, list{21}, list{22}, list{23}, list{24});
    end
end