function truncation_error(mps, d, w)
    error('This makes no sense. Do not use this function please');
    arguments
        mps
        d = 1:depth(mps)
        w = 1:period(mps)
    end
    
    for dd = 1:length(d)
        for ww = 1:length(w)
            [svals, charges] = schmidt_values(mps(dd), w(ww));
            truncation = 1;
            for i = 1:length(svals)
                truncation = truncation - (sum(cell2mat(svals(i))))^2;
            end
            disp(truncation);
        end
    end
end
