function clean_data = rm_extras(memory,start,stop)
    [iter, run] = size(memory);
    clean_data = cell(iter,run);
    for i = 1:iter
        for r = 1:run
            extract = memory{i,r};
            extract(:,start:stop)=[];
            clean_data{i,r} = extract;
        end
    end
    return
end