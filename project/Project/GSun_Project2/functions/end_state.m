function end_results = end_state(cleaned_data)
    [iter, run] = size(cleaned_data);
    end_results = cell(iter,run);
    for i = 1:iter
        for r = 1:run
            extract = cleaned_data{i,r}(end,:);
            end_results{i,r} = extract;
        end
    end
    return
end
    
