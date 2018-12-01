function collection = collect_var_results(type_index, cell_array)
    % Type index
    % 1 for time
    % 2 for p
    % 3 for e
    % 4 for s
    % 5 for es
    % Collects only end results
    [iter, run] = size(cell_array);
    collection = cell(1,run);
    for r = 1:run
        fill = [];
        for i = 1:iter
            fill = [fill;cell_array{i, r}(type_index)];
        end
        collection{1,r} = fill;
    end
end
        