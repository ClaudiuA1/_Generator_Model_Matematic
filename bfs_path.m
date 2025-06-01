function path = bfs_path(tree, start, target, num_nodes)
    visited = false(1, num_nodes);
    parent = zeros(1, num_nodes);
    queue = start;
    visited(start) = true;
    found = false;
    
    while ~isempty(queue) && ~found
        current = queue(1);
        queue(1) = [];
        neighbors = tree{current};
        for n = neighbors
            if ~visited(n)
                visited(n) = true;
                parent(n) = current;
                queue(end+1) = n;
                if n == target
                    found = true;
                    break;
                end
            end
        end
    end
    
    % Reconstruct path
    path = [];
    if found
        current = target;
        while current ~= start
            path = [current, path];
            current = parent(current);
        end
        path = [start, path];
    end
end