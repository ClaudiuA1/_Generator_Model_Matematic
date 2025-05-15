function path = findPath(G, startNode, endNode)
    % G: cell array de adiacență
    % startNode, endNode: noduri între care cauți drumul
    N = numel(G);
    visited = false(1, N);
    parent = zeros(1, N);
    queue = startNode;
    visited(startNode) = true;
    found = false;
    while ~isempty(queue)
        u = queue(1); queue(1) = [];
        if u == endNode
            found = true;
            break;
        end
        for v = G{u}
            if ~visited(v)
                visited(v) = true;
                parent(v) = u;
                queue(end+1) = v;
            end
        end
    end
    if ~found
        path = [];
        return;
    end
    % Reconstruiește path-ul
    path = endNode;
    while path(1) ~= startNode
        path = [parent(path(1)), path];
    end
end
