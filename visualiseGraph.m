function visualiseGraph(branches)
    % branches: cell array { {label, n1, n2}, ... }

    % 1) Extragem capetele muchiilor și etichetele
    numB = numel(branches);
    s = zeros(1, numB);
    t = zeros(1, numB);
    edgeLabels = cell(1, numB);
    for k = 1:numB
        edgeLabels{k} = branches{k}{1};
        s(k) = branches{k}{2};
        t(k) = branches{k}{3};
    end

    % 2) Construim graful
    G = graph(s, t);

    % 3) Plotăm graful, folosind edgeLabels și NodeLabel implicit (1:N)
    figure;
    p = plot(G, ...
        'Layout', 'layered', ...
        'EdgeLabel', edgeLabels, ...
        'NodeLabel', arrayfun(@num2str, 1:numnodes(G), 'UniformOutput', false), ...
        'LineWidth', 1.5, ...
        'MarkerSize', 8);
    title('Topologia circuitului');
end
