function visualiseDirectedCircuit(branches)
    % branches: cell array { {label, n1, n2}, ... }
    %
    % Construim un graf orientat în care fiecare muchie
    % merge din branches{k}{2} către branches{k}{3}, iar
    % eticheta aferentă e compusă, în caz că există mai
    % multe componente cu aceeași direcție între aceleași noduri.

    % 1) Extragem n1, n2 și etichetele brute din branches
    numB = numel(branches);
    s = zeros(numB,1);      % nodul „de la”
    t = zeros(numB,1);      % nodul „către”
    rawLabels = cell(numB,1);
    for k = 1:numB
        rawLabels{k} = branches{k}{1};  % eticheta componentei
        s(k) = branches{k}{2};          % branches{k}{2} = n1 (origine)
        t(k) = branches{k}{3};          % branches{k}{3} = n2 (destinație)
    end

    % 2) Creăm chei unice pentru perechile orientate n1→n2
    keyPairs = cell(numB,1);
    for k = 1:numB
        % Păstrăm sensul: 'n1_n2'
        keyPairs{k} = sprintf('%d_%d', s(k), t(k));
    end

    % 3) Grupăm indexurile ramurilor după keyPairs
    uniqueKeys = unique(keyPairs);
    numE = numel(uniqueKeys);
    s2 = zeros(numE,1);
    t2 = zeros(numE,1);
    edgeLabels2 = cell(numE,1);

    for i = 1:numE
        key = uniqueKeys{i};
        idxs = find(strcmp(keyPairs, key));  % ramurile care au aceeași cheie exact
        
        % Parsează din „key” cei doi noduri (fără a reordona)
        parts = sscanf(key, '%d_%d');
        u = parts(1);  % nodul sursă
        v = parts(2);  % nodul destinație

        % 4) Construim eticheta compusă din toate ramurile idxs
        labelsGroup = rawLabels(idxs);
        combinedLabel = strjoin(labelsGroup, ' + ');

        % 5) Salvăm muchia orientată u→v și eticheta rezultată
        s2(i) = u;
        t2(i) = v;
        edgeLabels2{i} = combinedLabel;
    end

    % 6) Construim graful orientat (digraph)
    Gd = digraph(s2, t2);

    % 7) Atribuim etichetele drept proprietate EdgeLabel
    % (le vom folosi la plot)
    Gd.Edges.Label = edgeLabels2(:);

    % 8) Plotăm graful orientat
    figure;
    p = plot(Gd, ...
        'Layout', 'layered', ...
        'EdgeLabel', Gd.Edges.Label, ...
        'NodeLabel', arrayfun(@num2str, 1:numnodes(Gd), 'UniformOutput', false), ...
        'LineWidth', 1.5, ...
        'MarkerSize', 8, ...
        'ArrowSize', 10);  % să se vadă bine săgețile
    title('Graf orientat al circuitului (sens branches{2}→branches{3})');
end
