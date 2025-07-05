function [KCL,KVL] = extractKirchhoff(mat)
% extractKirchhoff - Extrage ecuațiile Kirchhoff (paralel, KCL, KVL)
% dintr-o matrice de circuit “mat” (cell array m×n cu '0','1','R1','C1','Ve', ...).
%
%   mat      : cell array m×n în care:
%                - '0' = celulă liberă
%                - '1' = fir de conexiune
%                - 'R1','C1','L1','Ve','E1',… = componente (ramuri)
%   ecuatii  : struct cu câmpurile:
%                .parallel : cell array de ecuații de felul 'Ua - Ub = 0'
%                .KCL      : cell array de ecuații KCL (curenți)
%                .KVL      : cell array de ecuații KVL (tensiuni)

[rows, cols] = size(mat);
nodeMap = zeros(rows, cols);
nodeCount = 0;
%% 1. Flood-fill pentru identificarea nodurilor (doar fire '1')
for i=1:rows
    for j=1:cols
        if strcmp(mat{i,j}, '1') && nodeMap(i,j)==0 
            nodeCount = nodeCount + 1;
            cells = floodFill(mat, i, j, rows, cols);
            for k=1:size(cells,1)
                nodeMap(cells(k,1), cells(k,2)) = nodeCount;
            end
        end
        
    end
end

%% 2. Identificarea ramurilor (cu direcție canonică)
branches = {};  % fiecare: {label, node_min, node_max}
for i=1:rows
    for j=1:cols
        label = mat{i,j};
        if ~ischar(label), label = char(label); end % Asigură că e char
        if ~strcmp(label,'0') && ~strcmp(label,'1')
            neigh = getAdjacentNodes(nodeMap, i, j, rows, cols);
            if numel(neigh)==2
                % REGULA CANONICĂ: Stocăm mereu nodurile în ordine crescătoare
                node_min = min(neigh);
                node_max = max(neigh);
                branches{end+1} = {label, node_min, node_max};
            end
        end
    end
end
% Eliminăm ramurile duplicate care pot apărea din cauza componentelor întinse
[~, unique_indices] = unique(cellfun(@(c) c{1}, branches, 'UniformOutput', false));
branches = branches(unique_indices);

visualiseGraph(branches);
%% 3. KCL: Legea curentului pentru fiecare nod (metoda robustă)
KCL = {};
for n = 1:nodeCount
    % Găsim toate ramurile conectate la nodul 'n'
    connected_branches = {};
    for k = 1:numel(branches)
        if branches{k}{2} == n || branches{k}{3} == n
            % Excludem sursele de tensiune pure (V, E) din KCL
            label = branches{k}{1};
            if ~startsWith(label, 'V') && ~startsWith(label, 'E')
                connected_branches{end+1} = ['I_' label];
            end
        end
    end
    
    % Dacă nodul are cel puțin o conexiune validă, scriem ecuația
    if numel(connected_branches) > 1
        % Ecuatia este: I_comp1 + I_comp2 + ... = 0
        equation_str = [strjoin(connected_branches, ' + ') ' = 0'];
        KCL{end+1} = equation_str;
    end
end


%% 1) Identificăm ramurile paralele și forțăm păstrarea unei singure ’coarde’ pentru paralel
numB=numel(branches);
% a) Construim o cheie „nesortată” a perechii (n1,n2)
keyUnordered = cell(numB,1);
for b = 1:numB
    n1 = branches{b}{2};
    n2 = branches{b}{3};
    if n1 < n2
        keyUnordered{b} = sprintf('%d_%d', n1, n2);
    else
        keyUnordered{b} = sprintf('%d_%d', n2, n1);
    end
end

% b) Grupăm ramurile care apar sub aceeași cheie
parallelGroups = cell(0);
uniqueKeys = unique(keyUnordered);
for iKey = 1:numel(uniqueKeys)
    idxs = find(strcmp(keyUnordered, uniqueKeys{iKey}));
    if numel(idxs) > 1
        % Avem ≥2 ramuri între aceleași două noduri → sunt paralele
        parallelGroups{end+1} = idxs;  %#ok<AGROW>
    end
end

%% 2) Construim graful neorientat complet (lista de adiacență) 
Gfull = cell(nodeCount, 1);
for b = 1:numB
    n1 = branches{b}{2};
    n2 = branches{b}{3};
    Gfull{n1}(end+1) = n2;
    Gfull{n2}(end+1) = n1;
end

%% 3) BFS pentru a găsi arborele spanning (parentFull)
parentFull = zeros(1, nodeCount);
visited    = false(1, nodeCount);
queue = 1;
visited(1) = true;
headq = 1;
while headq <= numel(queue)
    u = queue(headq);
    headq = headq + 1;
    for v = Gfull{u}
        if ~visited(v)
            visited(v) = true;
            parentFull(v) = u;
            queue(end+1) = v;
        end
    end
end

%% 4) Identificăm toate ramurile care sunt în arbore (inTree = true/false)
inTree = false(1, numB);
for b = 1:numB
    n1 = branches{b}{2};
    n2 = branches{b}{3};
    if parentFull(n2) == n1 || parentFull(n1) == n2
        inTree(b) = true;
    end
end
% Construim lista de muchii a arborelui (numai ramurile cu inTree = true)
sTree = [];
tTree = [];
for b = 1:numB
    if inTree(b)
        sTree(end+1) = branches{b}{2};   % nodul de plecare
        tTree(end+1) = branches{b}{3};   % nodul de sosire
    end
end

% Creăm un obiect graph și îl afișăm
Gtree = graph(sTree, tTree);
figure;
plot(Gtree, ...
     'Layout', 'layered', ...
     'NodeLabel', arrayfun(@num2str, 1:nodeCount, 'UniformOutput', false), ...
     'LineWidth', 1.5, ...
     'MarkerSize', 8);
title('Arborele spanning extras din Gfull');

%% 5) Din fiecare grup de paralele, forțăm ca exact o ramură să rămână „în arbore”,
%    iar pe restul le marcăm drept „chords” (coarde). În acest fel nu mai generăm
%    bucle redundante între C1 și C2, de exemplu.
for g = 1:numel(parallelGroups)
    idxs = parallelGroups{g};
    % De la aceste idxs, fie ramura care este deja inTree = true o păstrăm,
    % fie dacă toate idxs sunt false, îi setăm pe idxs(1) în arbore și restul false.
    % Apoi pe idxs(2:end) îi forțăm să fie coarde (inTree = false).
    
    % (1) Căutăm dacă există deja una „în arbore”
    foundInTree = false;
    for k = idxs
        if inTree(k)
            foundInTree = true;
            break;
        end
    end
    if ~foundInTree
        % Niciuna nu era în arbore → punem pe prima ca în arbore
        inTree(idxs(1)) = true;
    end
    % (2) pe restul idxs(2:end) le forțăm coarde
    for k = idxs
        if k ~= idxs(1)
            inTree(k) = false;
        end
    end
end

% În acest moment, inTree(b) = true dacă ramura b face parte din arbore,
% și false dacă ramura b e coardă (inclusiv paralel “suplimentar”).

%% 6) Determinăm lista de coarde (chords)
chords = find(~inTree);

% Creăm lista de adiacență DOAR pentru ramurile din arborele de acoperire
treeAdj = cell(nodeCount, 1);
for bidx = 1:numel(branches)
    if inTree(bidx)
        n1 = branches{bidx}{2};
        n2 = branches{bidx}{3};
        treeAdj{n1}(end+1) = n2;
        treeAdj{n2}(end+1) = n1;
    end
end

% OPTIMIZARE: Sortăm coardele după lungimea geometrică pentru a favoriza buclele mici

chord_lengths = [];
for k = 1:numel(chords)
    chord_idx = chords(k);
    
    % Găsim nodurile corzii
    n1 = branches{chord_idx}{2};
    n2 = branches{chord_idx}{3};
    
    % Găsim coordonatele de pe grilă ale blocurilor de la capetele corzii
    % (Aici presupunem că fiecare nod are un bloc asociat - ar trebui să fie
    % adevărat pentru capetele unei coarde, care sunt componente)
    
    % Cautam blocurile care definesc nodurile n1 si n2
    % Aceasta este o căutare inversă; poate dura puțin.
    % coord1 = []; coord2 = [];
    % keys = blockGridCoords.keys;
    % for i_key = 1:numel(keys)
    %     coords = blockGridCoords(keys{i_key});
    %     % Trebuie să găsim ce bloc este conectat la nodul n1 și n2
    %     % Acest pas este complicat. O alternativă mai simplă...
    % end
    % 
    % O METODĂ MULT MAI SIMPLĂ:
    % Găsim direct coordonatele nodurilor din matricea 'nodeMap'
    % (prima apariție a numărului nodului)
    [r1, c1] = find(nodeMap == n1, 1, 'first');
    [r2, c2] = find(nodeMap == n2, 1, 'first');
    
    if ~isempty(r1) && ~isempty(r2)
        % Calculăm distanța euclidiană la pătrat (e suficientă pentru sortare)
        dist_sq = (r1-r2)^2 + (c1-c2)^2;
        chord_lengths(end+1) = dist_sq;
    else
        chord_lengths(end+1) = inf; % Coardă fără coordonate, o punem la final
    end
end

% Sortăm indicii corzilor pe baza lungimilor calculate
[~, sorted_indices] = sort(chord_lengths);
chords = chords(sorted_indices);


%% 7) Construim buclele KVL – câte o ecuație pentru fiecare coardă
KVL = {};
for ci = 1:numel(chords)
    chord_idx = chords(ci);
    
    % Nodurile corzii (direcția de referință este n_min -> n_max)
    chord_n_min = branches{chord_idx}{2};
    chord_n_max = branches{chord_idx}{3};
    
    % Găsim drumul unic în arbore între nodurile corzii
    path = findPathInTree(treeAdj, chord_n_min, chord_n_max, nodeCount);
    if isempty(path), continue, end % Nu s-a găsit drum, ar fi o eroare în arbore
    
    % Bucla este formată din drumul din arbore (n_min -> n_max) și coarda (n_max -> n_min)
    
    % Termenul pentru coardă: drumul prin arbore e n_min->n_max, deci bucla
    % se închide prin coardă de la n_max înapoi la n_min. Acesta este sensul
    % INVERS direcției de referință a corzii. Deci, semnul este '-'.
    terms = {['-U_' branches{chord_idx}{1}]};

    % Parcurgem drumul din arbore de la n_min la n_max
    for k = 1:(numel(path)-1)
        path_node1 = path(k);
        path_node2 = path(k+1);
        
        % Găsim ramura din arbore care leagă aceste două noduri
        branch_in_path_idx = findBranch(branches, inTree, path_node1, path_node2);
        if isempty(branch_in_path_idx), continue, end
        
        branch_label = branches{branch_in_path_idx}{1};
        branch_n_min = branches{branch_in_path_idx}{2};
        
        % Stabilim semnul
        % Dacă drumul (path_node1 -> path_node2) este în aceeași direcție
        % cu direcția de referință a ramurii (branch_n_min -> branch_n_max)...
        if path_node1 == branch_n_min
            % ...atunci căderea de tensiune se adună (semn +)
            terms{end+1} = ['+U_' branch_label];
        else
            % ...altfel, mergem împotriva direcției de referință, deci se scade (semn -)
            terms{end+1} = ['-U_' branch_label];
        end
    end
    
    % Curățăm primul termen dacă are un '+' la început
    if startsWith(terms{1}, '+')
        terms{1} = terms{1}(2:end);
    end
    
    KVL{end+1} = [strjoin(terms, ' ') ' = 0'];
end

% Funcții ajutătoare de adăugat la sfârșitul fișierului sau ca funcții private




    %% 8. Returnăm și afișăm

    KCL = KCL(:);
    KVL = KVL(:);

    disp('Ecuații KCL:');        disp(KCL);
    disp('Ecuații KVL:');        disp(KVL);

    % extractKVL_general(mat);
end
