function ecuatii = extractKirchhoff(mat)
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

%% 2. Identificarea ramurilor (componente) și nodurile conectate
branches = {};  % fiecare: {label, node1, node2}
for i=1:rows
    for j=1:cols
        label = mat{i,j};
        if ~strcmp(label,'0') && ~strcmp(label,'1')
            neigh = getAdjacentNodes(nodeMap, i, j, rows, cols);
            if numel(neigh)==2
                branches{end+1} = {label, neigh(1), neigh(2)};
            end
        end
    end
end
 visualiseGraph(branches)
%% 3. KCL: legea curentului pentru fiecare nod
KCL = {};
for n = 1:nodeCount
    % 3.a) Găsim toate ramurile care ating nodul n (inclusiv surse)
    ramIndices = find(cellfun(@(b) (b{2} == n) || (b{3} == n), branches));
    if isempty(ramIndices)
        continue;  % nod fără ramuri
    end
    
    % 3.b) Excludem complet sursele de tensiune (prefix V sau E)
    isPassive = ~cellfun(@(b) startsWith(b{1}, 'V') || startsWith(b{1}, 'E'), branches(ramIndices));
    nonSrcIdx = ramIndices(isPassive);
    if isempty(nonSrcIdx)
        continue;  % după excluderea surselor nu mai rămâne nicio ramură pasivă
    end
    
    % 3.c) Împărțim ramurile pasive în două categorii:
    %      incomingIdx = acele ramuri care ajung în nod (b{3} == n)
    %      outgoingIdx = acele ramuri care pleacă din nod (b{2} == n)
    incomingMask = cellfun(@(b) (b{3} == n), branches(nonSrcIdx));
    incomingIdx = nonSrcIdx(incomingMask);
    outgoingIdx = nonSrcIdx(~incomingMask);
    
    % 3.d) Scriem ecuația în funcție de numărul de ramuri pasive
    switch numel(nonSrcIdx)
        case 1
            % --- EXACT O SINGURĂ ramură pasivă la acest nod ---
            k = nonSrcIdx;       % indexul acelei singure ramuri
            lbl = branches{k}{1};% eticheta, ex. 'R1'
            
            if ~isempty(incomingIdx) && isempty(outgoingIdx)
                % Situație: exact o singură ramură și ea e incoming
                % Ilbl = suma tuturor curenților outgoing (dar nu există outgoing)
                % => pur și simplu Ilbl = 0 nu are sens fizic în KCL
                % scriem: Ilbl = (nimic)  → implicit Ilbl = 0, dar mai bine ignorăm
                % DECIDEM: nu generăm ecuație în acest caz (de obicei nu apare
                % “doar o singură ramură pasivă” într-un circuit normal).
                continue;
                
            elseif isempty(incomingIdx) && ~isempty(outgoingIdx)
                % Situație: exact o singură ramură și ea e outgoing
                % => (suma incoming) = Ilbl, dar incoming e vidă
                % => 0 = Ilbl, iar Ilbl nu poate fi 0 fizic
                % DECIDEM: nu generăm ecuație
                continue;
                
            elseif ~isempty(incomingIdx) && ~isempty(outgoingIdx)
                % Practic cam imposibil: dacă există o singură ramură pasivă,
                % nu poate fi în același timp și incoming, și outgoing.
                continue;
            end
            
        otherwise
            % --- CEL PUȚIN DOUĂ ramuri pasive la acest nod ---
            % Ecuație generală: suma(incoming) = suma(outgoing)
            
            % Construim lista curenților incoming
            if isempty(incomingIdx)
                lhs = {'0'};
            else
                lhs = cellfun(@(b) ['I' b{1}], branches(incomingIdx), 'UniformOutput', false);
            end
            
            % Construim lista curenților outgoing
            if isempty(outgoingIdx)
                rhs = {'0'};
            else
                rhs = cellfun(@(b) ['I' b{1}], branches(outgoingIdx), 'UniformOutput', false);
            end
            
            % Concatenează cu „ + ” și adaugă la KCL
            KCL{end+1} = sprintf('%s = %s', ...
                strjoin(lhs, ' + '), ...
                strjoin(rhs, ' + '));
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

%% 7) Construim acum buclele KVL – exact câte o ecuație pentru fiecare cordă rămasă
KVL = {};
for ci = 1:numel(chords)
    c = chords(ci);
    n1 = branches{c}{2};
    n2 = branches{c}{3};

    % 7.a) Reconstruim drumul unic n1→n2 în arbore
    %     Mai întâi extragem subgraful arborelui (numai ramurile cu inTree = true)
    treeAdj = cell(nodeCount,1);
    for bidx = 1:numB
        if inTree(bidx)
            a = branches{bidx}{2};
            b_ = branches{bidx}{3};
            treeAdj{a}(end+1) = b_;
            treeAdj{b_}(end+1) = a;
        end
    end

    %     Apoi facem un BFS în subgraful arborelui pentru a găsi părinții temporari
    parTemp = zeros(1, nodeCount);
    vis2    = false(1, nodeCount);
    q2 = n1;
    vis2(n1) = true;
    head2 = 1;
    while head2 <= numel(q2)
        u2 = q2(head2);
        head2 = head2 + 1;
        if u2 == n2
            break;
        end
        for nb = treeAdj{u2}
            if ~vis2(nb)
                vis2(nb) = true;
                parTemp(nb) = u2;
                q2(end+1) = nb;
            end
        end
    end

    %     Reconstruim calea din arbore de la n2 înapoi la n1
    path = n2;
    cur  = n2;
    while cur ~= n1
        cur = parTemp(cur);
        path(end+1) = cur;
    end
    %     Inversăm ca să fie [n1, …, n2]
    path = fliplr(path);

    % 7.b) Colectăm index‐ii ramurilor din buclă: corda c + fiecare segment (a→b_) din arbore
    loopBranches = c;  
    for k = 1:(numel(path)-1)
        a  = path(k);
        b_ = path(k+1);
        % găsim ramura (bidx) care leagă exact a ↔ b_ în subgraful arbore
        for bidx = 1:numB
            x = branches{bidx}{2};
            y = branches{bidx}{3};
            if inTree(bidx) && ((x == a && y == b_) || (x == b_ && y == a))
                loopBranches(end+1) = bidx;
                break;
            end
        end
    end

    % 7.c) Construim ecuația KVL cu semnele corespunzătoare
    terms = {};

    %  7.c.1) Începem cu corda C: bucla este închisă de la n2 → n1 (invers faţă de n1→n2)
    lbl_c = branches{c}{1};
    % Direcția naturală a corzii e n1→n2, dar bucla o parcurge n2→n1 → semn „−”
    terms{end+1} = ['-U' lbl_c];

    %  7.c.2) Apoi parcurgem path de la n1→n2, adăugând fiecare ramură cu semnul
    for k = 1:(numel(path)-1)
        a  = path(k);
        b_ = path(k+1);
        for bidx = 1:numB
            x = branches{bidx}{2};
            y = branches{bidx}{3};
            if (x == a && y == b_)
                % ramura merge a→b_, iar bucla merge tot a→b_ → semn „+”
                lbl = branches{bidx}{1};
                terms{end+1} = ['+U' lbl];
                break;
            elseif (x == b_ && y == a)
                % ramura are sens natural b_→a, dar bucla merge a→b_ → semn „−”
                lbl = branches{bidx}{1};
                terms{end+1} = ['-U' lbl];
                break;
            end
        end
    end

    % 7.d) Încheiem ecuația de buclă
    KVL{end+1} = [ strjoin(terms) ' = 0' ];
end

    %% 7. Returnăm și afișăm

    ecuatii.KCL      = KCL(:);
    ecuatii.KVL      = KVL(:);

    disp('Ecuații KCL:');        disp(ecuatii.KCL);
    disp('Ecuații KVL:');        disp(ecuatii.KVL);

    % extractKVL_general(mat);
end
