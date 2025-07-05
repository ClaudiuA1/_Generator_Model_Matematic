function [KCL, KVL] = generateKirchhoffFromMatrix(mat)
% generateKirchhoffFromMatrix - Extrage un set complet, independent și
% robust de ecuații Kirchhoff dintr-o matrice de circuit.
%
% Intrări:
%   mat – M×N string array, reprezentarea circuitului.
%
% Ieșiri:
%   KCL – Cell array cu N-1 ecuații KCL independente.
%   KVL – Cell array cu M-N+1 ecuații KVL fundamentale, sortate
%         pentru a favoriza buclele mici.

    [M, N] = size(mat);

    %% 1. IDENTIFICAREA NODURILOR
    % Se folosește flood-fill pentru a eticheta fiecare zonă de conductor ('1')
    % cu un număr de nod unic. Rezultatul este 'nodeMap'.
    nodeMap = zeros(M, N);
    nodeCount = 0;
    for i = 1:M
        for j = 1:N
            if mat(i,j) == "1" && nodeMap(i,j) == 0
                nodeCount = nodeCount + 1;
                % Folosim o coadă pentru un BFS-fill, mai eficient
                q = [i j];
                nodeMap(i,j) = nodeCount;
                head = 1;
                while head <= size(q,1)
                    [r,c] = deal(q(head,1), q(head,2));
                    head = head+1;
                    for d = [-1 0; 1 0; 0 -1; 0 1]'
                        rr = r + d(1); cc = c + d(2);
                        if rr>=1 && rr<=M && cc>=1 && cc<=N && mat(rr,cc)=="1" && nodeMap(rr,cc)==0
                            nodeMap(rr,cc) = nodeCount;
                            q(end+1,:) = [rr cc]; %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end

    %% 2. EXTRAGEREA RAMURILOR (CU DIRECȚIE CANONICĂ)
    % Se identifică toate componentele și nodurile la care se conectează.
    % Se impune o direcție de referință: de la nodul cu index mai mic la cel mai mare.
    branches = {};
    for i = 1:M
        for j = 1:N
            lbl = char(mat(i,j));
            if ~strcmp(lbl,'1') && ~strcmp(lbl,'0')
                neigh = getAdjacentNodes(nodeMap, i, j, M, N);
                if numel(neigh) == 2
                    branches{end+1} = {lbl, min(neigh), max(neigh)}; %#ok<AGROW>
                end
            end
        end
    end
    % Eliminăm ramurile duplicate (dacă o componentă ocupă mai multe celule)
    [~, unique_indices] = unique(cellfun(@(c) c{1}, branches, 'UniformOutput', false));
    branches = branches(unique_indices);
    numB = numel(branches);
    visualiseGraph(branches);



%% 3) Generare KCL (N ecuații), cu logică de semne și excepții
 KCL = {};
    for n = 1:(nodeCount)
        idx = find(cellfun(@(b) b{2}==n || b{3}==n, branches));
        if isempty(idx), continue; end
        
        labels = cellfun(@(b)b{1}, branches(idx), 'UniformOutput',false);
        
        % Regula 1: Verificăm dacă nodul e conectat la o sursă de INTRARE
        is_input_node = any(ismember(labels, 'V0'));
        if is_input_node
            % Ignorăm complet nodurile de intrare pentru KCL
            continue;
        end
        
        % Regula 2: Filtrăm sursele de IEȘIRE (curent zero)
        is_output_branch = ismember(labels, 'V1');
        idx_for_kcl = idx(~is_output_branch); % Păstrăm tot ce nu e ieșire

        % Regula 4: Generăm ecuația doar dacă avem cel puțin 2 termeni
        if numel(idx_for_kcl) < 2
            continue;
        end
        
        terms = {};
        for k = 1:numel(idx_for_kcl)
            branch_idx = idx_for_kcl(k);
            branch_label = branches{branch_idx}{1};
            branch_n_min = branches{branch_idx}{2};
            
            if n == branch_n_min, terms{end+1}=['+I_' branch_label];
            else, terms{end+1}=['-I_' branch_label]; end
        end
        
        if ~isempty(terms)
            if startsWith(terms{1},'+'), terms{1}(1)=[]; end
            KCL{end+1} = [strjoin(terms, ' ') ' = 0'];
        end
    end
    KCL=KCL(:);

    %% 4. CREAREA GRAFULUI ȘI A ARBORELUI DE ACOPERIRE (SPANNING TREE)
    Gfull = cell(nodeCount,1);
    for b = 1:numB
        n1 = branches{b}{2};  n2 = branches{b}{3};
        Gfull{n1}(end+1) = n2;
        Gfull{n2}(end+1) = n1;
    end
    
    parent = zeros(1,nodeCount);
    visited = false(1,nodeCount);
    q = 1; visited(1)=true;
    head = 1;
    while head <= numel(q)
        u = q(head); head = head+1;
        for v = Gfull{u}
            if ~visited(v), visited(v)=true; parent(v)=u; q(end+1)=v; end
        end
    end

    %% 5. IDENTIFICAREA RAMURILOR DIN ARBORE (inTree) ȘI A COARDELOR (chords)
    inTree = false(1,numB);
    for b = 1:numB
        n1 = branches{b}{2}; n2 = branches{b}{3};
        if parent(n2)==n1 || parent(n1)==n2, inTree(b)=true; end
    end
    chords = find(~inTree);
    
%% 5.5. ETAPA 1 KVL: TRATAREA SPECIALĂ A RAMURILOR PARALELE
parallel_kvl_terms = {};
canonical_keys = cell(numB, 1);
for k = 1:numB
    canonical_keys{k} = sprintf('%d_%d', min(branches{k}{2:3}), max(branches{k}{2:3}));
end

unique_canonical_keys = unique(canonical_keys);
for i = 1:numel(unique_canonical_keys)
    key = unique_canonical_keys{i};
    idxs = find(strcmp(canonical_keys, key));
    
    if numel(idxs) > 1
        % Am găsit un grup de ramuri paralele
        ref_branch_label = branches{idxs(1)}{1};
        for j = 2:numel(idxs)
            current_branch_label = branches{idxs(j)}{1};
            equation = ['U_' ref_branch_label ' - U_' current_branch_label ' = 0'];
            parallel_kvl_terms{end+1} = equation; %#ok<AGROW>
        end
        
        % Logica corectată și îmbunătățită
        % Dintre ramurile paralele, le eliminăm din `chords` doar pe cele
        % pentru care am scris deja o ecuație (adică de la a doua încolo).
        % Păstrăm prima ramură (idxs(1)) în caz că este o coardă necesară
        % pentru o buclă fundamentală mare.
        if numel(idxs) > 1
            chords_to_remove = intersect(idxs(2:end), chords);
            chords = setdiff(chords, chords_to_remove);
        end
    end
end

    %% 6. PREGĂTIRE PENTRU KVL: CALCUL treeAdj ȘI SORTARE CORZI
    % 6.1. Construim lista de adiacență a arborelui O SINGURĂ DATĂ
    treeAdj = cell(nodeCount, 1);
    for b = find(inTree)
        n1 = branches{b}{2}; n2 = branches{b}{3};
        treeAdj{n1}(end+1) = n2; treeAdj{n2}(end+1) = n1;
    end
    
    % 6.2. Sortăm coardele după lungimea căii din arbore pentru a favoriza buclele mici
    chord_distances = zeros(1, numel(chords));
    for k = 1:numel(chords)
        c = chords(k);
        path = findPathInTree(treeAdj, branches{c}{2}, branches{c}{3}, nodeCount);
        chord_distances(k) = length(path); % Lungimea căii în noduri
    end
    [~, sorted_order] = sort(chord_distances);
    chords = chords(sorted_order);

%% 7. GENERARE KVL (CU SEMNE ROBUSTE)
KVL_complex = cell(numel(chords), 1);
kvl_count = 0;
for c = chords
    chord_n_min = branches{c}{2};
    chord_n_max = branches{c}{3};
    
    path = findPathInTree(treeAdj, chord_n_min, chord_n_max, nodeCount);
    
    % VERIFICARE CORECTATĂ:
    if length(path) < 2, continue; end
    
    terms = {['-U_' branches{c}{1}]};
    for k = 1:(numel(path) - 1)
        path_n1 = path(k);
        path_n2 = path(k+1);
        branch_idx = findBranchInTree(branches, inTree, path_n1, path_n2);
        if isempty(branch_idx), continue; end
        
        branch_label = branches{branch_idx}{1};
        branch_n_min = branches{branch_idx}{2};
        
        if path_n1 == branch_n_min
            terms{end+1} = ['+U_' branch_label];
        else
            terms{end+1} = ['-U_' branch_label];
        end
    end
    
    kvl_count = kvl_count + 1;
    KVL_complex{kvl_count} = [strjoin(terms, ' ') ' = 0'];
end
KVL_complex = KVL_complex(1:kvl_count);

% La final, combinăm cele două seturi de ecuații KVL
KVL = [parallel_kvl_terms(:); KVL_complex(:)];

% Afișare rezultate
disp('Ecuații KCL:'); disp(KCL);
disp('Ecuații KVL:'); disp(KVL);

end

% =========================================================================
%                         FUNCȚII AJUTĂTOARE
% =========================================================================

function nodes = getAdjacentNodes(nodeMap, i, j, rows, cols)
% Găsește numerele de nod unice din celulele adiacente (sus, jos, stânga, dreapta)
  nodes = [];
  % Verificăm fiecare direcție și adăugăm nodul dacă există
  if i>1 && nodeMap(i-1,j)>0, nodes(end+1) = nodeMap(i-1,j); end
  if i<rows && nodeMap(i+1,j)>0, nodes(end+1) = nodeMap(i+1,j); end
  if j>1 && nodeMap(i,j-1)>0, nodes(end+1) = nodeMap(i,j-1); end
  if j<cols && nodeMap(i,j+1)>0, nodes(end+1) = nodeMap(i,j+1); end
  nodes = unique(nodes);
end

function path = findPathInTree(treeAdj, startNode, endNode, nodeCount)
% Găsește calea unică între două noduri într-un arbore folosind BFS.
    q = java.util.LinkedList();
    q.add(startNode);
    visited = false(1, nodeCount);
    visited(startNode) = true;
    parent = zeros(1, nodeCount);
    path_found = false;

    while ~q.isEmpty()
        u = q.remove();
        if u == endNode, path_found = true; break; end
        for v_idx = 1:length(treeAdj{u})
            v = treeAdj{u}(v_idx);
            if ~visited(v), visited(v) = true; parent(v) = u; q.add(v); end
        end
    end
    
    if ~path_found, path = []; return; end
    
    path = endNode;
    curr = endNode;
    while curr ~= startNode, curr = parent(curr); path(end+1) = curr; end
    path = fliplr(path);
end

function idx = findBranchInTree(branches, inTree, n1, n2)
% Găsește indexul ramurii din arbore care conectează nodurile n1 și n2.
    idx = [];
    for k=find(inTree) % Caută mult mai eficient, doar printre ramurile din arbore
        b_n1 = branches{k}{2};
        b_n2 = branches{k}{3};
        if (b_n1 == n1 && b_n2 == n2) || (b_n1 == n2 && b_n2 == n1)
            idx = k;
            return;
        end
    end
end