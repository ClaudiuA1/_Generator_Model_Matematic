
function ecuatii = extractKirchhoff(mat)
% extractKirchhoff - Extrage ecuațiile Kirchhoff (KCL și KVL) dintr-o matrice de circuit
%   mat      : cell array m×n cu '0','1', și denumiri de componente ('R1','C1','Ve', etc.)
%   ecuatii  : structură cu două câmpuri: ecuatii.KCL și ecuatii.KVL (cell arrays de stringuri)

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

%% 3. KCL: legea curentului pentru fiecare nod cu >=2 ramuri
KCL = {};
for n=1:nodeCount
    inc = cellfun(@(b) b{2}==n || b{3}==n, branches);
    if sum(inc)>1
    % Excludem sursele de tensiune (Ve, E1, etc.) din legea curentului
    valid = branches(inc);
    valid = valid(~cellfun(@(b) startsWith(b{1}, 'Ve') || startsWith(b{1}, 'E'), valid));
    if numel(valid)>1
        terms = cellfun(@(b)['i' b{1}], valid, 'UniformOutput',false);
        KCL{end+1} = strjoin(terms, ' = ');
    end
end
end

%% 4. Pregătiri pentru KVL
numB = numel(branches);
% Graful complet
Gfull = cell(nodeCount,1);
for b=1:numB
    n1=branches{b}{2}; n2=branches{b}{3};
    Gfull{n1}(end+1)=n2; Gfull{n2}(end+1)=n1;
end

% Arbore pe graful complet (DFS)
parentFull = zeros(1,nodeCount);
visited = false(1,nodeCount);
stack = 1; visited(1)=true;
while ~isempty(stack)
    u = stack(end); stack(end) = [];
    for v = Gfull{u}
        if ~visited(v)
            visited(v)=true;
            parentFull(v)=u;
            stack(end+1)=v;
        end
    end
end

% Determinarea chords (ramuri nu în arborele minimal, dar aici doar semnalăm)
inTree = false(1,numB);
for b=1:numB
    n1=branches{b}{2}; n2=branches{b}{3};
    if parentFull(n2)==n1 || parentFull(n1)==n2
        inTree(b)=true;
    end
end
chords = find(~inTree);

%% 5. KVL: bucle fundamentale pe baza chords + arbore complet
KVL = {};
for c = chords
    n1 = branches{c}{2}; n2 = branches{c}{3};
    % Drumul n1->n2 în arbore complet
    path = n1; u = n1;
    while (u~=n2 &&u>0)
        u = parentFull(u);
        path(end+1) = u;
    end
    % Ramurile buclei: chord + segmente din arbore
    loopB = c;
    for k=1:length(path)-1
        a=path(k); bnode=path(k+1);
        for bidx=1:numB
            x=branches{bidx}{2}; y=branches{bidx}{3};
            if (x==a&& y==bnode) || (x==bnode && y==a)
                if inTree(bidx)
                    loopB(end+1)=bidx;
                end
                break;
            end
        end
    end
 terms = {};
for idx = loopB
    lbl = branches{idx}{1};
    if startsWith(lbl,'Ve') || startsWith(lbl,'E')
        terms{end+1} = ['-U' lbl];
    else
        terms{end+1} = ['U' lbl];
    end
end
    KVL{end+1} = [strjoin(terms,' + ') ' = 0'];
end

%% 6. Return și afișare
ecuatii.KCL = KCL;
ecuatii.KVL = KVL;
disp('Ecuații KCL:'); disp(KCL');
disp('Ecuații KVL:'); disp(KVL');
end



% function ecuatii = extractKirchhoff(mat)
% % extractKirchhoff - Extrage ecuațiile Kirchhoff (KCL și KVL) dintr-o matrice de circuit
% %   mat      : cell array m×n cu '0','1', și denumiri de componente ('R1','C1','Ve', etc.)
% %   ecuatii  : structură cu două câmpuri: ecuatii.KCL și ecuatii.KVL (cell arrays de stringuri)
% 
% [rows, cols] = size(mat);
% nodeMap = zeros(rows, cols);
% nodeCount = 0;
% 
% %% 1. Flood-fill pentru identificarea nodurilor (doar fire '1')
% for i=1:rows
%     for j=1:cols
%         if strcmp(mat{i,j}, '1') && nodeMap(i,j)==0
%             nodeCount = nodeCount + 1;
%             cells = floodFill(mat, i, j, rows, cols);
%             for k=1:size(cells,1)
%                 nodeMap(cells(k,1), cells(k,2)) = nodeCount;
%             end
%         end
%     end
% end
% 
% %% 2. Identificarea ramurilor (componente) și nodurile conectate
% branches = {};  % fiecare: {label, node1, node2}
% for i=1:rows
%     for j=1:cols
%         label = mat{i,j};
%         if ~strcmp(label,'0') && ~strcmp(label,'1')
%             neigh = getAdjacentNodes(nodeMap, i, j, rows, cols);
%             if numel(neigh)==2
%                 branches{end+1} = {label, neigh(1), neigh(2)};
%             end
%         end
%     end
% end
% 
% %% 3. KCL: legea curentului pentru fiecare nod cu >=2 ramuri
% KCL = {};
% for n=1:nodeCount
%     inc = cellfun(@(b) b{2}==n || b{3}==n, branches);
%     if sum(inc)>1
%     % Excludem sursele de tensiune (Ve, E1, etc.) din legea curentului
%     valid = branches(inc);
%     valid = valid(~cellfun(@(b) startsWith(b{1}, 'Ve') || startsWith(b{1}, 'E'), valid));
%     if numel(valid)>1
%         terms = cellfun(@(b)['i' b{1}], valid, 'UniformOutput',false);
%         KCL{end+1} = strjoin(terms, ' = ');
%     end
% end
% 
% end
% 
% %% 4. Pregătiri pentru KVL
% numB = numel(branches);
% % Graful complet
% Gfull = cell(nodeCount,1);
% for b=1:numB
%     n1=branches{b}{2}; n2=branches{b}{3};
%     Gfull{n1}(end+1)=n2; Gfull{n2}(end+1)=n1;
% end
% 
% % Arbore pe graful complet (DFS)
% parentFull = zeros(1,nodeCount);
% visited = false(1,nodeCount);
% stack = 1; visited(1)=true;
% while ~isempty(stack)
%     u = stack(end); stack(end) = [];
%     for v = Gfull{u}
%         if ~visited(v)
%             visited(v)=true;
%             parentFull(v)=u;
%             stack(end+1)=v;
%         end
%     end
% end
% 
% % Determinarea chords (ramuri nu în arborele minimal, dar aici doar semnalăm)
% inTree = false(1,numB);
% for b=1:numB
%     n1=branches{b}{2}; n2=branches{b}{3};
%     if parentFull(n2)==n1 || parentFull(n1)==n2
%         inTree(b)=true;
%     end
% end
% chords = find(~inTree);
% 
% %% 5. KVL: bucle fundamentale pe baza chords + arbore complet
% KVL = {};
% for c = chords
%     % 1) găseşti path-ul în arbore complet
% n1 = branches{c}{2};
% n2 = branches{c}{3};
%     path = findPath(Gfull, n1, n2);
%     if isempty(path)
%         error('Nu există drum complet între nodurile %d și %d.', n1, n2);
%     end
% 
% 
%     % 2) traversezi path-ul și semnezi
%     terms = {};
% for k = 1:length(path)-1
%     a = path(k); bnode = path(k+1);
%     for idx = 1:numB
%         br = branches{idx};
%         lbl = br{1};
%         n1  = br{2};
%         n2  = br{3};
%         if n1 == a && n2 == bnode
%             terms{end+1} = ['+U' lbl];
%             break;
%         elseif n1 == bnode && n2 == a
%             terms{end+1} = ['-U' lbl];
%             break;
%         end
%     end
% end
% 
%     eq = strjoin(terms, ' ');
%     if startsWith(eq, '+'), eq = extractAfter(eq,1); end
%     KVL{end+1} = [eq ' = 0'];
% end
% 
% 
% %% 6. Return și afișare
% ecuatii.KCL = KCL;
% ecuatii.KVL = KVL;
% disp('Ecuații KCL:'); disp(KCL');
% disp('Ecuații KVL:'); disp(KVL');
% end
% 
