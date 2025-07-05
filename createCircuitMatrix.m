function diagMat = createCircuitMatrix()
% createCircuitMatrixWithConnections - Creează o reprezentare matricială
% simplificată a diagramei din modelul Simulink specificat, incluzând
% atât blocurile R, L, C, cât și legăturile dintre ele (marcate cu '1').
%
% modelName - numele modelului Simulink (modelul trebuie să fie deschis)
%
% diagMat - matrice de caractere în care:
%           'R', 'L', 'C' reprezintă elementele,
%           '1' reprezintă legăturile (conexiuni),
%           '0' reprezintă celule libere.
%
% Exemplu de utilizare:
%   open_system('numele_modelului');
%   M = createCircuitMatrixWithConnections('numele_modelului');
%   disp(M);

% 1. Găsește blocurile R, L și C
modelName=bdroot()
rBlocks = find_system(modelName, 'Regexp', 'on', 'Name', '^Resistor');
lBlocks = find_system(modelName, 'Regexp', 'on', 'Name', '^Inductor');
cBlocks = find_system(modelName, 'Regexp', 'on', 'Name', '^Capacitor');
uyBlocks = find_system(modelName, 'Regexp', 'on', 'Name', '^Voltage');


allBlocks = [rBlocks; lBlocks; cBlocks; uyBlocks];

if isempty(allBlocks)
    error('Nu s-au găsit blocuri R, L sau C în modelul %s.', modelName);
end

% 2. Extrage pozițiile (centrul) și tipurile blocurilor
centers = [];   % va conține [X, Y] pentru fiecare bloc
types   = {};   % va conține litera corespunzătoare (ex: 'R', 'L', 'C')

for i = 1:length(allBlocks)
    try
        pos = get_param(allBlocks{i}, 'Position');  % Format: [left top right bottom]
    catch
        continue;  % dacă blocul nu are 'Position', se sare peste el
    end
    % Calculează centrul
    centerX = (pos(1) + pos(3)) / 2;
    centerY = (pos(2) + pos(4)) / 2;
    centers = [centers; centerX, centerY];
    
    % Presupunem că primul caracter din nume definește tipul
    blkName = get_param(allBlocks{i}, 'Name');
    if ~isempty(blkName)
        if(blkName(end)=='r' || blkName(end) == 'e')
            types{end+1} = [blkName(1) '0'];
        else
             types{end+1} = [blkName(1) blkName(end)];  % ex.: 'R', 'L', 'C'
        end
    else
        types{end+1} = "0";
    end
end

if isempty(centers)
    error('Niciun bloc valid cu parametrul "Position" nu a fost găsit.');
end


% 2.5: Colectarea tuturor coordonatelor (blocuri + linii)
all_line_points = [];
lineHandles = find_system(modelName, 'FindAll', 'on', 'type', 'line');
for i = 1:length(lineHandles)
    pts = get_param(lineHandles(i), 'Points');
    if ~isempty(pts)
        all_line_points = [all_line_points; pts];
    end
end

% Combină coordonatele centrelor blocurilor cu cele ale punctelor de pe linii
if isempty(all_line_points)
    all_coords = centers;
else
    all_coords = [centers; all_line_points];
end


% 3. Definirea grilei pe baza limitelor globale
gridSize = 40;  % dimensiunea celulei în pixeli

% Calculează limitele din TOATE coordonatele (blocuri + linii)
minX = floor(min(all_coords(:,1)) / gridSize);
maxX = ceil(max(all_coords(:,1)) / gridSize);
minY = floor(min(all_coords(:,2)) / gridSize);
maxY = ceil(max(all_coords(:,2)) / gridSize);

numCols = maxX - minX + 1;
numRows = maxY - minY + 1;


% 4. Plasează blocurile în matrice și salvează coordonatele lor
diagMat = repmat("0", numRows, numCols); % Am scos +1 pentru o grilă exactă

% Creează un Map pentru a stoca coordonatele blocurilor
blockGridCoords = containers.Map('KeyType','char','ValueType','any');

for i = 1:size(centers,1)
    colIdx = floor(centers(i,1) / gridSize) - minX + 1;
    rowIdx = floor(centers(i,2) / gridSize) - minY + 1;
    
    diagMat(rowIdx, colIdx) = types{i};
    
    % Salvează coordonatele în hartă, folosind calea completă a blocului ca cheie
    blockGridCoords(allBlocks{i}) = [colIdx, rowIdx];
end


% 5. Desenează conexiunile (corpul liniei + ancorarea la blocuri)
% lineHandles = find_system(modelName, 'FindAll', 'on', 'type', 'line');

for i = 1:length(lineHandles)
    lh = lineHandles(i);
    
    % Obține geometria completă a liniei
    pts = get_param(lh, 'Points');
    if isempty(pts) || size(pts,1) < 2, continue, end
    
    % --- Pasul A: Desenează corpul liniei (logica ta originală) ---
    for k = 1:size(pts,1)-1
        p1 = pts(k,:);
        p2 = pts(k+1,:);
        
        col1 = floor(p1(1) / gridSize) - minX + 1;
        row1 = floor(p1(2) / gridSize) - minY + 1;
        col2 = floor(p2(1) / gridSize) - minX + 1;
        row2 = floor(p2(2) / gridSize) - minY + 1;
        
        coords = bresenham(col1, row1, col2, row2);
        
        for j = 1:size(coords,1)
            rIdx = coords(j,2);
            cIdx = coords(j,1);
            if rIdx >= 1 && rIdx <= numRows && cIdx >= 1 && cIdx <= numCols
                if diagMat(rIdx, cIdx) == "0"
                    diagMat(rIdx, cIdx) = "1";
                end
            end
        end
    end
    
    % --- Pasul B: Ancorează capătul liniei de blocurile destinație ---
    dstBlockHandles = get_param(lh, 'DstBlockHandle');
    if all(dstBlockHandles == -1), continue, end
    
    % Coordonatele de pe grilă ale capătului liniei
    line_end_point = pts(end,:);
    line_end_col = floor(line_end_point(1) / gridSize) - minX + 1;
    line_end_row = floor(line_end_point(2) / gridSize) - minY + 1;
        
    for k = 1:length(dstBlockHandles)
        dstBlkH = dstBlockHandles(k);
        if dstBlkH == -1, continue, end
        
        % Găsește calea completă a blocului destinație
        dstBlkPath = getfullname(dstBlkH);
        
        % Verifică dacă avem coordonatele pentru acest bloc
        if isKey(blockGridCoords, dstBlkPath)
            % Obține coordonatele salvate ale blocului
            block_coords = blockGridCoords(dstBlkPath);
            block_col = block_coords(1);
            block_row = block_coords(2);
            
            % Desenează o linie de legătură de la capătul liniei la bloc
            anchor_coords = bresenham(line_end_col, line_end_row, block_col, block_row);
                  
        end
    end
end

end

