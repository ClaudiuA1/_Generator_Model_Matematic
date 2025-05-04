function diagMat = createCircuitMatrixWithConnections()
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
modelName='CalCulator_ModelMatematic_Test';
rBlocks = find_system(modelName, 'Regexp', 'on', 'Name', '^R');
lBlocks = find_system(modelName, 'Regexp', 'on', 'Name', '^I');
cBlocks = find_system(modelName, 'Regexp', 'on', 'Name', '^C');
uBlocks = find_system(modelName, 'Regexp', 'on', 'Name', '^V');
yBlocks = find_system(modelName, 'Regexp', 'on', 'Name', '^E');



allBlocks = [rBlocks; lBlocks; cBlocks; uBlocks;yBlocks];

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
        types{end+1} = [blkName(1) blkName(end)];  % ex.: 'R', 'L', 'C'
    else
        types{end+1} = "0";
    end
end

if isempty(centers)
    error('Niciun bloc valid cu parametrul "Position" nu a fost găsit.');
end

% 3. Definirea grilei
gridSize = 40;  % dimensiunea celulei în pixeli
minX = floor(min(centers(:,1)) / gridSize);
maxX = ceil(max(centers(:,1)) / gridSize);
minY = floor(min(centers(:,2)) / gridSize);
maxY = ceil(max(centers(:,2)) / gridSize);

numCols = maxX - minX + 1;
numRows = maxY - minY + 1;

% Inițial, matricea este completată cu '0' (celule libere)
diagMat = repmat("0", numRows, numCols);

% 4. Plasează blocurile în matrice
% Inițial, creează un cell array completat cu stringuri goale
%diagMat = repmat({'0'}, numRows, numCols);

for i = 1:size(centers,1)
    % Mapare din coordonate în pixeli la indici de grilă:
    colIdx = floor(centers(i,1) / gridSize) - minX + 1;
    rowIdx = floor(centers(i,2) / gridSize) - minY + 1;
    % Inversăm rândurile pentru a avea 0 sus (coincide cu coordonatele vizuale)
    rowIdx = numRows - rowIdx + 1;
    
    diagMat(rowIdx, colIdx) = types{i};
end

% 5. Desenează conexiunile dintre blocuri
% Obține toate liniile din model
lineHandles = find_system(modelName, 'FindAll', 'on', 'type', 'line');

for i = 1:length(lineHandles)
    try
        pts = get_param(lineHandles(i), 'Points');  % Matrice de puncte [x y]
    catch
        continue;
    end
    if isempty(pts) || size(pts,2) < 2
        continue;
    end
    
    % Pentru fiecare segment de linie (între punctele consecutive)
    for k = 1:size(pts,1)-1
        p1 = pts(k,:);
        p2 = pts(k+1,:);
        
        % Mapare la indici de grilă pentru p1
        col1 = floor(p1(1) / gridSize) - minX + 1;
        row1 = floor(p1(2) / gridSize) - minY + 1;
        row1 = numRows - row1 + 1;
        
        % Mapare la indici de grilă pentru p2
        col2 = floor(p2(1) / gridSize) - minX + 1;
        row2 = floor(p2(2) / gridSize) - minY + 1;
        row2 = numRows - row2 + 1;

        % Obține toate celulele de pe segmentul dintre (col1, row1) și (col2, row2)
        coords = bresenham(col1, row1, col2, row2);
        
        % Asigură-te că ultima coordonată (endpoint-ul) este inclusă
        if isempty(coords) || ~isequal(coords(end,:), [col2, row2])
            coords(end+1,:) = [col2, row2];
        end
        
        for j = 1:size(coords,1)
            rIdx = coords(j,2);
            cIdx = coords(j,1);
            % Verifică dacă coordonata este validă în matrice
            if rIdx >= 1 && rIdx <= numRows && cIdx >= 1 && cIdx <= numCols
                % Dacă celula este liberă ('0'), o marchez cu '1'
                if diagMat(rIdx, cIdx) == "0"
                    diagMat(rIdx, cIdx) = "1";
                end
            end
        end
    end
end

end

