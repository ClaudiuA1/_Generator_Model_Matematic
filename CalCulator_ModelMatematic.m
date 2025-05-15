%CalCulator_ModelMatematic
clc 
clear
close all

objects = createCircuitMatrix()

khirchoffs=extractKirchhoff(objects)


    



%%
% function [KVL_eq, KCL_eq] = extractKirchhoffEquations(mat)
%     % mat este un cell array în care fiecare celulă conține:
%     %   - "0" pentru celule libere
%     %   - "1" pentru linii de conexiune
%     %   - Un string precum "R1", "Ve", "Cr", "I1", "E2" etc. pentru componente
% 
%     % Pasul 1: Identificarea nodurilor folosind o funcție flood fill
%     [rows, cols] = size(mat);
%     nodeMap = zeros(rows, cols);  % fiecare celulă va primi un index de nod
%     nodeIndex = 0;
%     for i = 1:rows
%         for j = 1:cols
%             if ~strcmp(mat{i,j}, "0")  % celulă ocupată (componentă sau linie)
%                 if nodeMap(i,j) == 0
%                     nodeIndex = nodeIndex + 1;
%                     % floodFill: găsește toate celulele conectate de la (i,j)
%                     nodeCells = floodFill(mat, i, j);
%                     % Marchez aceste celule în nodeMap cu indexul nodeIndex
%                     for k = 1:size(nodeCells,1)
%                         nodeMap(nodeCells(k,1), nodeCells(k,2)) = nodeIndex;
%                     end
%                     % Opțional: stochează coordonatele nodului
%                     nodes{nodeIndex} = nodeCells;  %#ok<AGROW>
%                 end
%             end
%         end
%     end
% 
%     % Pasul 2: Identificarea ramurilor (componente) și a nodurilor la care sunt conectate.
%     % Pentru fiecare celulă ce conține o componentă (nu "0" și nu "1")
%     branches = {};  % fiecare ramură: {label, node1, node2}
%     for i = 1:rows
%         for j = 1:cols
%             if ~strcmp(mat{i,j}, "0") && ~strcmp(mat{i,j}, "1")
%                 % Această celulă este o componentă; aflăm la ce noduri e conectată
%                 adjacentNodes = getAdjacentNodes(nodeMap, i, j);
%                 % Presupunem că o componentă conectează două noduri
%                 if numel(adjacentNodes) == 2
%                     branches{end+1} = {mat{i,j}, adjacentNodes(1), adjacentNodes(2)};  %#ok<AGROW>
%                 end
%             end
%         end
%     end
% 
%     % Pasul 3: Construim o matrice de incidență pentru graf
%     % Numărăm nodurile și ramurile
%     numNodes = nodeIndex;
%     numBranches = length(branches);
%     Inc = zeros(numNodes, numBranches);
%     branchLabels = cell(numBranches, 1);
%     for b = 1:numBranches
%         branch = branches{b};
%         branchLabels{b} = branch{1};
%         n1 = branch{2};
%         n2 = branch{3};
%         % Conform convenției: curentul intră în n1 și iese din n2
%         Inc(n1, b) = 1;
%         Inc(n2, b) = -1;
%     end
% 
%     % Pasul 4: Determinăm buclele fundamentale
%     % Aceasta se poate face, de exemplu, prin metoda lui Kirchhoff cu matricea incidență.
%     % Pentru o rețea cu numNodes noduri, numărul de bucle independente este:
%     % numLoops = numBranches - (numNodes - 1);
%     % Apoi, folosind algoritmi de găsire a bazelor ciclurilor (ex. algoritmul lui Horton),
%     % putem obține o bază a buclelor.
%     % Pentru fiecare buclă, scriem ecuația KVL: suma tensiunilor pe ramuri = 0.
%     % Notă: tensiunile pe ramuri vor fi notate ca U<componenta>.
% 
%     % [Acest pas necesită implementare suplimentară – aici doar exemplificăm ideea]
%     KVL_eq = {};  % cell array de stringuri pentru ecuațiile tensiunii
%     % Exemplu: pentru fiecare buclă identificată, compunem ecuația:
%     % "U(Ve) - U(R1) - U(Cr) - U(I1) = 0"
%     % De aici se poate extrage relația: U(Ve) = U(R1) + U(Cr) + U(I1)
% 
%     % Pasul 5: Aplicăm KCL la fiecare nod (nodurile interne, nu cele de referință)
%     KCL_eq = {};  % cell array de stringuri pentru ecuațiile curentului
%     % Pentru fiecare nod, suma curenților care intră = suma curenților care ies.
%     % Dacă pentru o ramură etichetată "R1", curentul este notat i(R1),
%     % atunci pentru nodurile la care se întâlnește se va scrie:
%     % "i(R1) = i(Cr) = i(I1)" sau similar, după caz.
% 
%     % În funcție de rețeaua determinată, aceste ecuații se pot genera automat.
%     % Pentru exemplificare, presupunem că avem:
%     % KVL_eq{1} = "UVe = URr + UCr + UIr";
%     % KCL_eq{1} = "iRr = iCr = iIr";
%     KVL_eq{end+1} = "UVe = URr + UCr + UIr";
%     KCL_eq{end+1} = "iRr = iCr = iIr";
% 
%     % Combinăm rezultatele într-o structură
%     ecuatii.KVL = KVL_eq;
%     ecuatii.KCL = KCL_eq;
% 
%     % Afișăm ecuațiile (exemplu)
%     disp('Ecuații KVL:');
%     disp(KVL_eq');
%     disp('Ecuații KCL:');
%     disp(KCL_eq');
% end
% 
% % --- Funcții Helper ---
% 
% function cellsFound = floodFill(mat, i, j)
%     % Funcție recursivă simplificată pentru a găsi toate celulele conectate de la (i,j)
%     [rows, cols] = size(mat);
%     persistent visited
%     if isempty(visited)
%         visited = zeros(rows, cols);
%     end
%     cellsFound = [];
%     if i < 1 || i > rows || j < 1 || j > cols
%         return;
%     end
%     if visited(i, j) == 1 || strcmp(mat{i,j}, "0")
%         return;
%     end
%     visited(i, j) = 1;
%     cellsFound = [i, j];
%     % Verifică vecinii 4-dir
%     cellsFound = [cellsFound; floodFill(mat, i-1, j)]; %#ok<AGROW>
%     cellsFound = [cellsFound; floodFill(mat, i+1, j)]; %#ok<AGROW>
%     cellsFound = [cellsFound; floodFill(mat, i, j-1)]; %#ok<AGROW>
%     cellsFound = [cellsFound; floodFill(mat, i, j+1)]; %#ok<AGROW>
% end
% 
% function nodes = getAdjacentNodes(nodeMap, i, j)
%     % Returnează nodurile adiacente pentru celula (i,j) folosind nodeMap
%     [rows, cols] = size(nodeMap);
%     nodes = [];
%     if i > 1, nodes(end+1) = nodeMap(i-1, j); end %#ok<AGROW>
%     if i < rows, nodes(end+1) = nodeMap(i+1, j); end %#ok<AGROW>
%     if j > 1, nodes(end+1) = nodeMap(i, j-1); end %#ok<AGROW>
%     if j < cols, nodes(end+1) = nodeMap(i, j+1); end %#ok<AGROW>
%     % Elimină zerouri și duplicate
%     nodes = unique(nodes(nodes > 0));
% end
