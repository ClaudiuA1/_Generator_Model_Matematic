function KVL = extractKVL_general(mat)
% extractKVL_fromMatrix - Extrage ecuațiile KVL direct dintr-o matrice de tip cell M×N,
% unde:
%   '0'     => celulă liberă
%   '1'     => fir de conexiune
%   alt text (ex. 'R1','C1','Ve','E1') => componentă de circuit
%
% Output:
%   KVL  : cell array de string-uri, fiecare reprezentând o ecuație KVL a unei bucle simple

     [M, N] = size(mat);

    % ─── Înainte de a defini dfs, initializează aici stiva și matricile de stare ───
    visitedWire = false(M, N);   % marchează «1»-urile procesate
    onStack     = false(M, N);   % marchează «1»-urile aflate în stiva curentă
    stack       = zeros(0, 2);   % stivă de coordonate [r, c]

    KVL = {};  % aici vom aduna ecuațiile KVL

    % Direcții de explorare la intersecție: sus → dreapta → jos → stânga
    drow = [-1,  0, +1,  0];
    dcol = [ 0, +1,  0, -1];

    % Parcurgem fiecare poziție din matrice
    for i0 = 1:M
        for j0 = 1:N
            if strcmp(mat{i0,j0}, '1') && ~visitedWire(i0,j0)
                % Pornim DFS de la poziția (i0,j0) care conține '1' nevizitat
                dfs(i0, j0);
            end
        end
    end

    % Afișăm ecuațiile KVL extrase
    disp('=== Ecuații KVL extrase din matrice ===');
    disp(KVL(:));

end
function dfs(r, c)
        % Adăugăm (r,c) în stivă
        stack(end+1, :) = [r, c];
        onStack(r, c)   = true;

        % Parcurgem în cele patru direcții, în ordinea sus→dreapta→jos→stânga
        for dir = 1:4
            rr = r + drow(dir);
            cc = c + dcol(dir);

            % Verificăm limitele matricei
            if rr < 1 || rr > M || cc < 1 || cc > N
                continue;
            end
            if ~strcmp(mat{rr,cc}, '1')
                continue;  % nu e fir, nu putem merge pe aici
            end

            if ~onStack(rr, cc)
                % Dacă încă nu e în stivă (nu am revenit înapoi), continuăm DFS recursiv
                dfs(rr, cc);
            else
                % Dacă (rr,cc) e deja în stivă, am găsit un ciclu
                idxInStack = find( (stack(:,1)==rr) & (stack(:,2)==cc), 1 );
                if ~isempty(idxInStack)
                    % Extragem sub-stiva care formează ciclul
                    cycleCells = stack(idxInStack:end, :);  % vector K×2

                    % Construim ecuația KVL pentru acest ciclu
                    kvlEq = buildKVL_fromCycle(cycleCells);

                    % Adăugăm doar dacă nu apare deja
                    if ~any(strcmp(KVL, kvlEq))
                        KVL{end+1} = kvlEq;
                    end
                end
            end
        end

        % Am terminat cu (r,c): îl scoatem din stivă și îl marcăm vizitat
        onStack(r,c)     = false;
        visitedWire(r,c) = true;
        stack(end, :)    = [];  % „pop” din stivă
end


function dfs(r, c)
        % Adăugăm (r,c) în stivă
        stack(end+1, :) = [r, c];
        onStack(r, c)   = true;

        % Parcurgem în cele patru direcții, în ordinea sus→dreapta→jos→stânga
        for dir = 1:4
            rr = r + drow(dir);
            cc = c + dcol(dir);

            % Verificăm limitele matricei
            if rr < 1 || rr > M || cc < 1 || cc > N
                continue;
            end
            if ~strcmp(mat{rr,cc}, '1')
                continue;  % nu e fir, nu putem merge pe aici
            end

            if ~onStack(rr, cc)
                % Dacă încă nu e în stivă (nu am revenit înapoi), continuăm DFS recursiv
                dfs(rr, cc);
            else
                % Dacă (rr,cc) e deja în stivă, am găsit un ciclu
                idxInStack = find( (stack(:,1)==rr) & (stack(:,2)==cc), 1 );
                if ~isempty(idxInStack)
                    % Extragem sub-stiva care formează ciclul
                    cycleCells = stack(idxInStack:end, :);  % vector K×2

                    % Construim ecuația KVL pentru acest ciclu
                    kvlEq = buildKVL_fromCycle(cycleCells);

                    % Adăugăm doar dacă nu apare deja
                    if ~any(strcmp(KVL, kvlEq))
                        KVL{end+1} = kvlEq;
                    end
                end
            end
        end

        % Am terminat cu (r,c): îl scoatem din stivă și îl marcăm vizitat
        onStack(r,c)     = false;
        visitedWire(r,c) = true;
        stack(end, :)    = [];  % „pop” din stivă
end


    function eq = buildKVL_fromCycle(cycleCells)
        % cycleCells: K×2, lista de coordonate [r,c] a buclei detectate (firuri).
        % Vom găsi toate componentele adiacente acestei bucle și vom genera
        % ecuația de tip KVL: „+Uxxx + Uyyy - Uzzz = 0”.

        K = size(cycleCells, 1);

        encountered = {};  % etichete componente unice (fără duplicare)

        % Pentru fiecare muchie din ciclu: (r,c)→(r2,c2)
        for t = 1:(K-1)
            r  = cycleCells(t,   1);
            c  = cycleCells(t,   2);
            r2 = cycleCells(t+1, 1);
            c2 = cycleCells(t+1, 2);

            % Definim cei 8 vecini ai segmentului cele două celule
            neighPositions = [ r-1, c;   r+1, c;   r, c-1;   r, c+1;   ...
                               r2-1, c2; r2+1, c2; r2, c2-1; r2, c2+1 ];

            for p = 1:size(neighPositions,1)
                rr = neighPositions(p,1);
                cc = neighPositions(p,2);
                if rr < 1 || rr > M || cc < 1 || cc > N
                    continue;
                end
                label = mat{rr,cc};
                if ~strcmp(label, '0') && ~strcmp(label, '1')
                    % Am găsit o componentă (ex. 'R1','C2','Ve','E1' etc.)
                    if ~any(strcmp(encountered, label))
                        encountered{end+1} = label;  %#ok<AGROW>
                    end
                end
            end
        end

        % Construim lista cu termeni pentru KVL
        termsList = {};
        for idx = 1:numel(encountered)
            lbl = encountered{idx};
            if startsWith(lbl, 'Ve') || startsWith(lbl, 'E')
                % Sursele de tensiune le marcăm cu semnul „-”
                termsList{end+1} = ['-U' lbl];  %#ok<AGROW>
            else
                % Pentru celelalte componente folosim semnul „+”
                termsList{end+1} = ['+U' lbl];  %#ok<AGROW>
            end
        end

        if isempty(termsList)
            eq = '';  % nu ar trebui să fie cazul, o buclă are măcar o componentă
        else
            eq = [strjoin(termsList, ' + ') ' = 0'];
        end
    end
    % ========================= sfârșit buildKVL_fromCycle =================

