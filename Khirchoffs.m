function ecuatii = Kirchoff(mat)
    ecuatii = {}; % Inițializăm lista de ecuații

    [rows, cols] = size(mat);

    % Inițializăm ecuațiile pentru tensiune și curent
    tensiune_eq = "";
    curent_eq = "";

    % Parcurgem matricea pentru a găsi componentele circuitului
    for i = 1:rows
        for j = 1:cols
            elem = mat{i, j};
            
            % Verificăm dacă este o componentă (nu 0 sau 1)
            if ~strcmp(elem, "0") && ~strcmp(elem, "1")
                
                % Construim ecuația de tensiune (Kirchhoff 2)
                if i > 1 && ~strcmp(mat{i-1, j}, "0") && ~strcmp(mat{i-1, j}, "1") % Conexiune sus
                    tensiune_eq = strcat(tensiune_eq, "U", elem, " + ");
                end
                if i < rows && ~strcmp(mat{i+1, j}, "0") && ~strcmp(mat{i+1, j}, "1") % Conexiune jos
                    tensiune_eq = strcat(tensiune_eq, "U", elem, " + ");
                end
                if j > 1 && ~strcmp(mat{i, j-1}, "0") && ~strcmp(mat{i, j-1}, "1") % Conexiune stânga
                    tensiune_eq = strcat(tensiune_eq, "U", elem, " + ");
                end
                if j < cols && ~strcmp(mat{i, j+1}, "0") && ~strcmp(mat{i, j+1}, "1") % Conexiune dreapta
                    tensiune_eq = strcat(tensiune_eq, "U", elem, " + ");
                end

                % Construim ecuația de curent (Kirchhoff 1)
                curent_eq = strcat(curent_eq, "i", elem, " = ");
            end
        end
    end

    % Eliminăm ultimul "+" sau "=" din ecuații
    if strlength(tensiune_eq) > 3
        tensiune_eq = extractBefore(tensiune_eq, strlength(tensiune_eq) - 2);
    end
    if strlength(curent_eq) > 3
        curent_eq = extractBefore(curent_eq, strlength(curent_eq) - 2);
    end

    % Adăugăm ecuațiile în listă dacă nu sunt goale
    if ~isempty(tensiune_eq)
        ecuatii{end+1} = tensiune_eq;
    end
    if ~isempty(curent_eq)
        ecuatii{end+1} = curent_eq;
    end

    % Afișăm ecuațiile rezultate
    for e = 1:length(ecuatii)
        disp(ecuatii{e});
    end
end
