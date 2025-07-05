function [state_vars, comps] = getStateVars(kirchoffs)
% getStateVars - Versiune robustă care extrage componentele și apoi
% determină variabilele de stare pe baza listei de componente.

    % --- ETAPA 1: Extrage numele tuturor componentelor din ecuații ---
    all_equations = [kirchoffs.KCL; kirchoffs.KVL];
    all_component_names = {};
    for i = 1:numel(all_equations)
        % Găsește orice nume de componentă care urmează după U_ sau I_
        matches = regexp(all_equations{i}, '(?<=[UI]_)[A-Z]\w*', 'match');
        all_component_names = [all_component_names, matches];
    end
    unique_names = sort(unique(all_component_names));
    
    % Ignorăm sursele de tensiune de ieșire sau alte elemente speciale
    % care nu sunt componente fizice reale în sensul R, L, C, I_src, V_src
    unique_names = unique_names(~startsWith(unique_names, 'V1')); % Exemplu de filtrare

    % --- ETAPA 2: Construiește lista `comps` ȘI identifică stările în același timp ---
    
    % Pre-alocare pentru eficiență
    comps = cell(numel(unique_names), 1);
    capacitor_vars = {};
    inductor_vars = {};

    for i = 1:numel(unique_names)
        name = unique_names{i};
        type = name(1); % Prima literă determină tipul
        
        % Adăugăm componenta în lista `comps`
        component_struct = struct('name', name, 'param', name);
        
        switch type
            case 'R'
                component_struct.type = 'R';
            case 'C'
                component_struct.type = 'C';
                % Este condensator, deci starea este TENSIUNEA
                capacitor_vars{end+1} = ['U_' name]; %#ok<AGROW>
            case 'I'
                component_struct.type = 'L'; % Convenția noastră: 'I' în nume => tip 'L' (bobină)
                % Este bobină, deci starea este CURENTUL
                inductor_vars{end+1} = ['I_' name]; %#ok<AGROW>
            case 'V'
                component_struct.type = 'Vsrc';
            % Poți adăuga și 'I' pentru surse de curent dacă e cazul
        end
        comps{i} = component_struct;
    end

    % --- ETAPA 3: Asamblează lista finală de state_vars ---
    
    % Sortăm și combinăm pentru a obține lista finală ordonată
    state_vars = [sort(unique(capacitor_vars)), sort(unique(inductor_vars))];

end