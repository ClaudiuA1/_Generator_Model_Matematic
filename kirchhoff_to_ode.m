function [A,B,C,D,state_syms,input_syms] = kirchhoff_to_ode( ...
    khirchoffs, comps, state_vars, input_vars, output_vars)
% kirchhoff_to_ode
%   Transformă orice set de ecuații KCL/KVL R-L-C direct în ẋ = A x + B u,
%   y = C x + D u, prin solve. Nu includem în Z decât curenții prin rezistori.
%
% Intrări:
%   khirchoffs.KCL = {...}, khirchoffs.KVL = {...}
%   comps = { struct('name','Rr','type','R','param','Rr'), ... }
%   state_vars  = {'U_Cr','I_Ir', …}   (voltaje pe C și curenți prin L)
%   input_vars  = {'Ve', …}            (surse de tensiune/curent)
%   output_vars = {'V1', …}            (noduri de ieșire, de ex. 'V1')
%
% Ieșiri:
%   A,B          – matricele state-space
%   C,D          – matricele de ieșire
%   state_syms   – [U_Cr(t); I_Ir(t); …]
%   input_syms   – [Ve(t); …]

%% 0. Declarații simbolice de bază
syms t;
n_states = numel(state_vars);
n_inputs = numel(input_vars);

% 0.1. Vector simbolic al stărilor X = [ U_Cr(t); I_Ir(t); … ]
X = sym(zeros(n_states,1));
for k = 1:n_states
    X(k) = str2sym([ state_vars{k} '(t)' ]);
end

% 0.2. Vector simbolic al intrărilor U = [ Ve(t); … ]
U = sym(zeros(n_inputs,1));
for j = 1:n_inputs
    U(j) = str2sym([ input_vars{j} '(t)' ]);
end

% 0.3. Definim simbolic parametrii ca reali: Rr, Crr, Lr, Ve, etc.
param_list = cellfun(@(c) c.param, comps, 'UniformOutput', false);
for idx = 1:numel(comps)
    sym(param_list{idx}, 'real');
end

%% 1. Construim ecuațiile Kirchhoff ca simbolice

all_eqs_str = [khirchoffs.KCL(:); khirchoffs.KVL(:)];
N_eq = numel(all_eqs_str);
eqs_sym = sym(zeros(N_eq,1));

for i = 1:N_eq
    s = all_eqs_str{i};
    tokens = regexp(s, '([A-Za-z]\w*)', 'tokens');
    uniq   = unique([tokens{:}]);
    for v = uniq
        varname = v{1};
        if ~any(strcmp(param_list, varname))
            s = regexprep(s, ['\<', varname, '\>'], [varname '(t)']);
        end
    end
    eqs_sym(i) = str2sym(s);
end

%% 2. Substituim legile constitutive (R, C, L, Vsrc, Isrc)

subs_from = sym.empty;
subs_to   = sym.empty;
for idx = 1:numel(comps)
    cmp = comps{idx};
    switch upper(cmp.type)
        case 'R'
            % U_Rx(t) = R * I_Rx(t)
            Uc = str2sym(['U_' cmp.name '(t)']);
            Ic = str2sym(['I_' cmp.name '(t)']);
            R  = sym(cmp.param, 'real');
            subs_from(end+1) = Uc;
            subs_to  (end+1) = R * Ic;
        case 'C'
            % I_Cx(t) = C * diff(U_Cx(t),t)
            Ic = str2sym(['I_' cmp.name '(t)']);
            Uc = str2sym(['U_' cmp.name '(t)']);
            Cc = sym(cmp.param, 'real');
            subs_from(end+1) = Ic;
            subs_to  (end+1) = Cc * diff(Uc, t);
        case 'L'
            % U_Lx(t) = L * diff(I_Lx(t),t)
            Uc = str2sym(['U_' cmp.name '(t)']);
            Ic = str2sym(['I_' cmp.name '(t)']);
            Ll = sym(cmp.param, 'real');
            subs_from(end+1) = Uc;
            subs_to  (end+1) = Ll * diff(Ic, t);
        case 'VSRC'
            % U_Vx(t) = Ve(t)
            Uc = str2sym(['U_' cmp.name '(t)']);
            subs_from(end+1) = Uc;
            subs_to  (end+1) = str2sym([cmp.param '(t)']);  % Ve(t)
        case 'ISRC'
            % I_Ix(t) = Ie(t)
            Ic = str2sym(['I_' cmp.name '(t)']);
            subs_from(end+1) = Ic;
            subs_to  (end+1) = str2sym([cmp.param '(t)']);  % Ie(t)
        otherwise
            error('Tip componentă necunoscut: %s', cmp.type);
    end
end

eqs_subs = subs(eqs_sym, subs_from, subs_to);

%% 3. Construim lista EXPLICITĂ a necunoscutelor algebrice Z

% 3.1. Pentru fiecare rezistor, necunoscuta algebrică e I_<name>(t)
Z = sym.empty;
for idx = 1:numel(comps)
    cmp = comps{idx};
    if strcmpi(cmp.type, 'R')
        Z(end+1) = str2sym(['I_' cmp.name '(t)']);
    end
end

% 3.2. Nu includem tensiunea de ieșire (“U_V1(t)”), 
%      pentru că tocmai ecuația „U_Cr(t) - U_V1(t)=0” leagă direct U_V1 de U_Cr.
%      Așadar, Z conține **doar** curenții prin rezistori.
n_alg = numel(Z);

%% 4. Pregătim derivările stărilor și apelul la solve

% 4.1. Derivatele stărilor: Xdot = [ diff(x1(t),t); diff(x2(t),t); … ]
Xdot_syms = arrayfun(@(x) diff(x, t), X);

% 4.2. Vectorul complet de necunoscute: [Xdot; Z]
unknowns = [ Xdot_syms; Z(:) ];

% 4.3. Rezolvăm simultan pentru derivări și curenți prin rezistori
try
    sol = solve(eqs_subs, unknowns, 'IgnoreAnalyticConstraints', true);
catch ME
    disp('--- Eroare în rezolvare simbolică ---');
    disp('Ecuațiile după substituții (eqs_subs):');
    disp(eqs_subs);
    disp('Necunoscutele cerute la solve (unknowns):');
    disp(unknowns);
    error('Nu am reușit să rezolv simbolic. Verifică state_vars, KCL/KVL și compușii.');
end

%% 5. Extragem expresiile pentru derivările stărilor dx = [dot(x1); dot(x2); …]

dx = sym(zeros(n_states,1));
for k = 1:n_states
    fld = char(Xdot_syms(k));       % ex: 'D(U_Cr(t),t)'
    dx(k) = simplify(sol.(fld));    % sol.D(U_Cr(t),t)
end
% Acum dx(k) e de forma f_k(x,u), fără I_Rr(t) în final.

%% 6. Construim matricile A și B din dx

A = sym(zeros(n_states, n_states));
B = sym(zeros(n_states, n_inputs));

for k = 1:n_states
    expr = dx(k);  % lineară în X(j)(t) și U(m)(t)
    % coeficient pe stări X(j)
    for j = 1:n_states
        [Cterms, Tterms] = coeffs(expr, X(j));
        idx = find(Tterms == X(j), 1);
        if isempty(idx)
            A(k,j) = 0;
        else
            A(k,j) = Cterms(idx);
        end
    end
    % coeficient pe intrări U(m)
    for m = 1:n_inputs
        [Cterms, Tterms] = coeffs(expr, U(m));
        idx = find(Tterms == U(m), 1);
        if isempty(idx)
            B(k,m) = 0;
        else
            B(k,m) = Cterms(idx);
        end
    end
end

A = simplify(A);
B = simplify(B);

%% 7. Construim C și D din output_vars

if isempty(output_vars)
    C = sym.empty;
    D = sym.empty;
    state_syms = X;
    input_syms = U;
    return;
end

n_out = numel(output_vars);
C = sym(zeros(n_out, n_states));
D = sym(zeros(n_out, n_inputs));

for k = 1:n_out
    yk = output_vars{k};  % ex 'V1'
    % Dacă yk e simplu (fără + - * /), deci e un nod de ieșire:
    if isempty(regexp(yk, '[\+\-\*/]', 'once'))
        % y = U_V1(t)
        y_sym = str2sym(['U_' yk '(t)']);
    else
        % yk e o expresie liniară în stări/intrări
        expr = yk;
        tokens = regexp(expr, '([A-Za-z]\w*)', 'tokens');
        uniq   = unique([tokens{:}]);
        for v = uniq
            nm = v{1};
            if any(strcmp(state_vars, nm))
                expr = regexprep(expr, ['\<',nm,'\>'], [nm '(t)']);
            elseif any(strcmp(input_vars, nm))
                expr = regexprep(expr, ['\<',nm,'\>'], [nm '(t)']);
            end
        end
        y_sym = str2sym(expr);
    end

    for j = 1:n_states
        [Cterms, Tterms] = coeffs(y_sym, X(j));
        idx = find(Tterms == X(j), 1);
        if isempty(idx)
            C(k,j) = 0;
        else
            C(k,j) = Cterms(idx);
        end
    end
    for m = 1:n_inputs
        [Cterms, Tterms] = coeffs(y_sym, U(m));
        idx = find(Tterms == U(m), 1);
        if isempty(idx)
            D(k,m) = 0;
        else
            D(k,m) = Cterms(idx);
        end
    end
end

C = simplify(C);
D = simplify(D);

state_syms = X;
input_syms = U;
end
