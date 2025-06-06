function [A,B,C,D,state_syms,input_syms] = kirchhoff_to_statespace( ...
    khirchoffs, comps, state_vars, input_vars, output_vars)
% kirchhoff_to_statespace_linear_eliminate  
%   Generează matricile A,B,C,D pentru un circuit liniar R-L-C, 
%   eliminând necunoscutele algebrice prin rezolvare liniară (fără solve).
%
% Intrări:
%   khirchoffs   - struct cu două câmpuri:
%                    KCL: cell array de şiruri (ecuaţii KCL)
%                    KVL: cell array de şiruri (ecuaţii KVL)
%   comps        - cell array de structuri cu câmpurile:
%                    .name  (şir, ex 'R1', 'L2', 'C3', 'Ve' etc.)
%                    .type  ('R','C','L','Vsrc','Isrc')
%                    .param (şir cu parametru simbolic, ex 'R1','L2','C3','Ve')
%   state_vars   - cell array de şiruri cu numele variabilelor de stare, ex {'U_C3','I_L2'}
%   input_vars   - cell array de şiruri cu numele intrărilor, ex {'Ve','Ie'}. 
%   output_vars  - cell array de şiruri cu numele ieșirilor, ex {'Vout'} sau {'U_C3'}.
%
% Ieșiri:
%   A,B           - matricile state-space din ẋ = A*x + B*u
%   C,D           - matricile y   = C*x + D*u
%   state_syms    - vector simbolic al stărilor ca funcţii de timp [U_C3(t); I_L2(t); …]
%   input_syms    - vector simbolic al intrărilor ca funcţii de timp [Ve(t); Ie(t); …]
%
% Schema algoritmului:
%   1) Transformăm ecuațiile KCL/KVL în ecuații simbolice (str2sym + regexprep).
%   2) Substituim legile componentelor R, L, C, Vsrc, Isrc cu subs. 
%   3) Identificăm necunoscutele algebrice (Z), stările (X), derivările stărilor (Xdot) și intrările (U).
%   4) Folosim equationsToMatrix peste eqs_subs cu [Z; Xdot; X; U] pentru a extrage coeficienții matricelor:
%        [ A_Z   A_Xdot   A_X   A_U ] * [Z; Xdot; X; U] = 0. 
%   5) Din prima coloană (corespunzătoare lui Z) extragem A_Z și restul blocurilor. 
%      Rezultă A_alg*Z + B_dotX*Xdot + B_X*X + B_U*U = 0.  
%      Ca atare Z = −A_alg^{-1} * (B_dotX*Xdot + B_X*X + B_U*U).  
%   6) Substituim această expresie a lui Z înapoi în eqs_subs (folosind subs), rămânând
%      doar cu ecuaţii diferențiale în X, Xdot și U:  0 =  M_mat*Xdot + N_mat*X + P_mat*U.  
%   7) Extragem, prin coeffs, M_mat, N_mat, P_mat, apoi calculăm A = −M^{-1}*N și B = −M^{-1}*P.
%   8) Generăm C,D din output_vars cum am făcut anterior (coeficienţi liniari).
%
% Atenție:
%  - În comps NU pune elemente de tip „Vout”! Pentru ieșiri folosește doar output_vars.
%  - state_vars trebuie să conţină EXACT toate U_Cx și I_Lx din circuit.
%  - input_vars trebuie să conţină EXACT toate sursele (Ve, Ie etc.), fără prefixul 'U_' sau 'I_'.

%% 1. Declarații inițiale și simboluri

syms t;  % variabila de timp

n_states = numel(state_vars);
n_inputs = numel(input_vars);

% 1.1. Definim stările X(k) = U_Cx(t) sau I_Lx(t)
X = sym(zeros(n_states,1));
for k = 1:n_states
    X(k) = str2sym([ state_vars{k} '(t)' ]);
end

% 1.2. Definim intrările U(j) = Ve(t) sau Ie(t)
U = sym(zeros(n_inputs,1));
for j = 1:n_inputs
    U(j) = str2sym([ input_vars{j} '(t)' ]);
end

% 1.3. Definim parametrii ca simbolic real (R1, L2, C3, Ve etc.)
param_list = cellfun(@(c) c.param, comps, 'UniformOutput', false);
for idx = 1:numel(comps)
    sym(param_list{idx}, 'real');
end

%% 2. Construim ecuațiile Kirchhoff simbolice (KCL + KVL)

all_eqs_str = [khirchoffs.KCL(:); khirchoffs.KVL(:)];
N_eq = numel(all_eqs_str);
eqs_sym = sym(zeros(N_eq,1));

for i = 1:N_eq
    s = all_eqs_str{i};
    tokens = regexp(s, '([A-Za-z]\w*)', 'tokens');
    unique_vars = unique([tokens{:}]);
    for v = unique_vars
        varname = v{1};
        if ~any(strcmp(param_list, varname))
            % Dacă varname nu e parametru, îl transformăm în funcție de timp:
            s = regexprep(s, ['\<', varname, '\>'], [varname '(t)']);
        end
    end
    eqs_sym(i) = str2sym(s);
end

%% 3. Substituim legile constitutive ale componentelor

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
            subs_to  (end+1) = str2sym([cmp.param '(t)']);

        case 'ISRC'
            % I_Ix(t) = Ie(t)
            Ic = str2sym(['I_' cmp.name '(t)']);
            subs_from(end+1) = Ic;
            subs_to  (end+1) = str2sym([cmp.param '(t)']);

        otherwise
            error('Tip componentă necunoscut: %s', cmp.type);
    end
end

eqs_subs = subs(eqs_sym, subs_from, subs_to);
% Acum eqs_subs conține doar stări (cu var diff(...) pentru C şi L),
% intrări (Ve(t), Ie(t)) şi necunoscute algebrice (I_Rx(t), U_Vout(t) etc.).

%% 4. Identificăm necunoscutele algebrice (Z), stările (X), derivatele (Xdot) și intrările (U)

all_syms = symvar(eqs_subs);  % toți simbolii din ecuații

% 4.1. Stările și derivările lor:
state_syms = X;                              % [U_C1(t); I_L2(t); …]
Xdot_syms  = arrayfun(@(x) diff(x,t), X);    % [diff(U_C1(t),t); diff(I_L2(t),t); …]

% 4.2. Intrările simbolice:
input_syms = U;                              % [Ve(t); Ie(t); …]

% 4.3. Construim lista tuturor „variabilelor cunoscute”:
known_syms = [ state_syms; Xdot_syms; input_syms ];

% 4.4. Necunoscutele algebrice Z = all_syms \ known_syms
algebraic_syms = setdiff(all_syms(:), known_syms(:));

% Dacă nu există necunoscute algebrice, Z este vid:
n_alg = numel(algebraic_syms);
if n_alg == 0
    Z = sym.empty;
else
    Z = algebraic_syms;  % vector simbolic cu necunoscutele algebrice ex. [I_Rr(t); U_V1(t); …]
end

%% 5. Rezolvăm liniar pentru Z

% 5.1. Construim lista completă de necunoscute pentru equationsToMatrix
%   Ordinea: [ Z; Xdot; X; U ]
all_vars_for_ETM = [Z; Xdot_syms; X; U];

% 5.2. Extragem cu equationsToMatrix matricea A_full și b_full, 
%      unde:  A_full * [Z; Xdot; X; U] + b_full = 0.
[A_full, b_full] = equationsToMatrix(eqs_subs, all_vars_for_ETM);

% 5.3. Împărțim A_full în blocuri:
%      A_full = [ A_Z      A_Xdot      A_X      A_U ] 
%   Unde A_Z este coloanele corespunzătoare lui Z, 
%         A_Xdot coloanele pentru derivările stărilor, 
%         A_X coloanele pentru stări, 
%         A_U coloanele pentru intrări.
n_dotx = numel(Xdot_syms);
n_x    = n_states;
n_u    = n_inputs;

% Indici:
idx_Z    = 1:n_alg;
idx_dotX = n_alg + (1:n_dotx);
idx_X    = n_alg + n_dotx + (1:n_x);
idx_U    = n_alg + n_dotx + n_x   + (1:n_u);

A_Z    = A_full(:, idx_Z);
A_dotX = A_full(:, idx_dotX);
A_X    = A_full(:, idx_X);
A_U    = A_full(:, idx_U);

% 5.4. Rezolvăm A_Z * Z + A_dotX * Xdot + A_X * X + A_U * U = 0  
%      ⇒ Z = - A_Z^{-1} * (A_dotX*Xdot + A_X*X + A_U*U)
if n_alg > 0
    if rank(A_Z) < n_alg
        error('A_Z nu este inversibilă (%d × %d), nu pot elimina simbolic necunoscutele algebrice.', n_alg, n_alg);
    end
    % Calculăm în mod simbolic Z în funcție de (Xdot, X, U):
    Z_expr = simplify(- inv(A_Z) * (A_dotX*Xdot_syms + A_X*X + A_U*U));
    %  Z_expr este un vector cu n_alg rânduri. 
    %  Z(i) = Z_expr(i), pentru fiecare necunoscut algebric.
else
    Z_expr = sym.empty; 
end

%% 6. Substituim expresiile lui Z înapoi în ecuațiile inițiale

if n_alg > 0
    subs_from2 = Z;
    subs_to2   = Z_expr;
    eqs_diff   = subs(eqs_subs, subs_from2, subs_to2);
else
    eqs_diff = eqs_subs;
end
% Acum eqs_diff conține numai stări, derivări de stări și intrări:  
%     0 = f( Xdot, X, U ).

%% 7. Extragem M_mat, N_mat, P_mat din eqs_diff

% 7.1. Inițializăm matricile cu size = (#ecuații) × (#stări / #intrări)
M_mat = sym(zeros(N_eq, n_states));  % coef. în faţa diff(X,t)
N_mat = sym(zeros(N_eq, n_states));  % coef. în faţa X
P_mat = sym(zeros(N_eq, n_inputs));  % coef. în faţa U

for i = 1:N_eq
    expr = eqs_diff(i);  % ecuația i, de forma M_i*Xdot + N_i*X + P_i*U = 0

    % 7.a. Coef. pentru diff(X(k),t)
    for k = 1:n_states
        dvar = Xdot_syms(k);
        [Cterms, Tterms] = coeffs(expr, dvar);
        idx_match = find(Tterms == dvar, 1);
        if isempty(idx_match)
            M_mat(i,k) = 0;
        else
            M_mat(i,k) = Cterms(idx_match);
        end
    end

    % 7.b. Coef. pentru X(k)
    for k = 1:n_states
        var = X(k);
        [Cterms, Tterms] = coeffs(expr, var);
        idx_match = find(Tterms == var, 1);
        if isempty(idx_match)
            N_mat(i,k) = 0;
        else
            N_mat(i,k) = Cterms(idx_match);
        end
    end

    % 7.c. Coef. pentru U(j)
    for j = 1:n_inputs
        varu = U(j);
        [Cterms, Tterms] = coeffs(expr, varu);
        idx_match = find(Tterms == varu, 1);
        if isempty(idx_match)
            P_mat(i,j) = 0;
        else
            P_mat(i,j) = Cterms(idx_match);
        end
    end
end

%% 8. Construim A și B

if rank(M_mat) < n_states
    error(['Matricea M_mat nu este inversibilă (rang < %d). Verifică:\n' ...
           '  - ai scris suficient de multe ecuații KCL/KVL?\n' ...
           '  - ai inclus toate stările (U_Cx, I_Lx) în state_vars?'], n_states);
end

A = simplify(- inv(M_mat) * N_mat);
B = simplify(- inv(M_mat) * P_mat);

%% 9. Construim C și D din output_vars

if isempty(output_vars)
    C = sym.empty;
    D = sym.empty;
    state_syms = X;
    input_syms = U;
    return;
end

n_outputs = numel(output_vars);
C = sym(zeros(n_outputs, n_states));
D = sym(zeros(n_outputs, n_inputs));

for k = 1:n_outputs
    yk = output_vars{k};

    % 9.a. Dacă yk e un singur nume (fără + - * /), interpretăm că este „nodul” 
    %      de tensiune la iesire → yk_sym = 'U_<yk>(t)'.
    if isempty(regexp(yk, '[\+\-\*/]', 'once'))
        yk_sym = str2sym(['U_' yk '(t)']);
    else
        % Dacă e o combinație liniară, adăugăm '(t)' lângă fiecare stare/intrare
        expr = yk;
        tokens = regexp(expr, '([A-Za-z]\w*)', 'tokens');
        unique_vars = unique([tokens{:}]);
        for v = unique_vars
            varname = v{1};
            if any(strcmp(state_vars, varname))
                expr = regexprep(expr, ['\<', varname, '\>'], [varname '(t)']);
            elseif any(strcmp(input_vars, varname))
                expr = regexprep(expr, ['\<', varname, '\>'], [varname '(t)']);
            end
        end
        yk_sym = str2sym(expr);
    end

    % 9.b. Extragem coeficienții liniari din yk_sym = C_row * X + D_row * U
    for i = 1:n_states
        var = X(i);
        [Cterms, Tterms] = coeffs(yk_sym, var);
        idx_match = find(Tterms == var, 1);
        if isempty(idx_match)
            C(k, i) = 0;
        else
            C(k, i) = Cterms(idx_match);
        end
    end
    for j = 1:n_inputs
        varu = U(j);
        [Cterms, Tterms] = coeffs(yk_sym, varu);
        idx_match = find(Tterms == varu, 1);
        if isempty(idx_match)
            D(k, j) = 0;
        else
            D(k, j) = Cterms(idx_match);
        end
    end
end

C = simplify(C);
D = simplify(D);

state_syms = X;
input_syms = U;
end
