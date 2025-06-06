function [A,B,C,D,state_syms,input_syms] = kirchhoff_to_statespace_linear_eliminate( ...
    khirchoffs, comps, state_vars, input_vars, output_vars)
% kirchhoff_to_statespace_linear_eliminate  
%   Generează matricile A,B,C,D ale unui circuit R-L-C cu intrare/ieșire,
%   eliminând necunoscutele algebrice (curenți prin rezistori şi tensiuni de ieșire)
%   în mod explicit, fără a se baza pe symvar sau atoms(symfun).
%
% Intrări:
%   khirchoffs   - struct cu:
%                    .KCL (cell array de șiruri, ex {'I_Rr=I_Ir+I_Cr'})
%                    .KVL (cell array de șiruri, ex {'U_Rr+U_Ir-U_Ve=0', 'U_Cr-U_Ir=0', 'U_Cr-U_V1=0'})
%   comps        - cell array de structuri, fiecare cu:
%                    .name  (șir, ex 'Rr','Ir','Cr','Ve' etc.)
%                    .type  ('R','C','L','Vsrc','Isrc')
%                    .param (șir cu numele simbolic, ex 'Rr','Crr','Lr','Ve')
%   state_vars   - cell array de șiruri cu stările, ex {'U_Cr','I_Ir'}
%   input_vars   - cell array de șiruri cu intrările, ex {'Ve','Ie'}
%   output_vars  - cell array de șiruri cu ieșirile, ex {'V1'} sau {'U_Cr'} (fără prefix 'U_')
%
% Ieșiri:
%   A,B          - matricile state-space (ẋ = A x + B u)
%   C,D          - matricile pentru ieșire (y = C x + D u)
%   state_syms   - vector simbolic cu stările (funcții de timp) [U_Cr(t); I_Ir(t); …]
%   input_syms   - vector simbolic cu intrările (funcții de timp) [Ve(t); Ie(t); …]
%
% Algoritm principal:
%   1) Convertim ecuațiile KCL/KVL din șiruri în ecuații simbolice, introducând '(t)' 
%      pentru toate variabilele care nu sunt parametri.
%   2) Substituim legile componentelor R, C, L, Vsrc, Isrc (ex.: U_R = R I_R, I_C = C dU_C/dt, etc.).
%   3) Extragem **explicit** necunoscutele algebrice Z:
%       – Pentru fiecare componentă de tip 'R', necunoscuta este I_<numeR>(t).
%       – Pentru fiecare element din output_vars, necunoscuta este U_<output>(t).
%   4) Pentru fiecare ecuație, extragem coeficienții liniari cu `coeffs(expr, var)` pentru:
%       – fiecare necunoscut algebric Z(k)
%       – fiecare derivată diff(X(j), t) 
%       – fiecare stare X(m)
%       – fiecare intrare U(n)
%      și construim blocurile A_Z, A_dotX, A_X, A_U.
%   5) Rezolvăm liniar A_Z * Z + A_dotX * Xdot + A_X * X + A_U * U = 0  ⟹  
%       Z = – A_Z⁻¹ (A_dotX * Xdot + A_X * X + A_U * U).
%   6) Substituim expresiile pentru Z înapoi în ecuații → obținem ecuații diferențiale
%       0 = M_mat * Xdot + N_mat * X + P_mat * U.
%   7) Extragem M_mat, N_mat, P_mat cu coeffs → A = – M⁻¹ N, B = – M⁻¹ P.
%   8) Construim C, D folosind tot `coeffs`, din output_vars (ex.: y = U_V1(t) = C*X + D*U).
%
% Atenție:
%  - În lista `comps`, **nu** pune «Vout» sau «V1». Ieșirea se declară în `output_vars`.
%  - `state_vars` trebuie să conțină *exact* toate U_Cx și I_Lx din circuit.
%  - `input_vars` trebuie să conțină *exact* toate sursele (Ve, Ie etc.), fără prefix.
%  - `output_vars` conține numele nodurilor de ieșire (fără „U_”), ex: {'V1'}, iar în KVL
%    trebuie să existe o ecuație „U_Cnod - U_V1 = 0” care leagă ieșirea de o stare.

%% 1. Declarații inițiale și definirea stărilor/intrărilor

syms t;  % variabila de timp

n_states = numel(state_vars);
n_inputs = numel(input_vars);

% 1.1. Construim vectorul X al stărilor, ca funcții de timp „*(t)”
X = sym(zeros(n_states,1));
for k = 1:n_states
    X(k) = str2sym([ state_vars{k} '(t)' ]);
end

% 1.2. Construim vectorul U al intrărilor, ca funcții de timp „*(t)”
U = sym(zeros(n_inputs,1));
for j = 1:n_inputs
    U(j) = str2sym([ input_vars{j} '(t)' ]);
end

% 1.3. Definim simbolic parametrii ca fiind reali (Rr, Crr, Lr, Ve etc.)
param_list = cellfun(@(c) c.param, comps, 'UniformOutput', false);
for idx = 1:numel(comps)
    sym(param_list{idx}, 'real');
end

%% 2. Transformăm ecuațiile KCL + KVL în ecuații simbolice

all_eqs_str = [khirchoffs.KCL(:); khirchoffs.KVL(:)];
N_eq = numel(all_eqs_str);
eqs_sym = sym(zeros(N_eq,1));

for i = 1:N_eq
    s = all_eqs_str{i};
    % Găsim toate "cuvintele" (alfanumerice) din șir
    tokens = regexp(s, '([A-Za-z]\w*)', 'tokens');
    unique_vars = unique([tokens{:}]);
    for v = unique_vars
        varname = v{1};
        % Dacă "varname" nu este un parametru în comps, atunci e o necunoscută și
        % trebuie să devină o funcție de timp: varname → varname(t)
        if ~any(strcmp(param_list, varname))
            s = regexprep(s, ['\<', varname, '\>'], [varname '(t)']);
        end
    end
    % Convertim șirul final în ecuație simbolică
    eqs_sym(i) = str2sym(s);
end

%% 3. Substituim legile constitutive (R, L, C, Vsrc, Isrc)

subs_from = sym.empty;
subs_to   = sym.empty;

for idx = 1:numel(comps)
    cmp = comps{idx};
    switch upper(cmp.type)
        case 'R'
            % Rezistor: U_Rx(t) = R * I_Rx(t)
            Uc = str2sym(['U_' cmp.name '(t)']);   % ex: U_Rr(t)
            Ic = str2sym(['I_' cmp.name '(t)']);   % ex: I_Rr(t)
            R  = sym(cmp.param, 'real');           % ex: Rr
            subs_from(end+1) = Uc;
            subs_to  (end+1) = R * Ic;

        case 'C'
            % Condensator: I_Cx(t) = C * diff(U_Cx(t), t)
            Ic = str2sym(['I_' cmp.name '(t)']);
            Uc = str2sym(['U_' cmp.name '(t)']);
            Cc = sym(cmp.param, 'real');  % ex: Crr
            subs_from(end+1) = Ic;
            subs_to  (end+1) = Cc * diff(Uc, t);

        case 'L'
            % Bobină: U_Lx(t) = L * diff(I_Lx(t), t)
            Uc = str2sym(['U_' cmp.name '(t)']);
            Ic = str2sym(['I_' cmp.name '(t)']);
            Ll = sym(cmp.param, 'real');  % ex: Lr
            subs_from(end+1) = Uc;
            subs_to  (end+1) = Ll * diff(Ic, t);

        case 'VSRC'
            % Sursă de tensiune: U_Vx(t) = Ve(t)
            Uc = str2sym(['U_' cmp.name '(t)']);   % ex: U_Ve(t)
            subs_from(end+1) = Uc;
            subs_to  (end+1) = str2sym([cmp.param '(t)']);  % ex: Ve(t)

        case 'ISRC'
            % Sursă de curent: I_Ix(t) = Ie(t)
            Ic = str2sym(['I_' cmp.name '(t)']);   % ex: I_Ie(t)
            subs_from(end+1) = Ic;
            subs_to  (end+1) = str2sym([cmp.param '(t)']);  % ex: Ie(t)

        otherwise
            error('Tip componentă necunoscut: %s', cmp.type);
    end
end

eqs_subs = subs(eqs_sym, subs_from, subs_to);
% --- Acum eqs_subs conține doar:
%     1) Stările (ca funcții de t, unele cu diff(...) pentru C și L),
%     2) Intrările (Ve(t), Ie(t)),
%     3) Necunoscutele algebrice rămase (I_Rx(t) pentru rezistori și U_Vout(t) dacă ieșirea e nod),
%     4) Parametrii (Rr, Crr, Lr, Ve etc. – dar aceștia vor fi tratați mai jos).

%% 4. Construim lista EXPLICITĂ a necunoscutelor algebrice Z

% 4.1. Pentru fiecare componentă de tip R, necunoscuta algebrică este I_<name>(t)
Z_R = sym.empty;
for idx = 1:numel(comps)
    cmp = comps{idx};
    if strcmpi(cmp.type, 'R')
        Z_R(end+1) = str2sym(['I_' cmp.name '(t)']);
    end
end

% 4.2. Pentru fiecare ieșire (output_vars), necunoscuta algebrică e U_<output>(t)
Z_out = sym.empty;
for k = 1:numel(output_vars)
    outname = output_vars{k};  % ex 'V1'
    Z_out(end+1) = str2sym(['U_' outname '(t)']);
end

% 4.3. Combinăm lista completă a lui Z:
Z = [Z_R(:); Z_out(:)];
n_alg = numel(Z);

% 4.4. În mod EXTREM de rar, dacă nu există rezistoare și nici ieșiri, Z poate fi vid:
if n_alg == 0
    % Nu avem necunoscute algebrice
    Z = sym.empty;
end

%% 5. Extragem coeficientii folosind jacobian (metoda robusta)

% 5.1. Construim lista Xdot
Xdot_syms = arrayfun(@(x) diff(x, t), X);

% 5.2. Asiguram ca ecuatiile sunt de forma 'expr = 0'
eqs_expr = lhs(eqs_subs) - rhs(eqs_subs);

% 5.3. Extragem matricele de coeficienti direct cu jacobian
% Daca sistemul e liniar, jacobianul expresiilor fata de variabile
% este chiar matricea coeficientilor.
if n_alg > 0
    A_Z = jacobian(eqs_expr, Z);
else
    A_Z = sym([]); % Matrice goala daca nu avem var. algebrice
end

A_dotX = jacobian(eqs_expr, Xdot_syms);
A_X    = jacobian(eqs_expr, X);
A_U    = jacobian(eqs_expr, U);

% 5.4. Rezolvam liniar pentru Z
disp('Matricea A_Z calculata cu jacobian:');
disp(A_Z);
disp('Necunoscutele algebrice Z:');
disp(Z);

if n_alg > 0
    % Verificam daca A_Z este patratica si inversibila.
    % Numarul de ecuatii trebuie sa fie cel putin egal cu numarul de necunoscute.
    if rank(A_Z) < n_alg
        error('Blocul A_Z (dim %d x %d) nu este inversibil. Verificati KCL/KVL.', size(A_Z,1), size(A_Z,2));
    end
    % LINIA CORECTATĂ (folosind operatorul backslash)
% BLOCUL NOU - Selectează un sub-sistem pătratic și consistent
RHS = - (A_dotX*Xdot_syms + A_X*X + A_U*U);

% BLOCUL NOU - Compatibil cu versiuni mai vechi de MATLAB

% 1. Apelam rref cu un singur argument de iesire
R = rref(A_Z');

% 2. Gasim manual coloanele pivot din matricea R
pivot_rows_idx = [];
[num_rows, ~] = size(R);
for i = 1:num_rows
    % Gasim prima pozitie nenula pe fiecare rand
    pivot_col = find(R(i, :), 1, 'first');
    if ~isempty(pivot_col)
        % Adaugam indexul coloanei la lista noastra de pivoti
        pivot_rows_idx(end+1) = pivot_col;
    end
end
% Sortam pentru a mentine o ordine consistenta, desi nu e strict necesar
pivot_rows_idx = sort(pivot_rows_idx);

if length(pivot_rows_idx) < n_alg
    error('Sistemul algebric este dependent. Nu se pot determina toate necunoscutele Z.');
end

% Selectăm doar rândurile independente pentru a forma un sistem pătratic
A_Z_sq = A_Z(pivot_rows_idx, :);
RHS_sq = RHS(pivot_rows_idx, :);

fprintf('S-a selectat un sub-sistem %d x %d pentru a rezolva Z.\n', size(A_Z_sq,1), size(A_Z_sq,2));

% Rezolvăm sistemul PĂTRATIC, care acum are o soluție unică și finită
Z_expr = simplify(A_Z_sq \ RHS_sq);
else
    Z_expr = sym.empty;
end
%% 6. Substituim expresiile pentru Z în ecuații, obţinem doar ecuații diferențiale

if n_alg > 0
    eqs_diff = subs(eqs_subs, Z, Z_expr);
else
    eqs_diff = eqs_subs;
end
% ⇒ Acum eqs_diff conține doar termeni în Xdot, X și U:
%    0 = f_i(Xdot, X, U)

%% 7. Extragem M_mat, N_mat si P_mat folosind jacobian

% Asiguram ca ecuatiile diferentiale sunt de forma 'expr = 0'
eqs_diff_expr = lhs(eqs_diff) - rhs(eqs_diff);

% Extragem matricele M, N, P
M_mat = jacobian(eqs_diff_expr, Xdot_syms);
N_mat = jacobian(eqs_diff_expr, X);
P_mat = jacobian(eqs_diff_expr, U);

%% 8. Selectam randurile independente si construim A si B

% Combinam matricile pentru a gasi randurile independente ale intregului sistem dinamic
Full_dyn_mat = [M_mat, N_mat, P_mat];

% Folosim metoda compatibila cu rref pentru a gasi randurile pivot
R_dyn = rref(Full_dyn_mat');
pivot_rows_dyn_idx = [];
[num_rows_dyn, ~] = size(R_dyn);
for i = 1:num_rows_dyn
    pivot_col_dyn = find(R_dyn(i, :), 1, 'first');
    if ~isempty(pivot_col_dyn)
        pivot_rows_dyn_idx(end+1) = pivot_col_dyn;
    end
end
pivot_rows_dyn_idx = sort(pivot_rows_dyn_idx);

% Verificam daca am gasit suficiente ecuatii de stare
if length(pivot_rows_dyn_idx) < n_states
    error('Nu s-au putut gasi suficiente ecuatii dinamice independente (%d gasite, %d necesare). Verifica KCL/KVL sau definitia starilor.', ...
        length(pivot_rows_dyn_idx), n_states);
end

% Selectam doar randurile independente pentru a forma matrici patratice/corecte
M_sq = M_mat(pivot_rows_dyn_idx, :);
N_sq = N_mat(pivot_rows_dyn_idx, :);
P_sq = P_mat(pivot_rows_dyn_idx, :);

% Acum M_sq este patratica si inversibila. Calculam A si B.
if rank(M_sq) < n_states
    % Aceasta eroare nu ar trebui sa apara daca logica e corecta,
    % dar o pastram ca masura de siguranta.
    error('Matricea M_mat selectata nu este inversibilă (rang < %d). Verifică state_vars.', n_states);
end

A = simplify(-inv(M_sq) * N_sq);
B = simplify(-inv(M_sq) * P_sq);
%% 9. Construim C și D din output_vars folosind jacobian

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

% Construim un vector simbolic al iesirilor, Y
Y_sym = sym(zeros(n_outputs, 1));
for k = 1:n_outputs
    yk_str = output_vars{k};
    % Aici trebuie sa definim iesirile ca expresii ale starilor, intrarilor, etc.
    % Presupunem ca iesirea este o tensiune de nod, ex: 'V1' -> U_V1(t)
    if isempty(regexp(yk_str, '[\+\-\*/]', 'once'))
        yk_sym = str2sym(['U_' yk_str '(t)']);
    else
        % Daca e o expresie mai complexa, ar trebui tratata similar cu pasul 2
        % (momentan nu este implementat in detaliu aici)
        yk_sym = str2sym(yk_str); % Simplificare
    end
    
    % Substituim necunoscutele algebrice (inclusiv iesirea U_V1(t))
    % cu expresiile lor in functie de X si U
    if n_alg > 0
        Y_sym(k) = subs(yk_sym, Z, Z_expr);
    else
        Y_sym(k) = yk_sym;
    end
end

% Acum Y_sym este exprimat doar in functie de X si U.
% Calculam C si D ca jacobienii lui Y fata de X si U.
C = simplify(jacobian(Y_sym, X));
D = simplify(jacobian(Y_sym, U));

state_syms = X;
input_syms = U;

end
