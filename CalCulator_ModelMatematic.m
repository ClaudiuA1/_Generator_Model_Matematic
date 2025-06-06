%CalCulator_ModelMatematic
clc 
clear
close all

mat = createCircuitMatrix()
%[KCL,KVL]=generate_kirchhoff_equations(mat)
% plot_circuit_graph(mat)


  % khirchoffs=extractKirchhoff(mat);
% 1) Definirea ecuațiilor Kirchhoff (stringuri)
khirchoffs.KCL = {
    'I_Rr1 = I_C1'
    'I_Rr1 = I_R2'
    'I_C1 = I_R2'
    'I_Rr = I_Ir'
    'I_Ir = I_Cr + I_R1'
    %  Observație: nu adăugăm „+ I_V1”, pentru că tu ai spus că
    %   V1 nu e rezistență, ci doar tensiune de ieșire.
};
khirchoffs.KVL = {
    'U_Cr + U_Rr + U_Ir - U_Ve = 0';  % bucla Rr–Lr–Ve
     'U_R2 - U_V1 = 0';  
     '  U_R2 + U_R1 + U_C1 - U_Cr = 0'   % bucla Lr–Cr
             % bucla Cr–V1 (nodul de ieșire)
};

% 2) Lista componentelor (fără V1 ca 'Vout')
comps = {
    struct('name','Rr','type','R','param','Rr');
    struct('name','R1','type','R','param','R1');
    struct('name','R2','type','R','param','R2');

    struct('name','Ir','type','L','param','Lr');
    struct('name','Cr','type','C','param','Cr');
    struct('name','C1','type','C','param','C1');

    struct('name','Ve','type','Vsrc','param','Ve');
};

% 3) Definirea stărilor și intrărilor
state_vars = {'U_Cr','U_C1','I_Ir'};  % tensiunea pe C și curentul prin L
input_vars = {'Ve'};           % Ve(t) – tensiunea de intrare
output_vars = {'V1'};          % V1(t) – ieșirea = tensiune nod

% 4) Apelul funcției
[A,B,C,D, Xs, Us] = kirchhoff_to_statespace_linear_eliminate( ...
    khirchoffs, comps, state_vars, input_vars, output_vars);

disp('A = '); disp(simplify(A));
disp('B = '); disp(simplify(B));
disp('C = '); disp(simplify(C));
disp('D = '); disp(simplify(D));
