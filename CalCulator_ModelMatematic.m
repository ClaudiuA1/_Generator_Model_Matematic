%CalCulator_ModelMatematic
clc 
clear
close all

mat = createCircuitMatrix()
%   [KCL, KVL]=extractKirchhoff(mat);
 [KCL, KVL]=generateKirchhoffFromMatrix(mat);


  kirchoffs.KCL=KCL;
 kirchoffs.KVL=KVL;

%%
% % 1) Definirea ecuațiilor Kirchhoff (stringuri)
% close all
% clc
% kirchoffs.KCL = { 
% 'I_C0 +I_I0 +I_I1 = 0' 
%     '-I_C0 -I_I1 +I_R1 = 0'
%     'I_I0 -I_R0 = 0'      
%     '-I_I2 -I_C1 +I_R2 = 0'
%     'I_C1 - I_R2 = 0'
% };
% kirchoffs.KVL = {
% 'U_C0 - U_I1 = 0'                        
%     'U_R1 - U_V1 = 0'                        
%     '-U_C1 +U_I2 +U_R2 = 0'                  
%     '-U_V0 +U_R0 -U_I0 +U_C0 +U_R1 -U_I2 = 0'
% 
% 
% };

[state_vars,comps]=getStateVars(kirchoffs);

%  Definirea intrărilor si iesirilor
% Constrangere: V0 va fi intotdeauna intrarea si V1 intotdeauna iesirea
% Daca nu se respecta aceasta constrangere trebui schimbat in sectiunea 3
% de generare KCL, si in getStateVars, filtrul de excludere a V1
 input_vars = {'V0'};           % V0(t) – tensiunea de intrare
 output_vars = {'V1'};          % V1(t) – ieșirea = tensiune nod

% 4) Apelul funcției
 [A,B,C,D] =kirchhoff_to_statespace(kirchoffs, comps, state_vars, input_vars, output_vars)
%%
% disp('A = '); disp(simplify(A));
% disp('B = '); disp(simplify(B));
% disp('C = '); disp(simplify(C));
% disp('D = '); disp(simplify(D));


syms R0 R1 C0 C1 R2 I0
R0_var=1000;
R1_var=500;
R2_var=750;
C0_var=100e-6;
C1_var=10e-6;
I0_var=1;

A=double(subs(A,{R0,R1,R2,C0,C1,I0},{R0_var,R1_var,R2_var,C0_var,C1_var,I0_var}));
B=double(subs(B,{R0,R1,R2,C0,C1,I0},{R0_var,R1_var,R2_var,C0_var,C1_var,I0_var}));
C=double(subs(C,{R0,R1,R2,C0,C1,I0},{R0_var,R1_var,R2_var,C0_var,C1_var,I0_var}));
D=double(D);


syms s
H=tf("s")

Hf= minreal(zpk(C*((H*eye(size(A))-A)^-1)*B+D))

step(feedback(Hf,1))