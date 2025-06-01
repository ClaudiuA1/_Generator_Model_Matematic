%CalCulator_ModelMatematic
clc 
clear
close all

mat = createCircuitMatrix()
%[KCL,KVL]=generate_kirchhoff_equations(mat)
% plot_circuit_graph(mat)


  khirchoffs=extractKirchhoff(mat)


    



