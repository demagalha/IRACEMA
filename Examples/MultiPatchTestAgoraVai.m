clear all
close all
clc

load RectangularLeissa15.mat
load RectangularLeissa15_2.mat
tic

Model.DegreeElevate(6,1);
Model.DegreeElevate(6,2);
Model.DegreeElevate(6,3);

Model2.DegreeElevate(6,1);
Model2.DegreeElevate(6,2);
Model2.DegreeElevate(6,3);


% Model.KnotRefine(0.05:0.05:0.95,1);
% Model.KnotRefine(0.05:0.05:0.95,2);
% 
% Model2.KnotRefine(0.05:0.05:0.95,1);
% Model2.KnotRefine(0.05:0.05:0.95,2);

YOUNG_MODULUS = 30*10^9;
RHO = 2.32*10^3;
POISSON = 0.2;
D = get_matprop_matrix(1,YOUNG_MODULUS,POISSON); 
[K, M, IEN] = MultiPatchAssemble(Model,Model2,D,RHO);
[V,W] = eigs(K,M,16,'sm');
omega = sqrt(diag(W));
toc