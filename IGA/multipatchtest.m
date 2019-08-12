clear all
close all
clc

load RectangularLeissa15.mat
load RectangularLeissa15_2.mat

%REFINEMENT
%{
Model.DegreeElevate(2,1);
Model.DegreeElevate(3,2);
Model.DegreeElevate(1,3);

Model2.DegreeElevate(2,1);
Model2.DegreeElevate(3,2);
Model2.DegreeElevate(1,3);
%}

Model.KnotRefine(0.05:0.05:0.95,1);
Model.KnotRefine(0.05:0.05:0.95,2);

Model2.KnotRefine(0.05:0.05:0.95,1);
Model2.KnotRefine(0.05:0.05:0.95,2);

%END OF REFINEMENT

%CONTROL POINTS
P1 = Model.get_point_cell;
P2 = Model2.get_point_cell;
%END OF CONTROL POINTS

%PATCH 1 MATERIAL PROPERTIES
YOUNG_MODULUS_1 = 30*10^9;
RHO_1 = 2.32*10^3;
POISSON_1 = 0.2;
D_1 = get_matprop_matrix(1,YOUNG_MODULUS_1,POISSON_1); 
%END OF PATCH 1 MATERIAL PROPERTIES

%PATCH 2 MATERIAL PROPERTIES
YOUNG_MODULUS_2 = 30*10^9;
RHO_2 = 2.32*10^3;
POISSON_2 = 0.2;
D_2 = get_matprop_matrix(1,YOUNG_MODULUS_2,POISSON_2);
%END OF PATCH 2 MATERIAL PROPERTIES

%PATCH 1 ASSEMBLY
[K_1,M_1,IEN_1] = Assemble(Model,D_1,RHO_1);
ID_1 = reshape(1:max(max(IEN_1))*3,3,max(max(IEN_1)));
%END OF PATCH 1 ASSEMBLY

%PATCH 2 ASSEMBLY
[K_2,M_2,IEN_2] = Assemble(Model2,D_2,RHO_2);
ID_2 = reshape(1:max(max(IEN_2))*3,3,max(max(IEN_2)));
%END OF PATCH 2 ASSEMBLY

%MASTER-SLAVE INDEX
%get the control points index to the common surface on each patch. (MS(1,:)
%returns [a b], being a the index of control point on solid 1 and the
%equivalent b for the solid 2

MS = MasterSlave(Model,Model2,6,1);

%END OF MASTER-SLAVE INDEX

ID_M = reshape(ID_1(:,MS(:,1)),numel(ID_1(:,MS(:,1))),1);
ID_S = reshape(ID_2(:,MS(:,2)),numel(ID_2(:,MS(:,2))),1);
OVERLAP_SIZE = numel(ID_M);

[sz1, ~] = size(K_1);
[sz2, ~] = size(K_2);
sz3 = sz1 - OVERLAP_SIZE;
sz4 = sz2 - OVERLAP_SIZE;
K = zeros(sz1+sz2-OVERLAP_SIZE);
M = K;
K(1:sz1,1:sz1) = K_1;
M(1:sz1,1:sz1) = M_1;
K(sz3+1:end,sz3+1:end) = K(sz3+1:end,sz3+1:end) +K_2;
M(sz3+1:end,sz3+1:end) = M(sz3+1:end,sz3+1:end) +M_2;
for i=numel(ID_M)
    mrow = ID_M(i);
    srow = ID_S(i);
    [~,mcol,m] = find(K_1(mrow,:)); % Rows and Columns of M DOF
    [~,scol,s] = find(K_2(srow,:)); % Rows and Columns of S DOF
    K(mrow,mcol) = K(mrow,mcol) -K_2(srow,scol);    
    M(mrow,mcol) = M(mrow,mcol) -M_2(srow,scol);
end

K = sparse(K);
M = sparse(M);
[V,W] = eigs(K,M,16,'sm');
omega = sqrt(diag(W));

