%clear all
%close all
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

Model.KnotRefine(0.1:0.1:0.9,1);
Model.KnotRefine(0.1:0.1:0.9,2);

Model2.KnotRefine(0.1:0.1:0.9,1);
Model2.KnotRefine(0.1:0.1:0.9,2);

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


%local_1 = reshape(ID_1(:,MS(:,1)),numel(ID_1(:,MS(:,1))),1);
%local_2 = reshape(ID_2(:,MS(:,2)),numel(ID_2(:,MS(:,2))),1);
%%where

for i=1:numel(MS,1)
    ID_M = ID_1(:,MS(i,1));
    ID_S = ID_2(:,MS(i,2));
    
    for ii=1:numel(ID_M)
        for jj=1:numel(ID_S)
            K_1(ID_M(ii),ID_M(jj)) = K_1(ID_M(ii),ID_M(jj)) + K_2(ID_S(ii),ID_S(jj));
            M_1(ID_M(ii),ID_M(jj)) = M_1(ID_M(ii),ID_M(jj)) + M_2(ID_S(ii),ID_S(jj));
        end
    end
end


temp = sort(MS(:,2),'descend');

for i=1:numel(temp)
    K_2(:,temp(i)) = [];
    K_2(temp(i),:) = [];
    
    M_2(:,temp(i)) = [];
    M_2(temp(i),:) = [];
end

K = sparse(size(K_1)+size(K_2));
M = sparse(size(K_1)+size(K_2));

[sz1 sz2] = size(K_1);
[sz3 sz4] = size(K_2);
K(1:sz1,1:sz2) = K_1;
K(sz1+1:sz1+sz3,sz2+1:sz2+sz4) = K_2;

M(1:sz1,1:sz2) = M_1;
M(sz1+1:sz1+sz3,sz2+1:sz2+sz4) = M_2;

    [autovector,ome] = eigs(K,M,10,'sm');
    omega = diag(sqrt(ome));
