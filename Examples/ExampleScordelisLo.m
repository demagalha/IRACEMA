clear all
close all
clc
%% Pre-processing

% for inter=1:5
load('ScordelisSurfTCima.mat');
Model = geo_extrusion(Scordelis,[0 0 -0.25]);
% if inter > 1
% Model.DegreeElevate(inter-1,1);
% Model.DegreeElevate(inter-1,2);
% end
% Model.DegreeElevate(inter-1,3);
r = 9;
interval = 0+1/r:1/r:1-1/r;
Model.DegreeElevate(3,1);
Model.DegreeElevate(3,2);
Model.DegreeElevate(1,3);
Model.KnotRefine(interval,1);
Model.KnotRefine(interval,2);
Model.KnotRefine(0.5,3);

% Material Properties
E = 4.32*(10^8);
vu = 0;
D = get_matprop_matrix(1,E,vu);

%% Assembly
f = @(x,y,z) [0; 0; -90/0.25];
[K,F,IEN] = ElastoAssemble(Model,D,f);

%% Boundary Conditions
P = Model.get_point_cell;
ID = reshape(1:max(max(IEN))*3,3,max(max(IEN)));
constNod = [];
sym1Nod = [];
sym2Nod = [];
for i=1:numel(P)
    if P{i}(2) == 25
        constNod = [constNod i];
    end
    if P{i}(2) == 0
        sym1Nod = [sym1Nod i];
    end
    if P{i}(1) == 0
        sym2Nod = [sym2Nod i];
    end
end
constraints = reshape(ID([1,3],constNod),numel(ID([1,3],constNod)),1);
symmetry1 = reshape(ID(1,sym1Nod),numel(ID(1,sym1Nod)),1);
symmetry2 = reshape(ID(2,sym2Nod),numel(ID(2,sym2Nod)),1);
bc = [constraints; symmetry1; symmetry2];
for i=1:numel(bc)
    K(bc(i),bc(i)) = 1e30;
end
K = sparse(K);
F = sparse(F);

%% Solution
d = K\F;
B = Model.get_point_cell;
u = cell(size(B));
comb = u;
scaling_factor = 1;

for i=1:size(ID,2)
    u{i} = scaling_factor*[full(d(ID(:,i)))' 0];
    comb{i} = B{i} +u{i};
end
DeformedModel = Geometry('volume',Model.pu,Model.U,Model.pv,Model.V,Model.pw,Model.W,comb);
% DeformedModel.plot_geo('coarse',0,1);
% shading interp
% Convergence(inter) = 
-u{1}(3)
% end