clear all
close all
clc

%% Pre-processing
R = 300;
t = 3;
L = 300;
Model = BuildQuarterCilinder(R,R-t,L);

Model.DegreeElevate(3,1);
Model.DegreeElevate(4,2);
Model.DegreeElevate(1,3);
Model.KnotRefine(1/20:1/20:1-1/20,1);
Model.KnotRefine(1/20:1/20:1-1/20,2);
Model.KnotRefine(0.5,3);
% Model.plot_geo('coarse',1,0);

% 'Local' Refinement for Force Boundary Condition
% Model.KnotRefine(0.01,1);
% Model.KnotRefine(0.01,3);
% Material Properties
E = 3e6;
vu = 0.3;
D = get_matprop_matrix(1,E,vu);

%% Assembly
f = @(x,y,z) [0*x; 0*y; 0*z];
[K,F,IEN] = ElastoAssemble(Model,D,f);

%% Boundary Conditions
P = Model.get_point_cell;
ID = reshape(1:max(max(IEN))*3,3,max(max(IEN)));
constNod = [];
sym1Nod = [];
sym2Nod = [];
fNod = [];
for i=1:numel(P)
    if P{i}(3) == 0
        sym1Nod = [sym1Nod i];
    end
    if P{i}(1) == 0
        sym2Nod = [sym2Nod i];
    end
    if P{i}(3) == L
        constNod = [constNod i];
    end
    if (P{i}(3) == 0) && (P{i}(1) == R) && (P{i}(4) == 1)
        fNod = [fNod i];
    end
end
    
bc = reshape(ID(1,constNod),numel(ID(1,constNod)),1);
sym1 = reshape(ID(3,sym1Nod),numel(ID(3,sym1Nod)),1);
sym2 = reshape(ID(1,sym2Nod),numel(ID(1,sym2Nod)),1);
fc = reshape(ID(1,fNod),numel(ID(1,fNod)),1);

bc = [bc; sym1; sym2];
for i=1:numel(bc)
    K(bc(i),bc(i)) = 1e30;
end
for i=1:numel(fc)
    F(fc(i)) = -1/4;
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
DeformedModel.plot_geo('coarse',0,1);
xlim([0 300])
ylim([0 300])
% shading interp

