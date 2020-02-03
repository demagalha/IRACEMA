clear all
close all
clc

%% Pre-processing
Model = BuildQuarterHemisphere(10,0.04);

% Refinement
Model.DegreeElevate(3,1);
Model.DegreeElevate(3,2);
Model.DegreeElevate(4,3);
Model.KnotRefine(0.5,1);
Model.KnotRefine(0.5,2);
Model.KnotRefine(0.5,3);

% Material Properties
E = 6.825*(10^7);
vu = 0.3;
D = get_matprop_matrix(1,E,vu);

%% Assembly
f = @(x,y,z) [0; 0; 0];
[K,F,IEN] = ElastoAssemble(Model,D,f);

%% Boundary Conditions
P = Model.get_point_cell;
ID = reshape(1:max(max(IEN))*3,3,max(max(IEN)));
constNod = [];
sym1Nod = [];
sym2Nod = [];
f1Nod = [];
f2Nod = [];
for i=1:numel(P)
    C(i,1:3) = P{i}(1:3);
end
R = max(C(:,1));
xmax = R;
zmax = R;
ymax = R;
for i=1:numel(P)
    if (P{i}(1)==0) && (P{i}(2) ==0) && (P{i}(3) == 10.02)
        constNod = [constNod i];
    elseif P{i}(2) == 0
        sym1Nod = [sym1Nod i];
    elseif P{i}(1) == 0
        sym2Nod = [sym2Nod i];
    end
    if (P{i}(1) == 0) && (P{i}(2) == 10.02)
        f1Nod = [f1Nod i];
    end
    if (P{i}(1) == 10.02) && (P{i}(2) == 0)
        f2Nod = [f2Nod i];
    end
end
constraints = reshape(ID(:,constNod),numel(ID(:,constNod)),1);
symmetry1 = reshape(ID(2,sym1Nod),numel(ID(2,sym1Nod)),1);
symmetry2 = reshape(ID(1,sym2Nod),numel(ID(1,sym2Nod)),1);
bc = [constraints; symmetry1; symmetry2];
for i=1:numel(bc)
    K(bc(i),bc(i)) = 1e30;
end
K = sparse(K);
f1 = reshape(ID(2,f1Nod),numel(ID(2,f1Nod)),1);
f2 = reshape(ID(1,f2Nod),numel(ID(1,f2Nod)),1);

for i=1:numel(f1)
    F(f1(i)) = -1;
end
for i=1:numel(f2)
    F(f2(i)) = 1;
end
F = sparse(F);

%% Solution
d = K\F;
B = Model.get_point_cell;
u = cell(size(B));
comb = u;
scaling_factor = 33;

for i=1:size(ID,2)
    u{i} = scaling_factor*[full(d(ID(:,i)))' 0];
    comb{i} = B{i} +u{i};
end
DeformedModel = Geometry('volume',Model.pu,Model.U,Model.pv,Model.V,Model.pw,Model.W,comb);
DeformedModel.plot_geo('coarse',0,1);
