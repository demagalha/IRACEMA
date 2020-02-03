clear all
close all
clc

%% Pre-processing
R = 300;
t = 3;
L = 600;
Model = BuildCilinder(R,R-t,L);

Model.DegreeElevate(3,1);
Model.DegreeElevate(4,2);
Model.DegreeElevate(2,3);
Model.KnotRefine(0.5,1);
Model.KnotRefine(0.5,2);
Model.KnotRefine([0.25:0.25:1],3);
% Model.plot_geo('coarse',1,0);
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
fNod = [];
for i=1:numel(P)
    if (P{i}(3) == 0 || P{i}(3) == 600)
        constNod = [constNod i];
    end
    if (P{i}(2) == 0 && P{i}(3) == R)
        fNod = [fNod i];
    end
end
    
bc = reshape(ID(:,constNod),numel(ID(:,constNod)),1);
fc = reshape(ID(1,fNod),numel(ID(1,fNod)),1);
for i=1:numel(bc)
    K(bc(i),bc(i)) = 1e30;
end
for i=1:numel(fc)
    if P{fNod(i)}(1) < 0
        F(fc(i)) = 1;
    else
        F(fc(i)) = -1;
    end
end
K = sparse(K);
F = sparse(F);
%% Solution
d = K\F;
B = Model.get_point_cell;
u = cell(size(B));
comb = u;
scaling_factor = 0.5e5;

for i=1:size(ID,2)
    u{i} = scaling_factor*[full(d(ID(:,i)))' 0];
    comb{i} = B{i} +u{i};
end
DeformedModel = Geometry('volume',Model.pu,Model.U,Model.pv,Model.V,Model.pw,Model.W,comb);
DeformedModel.plot_geo('coarse',0,1);
xlim([-500 500])
ylim([-500 500])
shading interp

