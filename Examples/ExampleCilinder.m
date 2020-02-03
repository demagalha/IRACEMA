clear all
close all
clc

%% Pre-processing
R = 1; % MidSurface Radius
t = 0.01; % Thickness
h = 5; % Height
Model = BuildQuarterCilinder(R+t/2,R-t/2,h); % Example Cilinder from Book
% Model.DegreeElevate(3,1);
% Model.DegreeElevate(4,2);
% Model.DegreeElevate(1,3);
Model.KnotRefine([0.1:0.1:0.9],3);
Model.KnotRefine([1/10 1/3 0.5 2/3 9/10],1);
Model.KnotRefine([1/10 1/3 0.5 2/3 9/10],2);
YOUNG_MODULUS =  3*10^6;
POISSON = 0.2;
D = get_matprop_matrix(1,YOUNG_MODULUS,POISSON);
RHO = 1;
% f = @(x,y,z) [x/(R-t); y/(R-t); 0*z]; %[ cos(theta) ; sin(theta); 0] theta = atan(y/x)
%% Assembly
f = @(x,y,z) InternalPressure(x,y,z,R);
[K,F,IEN] = ElastoAssemble(Model,D,f);

%% Boundary Conditions
P = Model.get_point_cell;
ID = reshape(1:max(max(IEN))*3,3,max(max(IEN)));
constNod = [];
sym1Nod = [];
sym2Nod = [];
% Get Indexes where there is constraint
% Get Indexes where there is symmetry plane
for i=1:numel(P)
    if (P{i}(3) == 0 || P{i}(3) == 5)
        constNod = [constNod i];
    end
    if P{i}(1) == 0
        sym1Nod = [sym1Nod i];
    end
    if P{i}(2) == 0
        sym2Nod = [sym2Nod i];
    end
end

constraints = reshape(ID(:,constNod),numel(ID(:,constNod)),1);
symmetry1 = reshape(ID(1,sym1Nod),numel(ID(1,sym1Nod)),1);
symmetry2 = reshape(ID(2,sym2Nod),numel(ID(2,sym2Nod)),1);
bc = [constraints; symmetry1; symmetry2];
% Apply no displacement constraints
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
scaling_factor = 1e5;

for i=1:size(ID,2)
    u{i} = scaling_factor*[full(d(ID(:,i)))' 0];
    comb{i} = B{i} +u{i};
end
DeformedModel = Geometry('volume',Model.pu,Model.U,Model.pv,Model.V,Model.pw,Model.W,comb);
DeformedModel.plot_geo('coarse',0,1);
shading interp
