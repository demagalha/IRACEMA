clear all
close all
clc

%% Pre-Processing
L = 4;
R = 1;
Model = BuildHoledPlate(L,R);

% Refinement
Model.DegreeElevate(3,1);
Model.DegreeElevate(4,2);
Model.KnotRefine(0.2:0.2:0.8,1);
Model.KnotRefine(0.2:0.2:0.8,2);

% Material Properties
E = 1e5;
vu = 0.3;
C = (E/(1-vu^2))*[1 vu 0; vu 1 0; 0 0 0.5*(1-vu)]; % Stiffness Tensor

%% Assembly
[INN, IEN, nel, nen] = Model.get_connectivity;
ID = reshape(1:max(max(IEN))*2,2,max(max(IEN)));
LM = zeros(2*nen,nel);
for i=1:nel
    LM(:,i) = reshape(ID(:,IEN(:,i)),2*nen,1);
end
% Model parameters
pu = Model.pu;
pv = Model.pv;
U = Model.U;
pu = Model.pu;
nu = length(U)-pu-2;
V = Model.V;
pv = Model.pv;
nv = length(V)-pv-2;
P = Model.get_point_cell;

% Gauss-Legendre Quadrature Points
[u, wu] = getGP(pu);
[v, wv] = getGP(pv);

N_QUAD_U = length(u);
N_QUAD_V = length(v);

% Alocate Memory for Stiffness and Force Arrays
N_DOF = numel(INN);
K = zeros(N_DOF);
F = zeros(N_DOF,1);
N_ELE_DOF = 2*nen;
for e=1:nel
    ni = INN(IEN(1,e),1);
    nj = INN(IEN(1,e),2);
    % Check if element has measure zero
    if (U(ni+1) == U(ni) || V(nj+1) == V(nj))
        continue
    end
    K_e = zeros(N_ELE_DOF);
%     F_e = zeros(2,nen);
    for i=1:N_QUAD_U
        for j=1:N_QUAD_V
            [R, dR, J] = Shape2D(Model,u(i),v(j),e,P,IEN,INN);
            Jmod = abs(J*wu(i)*wv(j));
            K_e = K_e + BuildKLocal2D(dR,Jmod,C);
%             F_e = F_e+BuildFLocal2D(R,Jmod,Model,u(i),ni,v(j),nj,f);
        end
    end
%     F_e = reshape(F_e,[nen*2,1]);
    % Assembly
    idx = LM(:,e);
    for i=1:N_ELE_DOF
        ii = idx(i);
%         F(ii) = F(ii) + F_e(i);
        for j=1:N_ELE_DOF
            jj = idx(j);
            K(ii,jj) = K(ii,jj) +K_e(i,j);
        end
    end
end

%% Boundary Conditions
P = Model.get_point_cell;

sym1Nod = [];
sym2Nod = [];
f1Nod = [];
f2Nod = [];

for i=1:numel(P)
    if P{i}(2) == 0
        sym1Nod = [sym1Nod i];
    end
    if P{i}(1) == L
        sym2Nod = [sym2Nod i];
    end
    if P{i}(2) == L
        f1Nod = [f1Nod i];
    end
    if P{i}(1) == 0
        f2Nod = [f2Nod i];
    end
end
T = 10;

sym1 = reshape(ID(2,sym1Nod),numel(ID(2,sym1Nod)),1);
sym2 = reshape(ID(1,sym1Nod),numel(ID(1,sym1Nod)),1);
bc = [sym1; sym2];
for i=1:numel(bc)
    K(bc(i),bc(i)) = 1e30;
end


f1 = reshape(ID(2,f1Nod),numel(ID(2,f1Nod)),1);
f1values = T*ones(size(f1));
f2 = reshape(ID(1,f2Nod),numel(ID(1,f2Nod)),1);
f2values = -T*ones(size(f2));
fbc = [f1, f1values; f2, f2values];
for i=1:length(fbc)
    F(fbc(i,1)) = fbc(i,2);
end
K = sparse(K);
F = sparse(F);

%% Solution
d = K\F;

B = Model.get_point_cell;
u = cell(size(B));
comb = u;
scaling_factor = 33;

for i=1:size(ID,2)
    u{i} = scaling_factor*[d(2*(i-1)+1), d(2*(i-1)+2), 0, 0];
    comb{i} = B{i} +u{i};
end

DeformedModel = Geometry('surf',Model.pu,Model.U,Model.pv,Model.V,comb);
DeformedModel.plot_geo('fine',0,1);