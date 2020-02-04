clear all
% close all
clc

%% Pre-processing
% Constructing the Model
U = [0 0 1 1];
V = [0 0 1 1];

P{1,1} = [-5e-3 0 0 1]; P{1,2} = [-5e-3 0.1e-3 0 1];
P{2,1} = [5e-3 0 0 1]; P{2,2} = [5e-3 0.1e-3 0 1];
Model = Geometry('surf',1,U,1,V,P);

% k-refinement of the Model
% Model.DegreeElevate(5,1);
Model.DegreeElevate(1,2);
Model.KnotRefine(1/20:1/20:1-1/20,1);
% Model.KnotRefine(1/4:1/4:1-1/4,2);

%Material Properties
E = 4.5*(10^9); % Young Modulus
vu = 0.4; % Poisson Ratio
C = (E/(1-vu^2))*[1 vu 0; vu 1 0; 0 0 0.5*(1-vu)]; % Stiffness Tensor

%Body Force in y direction
f = @(x) 10000*[0; x.^2 - 25];

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
    F_e = zeros(2,nen);
    for i=1:N_QUAD_U
        for j=1:N_QUAD_V
            [R, dR, J] = Shape2D(Model,u(i),v(j),e,P,IEN,INN);
            Jmod = abs(J*wu(i)*wv(j));
            K_e = K_e + BuildKLocal2D(dR,Jmod,C);
            F_e = F_e+BuildFLocal2D(R,Jmod,Model,u(i),ni,v(j),nj,f);
        end
    end
    F_e = reshape(F_e,[nen*2,1]);
    % Assembly
    idx = LM(:,e);
    for i=1:N_ELE_DOF
        ii = idx(i);
        F(ii) = F(ii) + F_e(i);
        for j=1:N_ELE_DOF
            jj = idx(j);
            K(ii,jj) = K(ii,jj) +K_e(i,j);
        end
    end
end


%% Boundary Conditions

constNod = [];
for i=1:numel(P)
    if (P{i}(1)==5e-3 || P{i}(1)==(-5e-3))
        constNod = [constNod i];
    end
end
bc = reshape(ID(:,constNod),numel(ID(:,constNod)),1);

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
scaling_factor = 1e4;

for i=1:size(ID,2)
    u{i} = scaling_factor*[d(2*(i-1)+1), d(2*(i-1)+2), 0, 0];
    comb{i} = B{i} +u{i};
end

DeformedModel = Geometry('surf',Model.pu,Model.U,Model.pv,Model.V,comb);