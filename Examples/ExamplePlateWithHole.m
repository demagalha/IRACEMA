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
Model.KnotRefine(1/11:1/11:1-1/11,1);
Model.KnotRefine(1/11:1/11:1-1/11,2);

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
T = -10;

sym1 = reshape(ID(2,sym1Nod),numel(ID(2,sym1Nod)),1);
sym2 = reshape(ID(1,sym1Nod),numel(ID(1,sym1Nod)),1);
bc = [sym1; sym2];
for i=1:numel(bc)
    K(bc(i),bc(i)) = 1e30;
end


f1 = reshape(ID(1,f1Nod),numel(ID(1,f1Nod)),1);
f1values = T*ones(size(f1));
f2 = reshape(ID(1,f2Nod),numel(ID(1,f2Nod)),1);
f2values = T*ones(size(f2));
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
uu = u;
comb = u;
scaling_factor = 33;

for i=1:size(ID,2)
    u{i} = scaling_factor*[full(d(ID(:,i)))' 0 0];
    uu{i} = [full(d(ID(:,i)))' 0 1];
    comb{i} = B{i} +u{i};
end

DeformedModel = Geometry('surf',Model.pu,Model.U,Model.pv,Model.V,comb);
% DeformedModel.plot_geo('fine',0,1);

u = 1;
v = 1;
U = Model.U;
pu = Model.pu;
V = Model.V;
pv = Model.pv;

su = FindSpanLinear(length(U)-pu-2,pu,u,U);
sv = FindSpanLinear(length(V)-pv-2,pv,v,V);

NU = DersBasisFun(su,u,pu,1,U);
N = NU(1,:);
dN = NU(2,:);
[~, si] = size(NU);

MV = DersBasisFun(sv,v,pv,1,V);
M = MV(1,:);
dM = MV(2,:);
[~,sj] = size(MV);
clear NU MV

P = Model.get_point_cell;
weight = Model.weight;
P = P(su-pu+1:su+1,sv-pv+1:sv+1);
weight = weight(su-pu+1:su+1,sv-pv+1:sv+1);
Q = 0;
dQdu = 0;
dQdv = 0;
dBdu = zeros(size(B));
dBdv = dBdu;
B = zeros(1,si*sj);
for idx=1:si*sj
    [i,j] = ind2sub([si,sj],idx);
    B(idx) = N(i)*M(j);
    dBdu(idx) = dN(i)*M(j);
    dBdv(idx) = N(i)*dM(j);
    Q = Q+B(idx)*weight(idx);
    dQdu = dQdu+ dBdu(idx)*weight(idx);
    dQdv = dQdv+ dBdv(idx)*weight(idx);
end
R = zeros(1,si*sj);
dRdu = zeros(size(R));
dRdv = zeros(size(R));
for idx=1:si*sj
    [i,j] = ind2sub([si,sj],idx);
    R(idx) = B(idx)/Q;
    dRdu(idx) = dN(i)*M(j) - (R(idx)/Q)*dQdu;
    dRdv(idx) = N(i)*dM(j) - (R(idx)/Q)*dQdv;
end

dR = zeros(2,length(R));
dRdx = dR;
dR(1,:) = dRdu;
dR(2,:) = dRdv;
dxdu = zeros(2);
for idx=1:si*sj
    [i,j] = ind2sub([si,sj],idx);
    for xx=1:2
        for yy=1:2
            dxdu(xx,yy) = dxdu(xx,yy) +P{i,j}(xx)*dR(yy,idx);
        end
    end
end
dudx = inv(dxdu);

for idx=1:si*sj
    for xx=1:2
        for yy=1:2
            dRdx(xx,idx) = dRdx(xx,idx) +dR(yy,idx)*dudx(yy,xx);
        end
    end
end
Jmod = det(dxdu);
[sz1, sz2] = size(uu);
uu = uu(su-pu+1:su+1,sv-pv+1:sv+1);
xx = zeros(numel(uu),1);
yy = xx;
for i=1:numel(uu)
    xx(i) = uu{i}(1)*uu{i}(4);
    yy(i) = uu{i}(2)*uu{i}(4);
end

sx = dRdx(1,:)*xx(:);
sy = dRdx(2,:)*yy(:);
sxy = dRdx(1,:)*yy(:) + dRdx(2,:)*xx(:);

sigma = C*[sx; sy; sxy]


    