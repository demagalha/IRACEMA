clear all
close all
clc

%% Geometry of the Beam
L = 1;
P1 = [0 0 0 1];
P2 = [L 0 0 1];
Model = Geometry('curve',1,[0 0 1 1],{P1,P2});

% Beam Properties
% Steel Beam
b = 150e-3;
h = 5e-3;
rho = 7900;
A = b*h;
I = b*h*h*h/12;
E = 200*(10^9);
a = sqrt((E*I)/(rho*A));

%% k-refinement of the Model
p = 2; % Number of p-refinements
r = 1000; % Number of h-refinements
interval = linspace(0,1,r+2);
interval = interval(2:end-1);
Model.DegreeElevate(p,1);
Model.KnotRefine(interval,1);

%% Assembly
[INN, IEN, nel nen] = Model.get_connectivity;
ID = reshape(1:max(max(IEN)),1,max(max(IEN)));
LM = zeros(nen,nel);
for i=1:nel
    LM(:,i) = reshape(ID(:,IEN(:,i)),nen,1);
end

% Model Parameters
pu = Model.pu;
U = Model.U;
P = Model.get_point_cell;

NDOF = numel(INN);
NELEDOF = nen;
K = zeros(NDOF);
M = K;
check1 = min(INN);
check2 = max(INN);
for e=1:nel
    ni = INN(IEN(1,e),1);
    if (U(ni+1) == U(ni))
        continue
    end
    K_e = zeros(nen);
    M_e = K_e;
    if ni == check1 || ni == check2
        [u, wu] = getGP(pu);
    else
        [u,wu] = getCG(pu);
    end
    NQUADU = length(u);
    for i=1:NQUADU
        [R, d2R, J] = BeamShape(Model,u(i),e,P, IEN, INN);
        Jmod = abs(J*wu(i));
        K_e = K_e + d2R'*E*I*d2R*Jmod;
        M_e = R*rho*A*R'*Jmod;
    end
    % Assemblage
    idx = LM(:,e);
    for i=1:NELEDOF
        ii = idx(i);
        for j=1:NELEDOF
            jj = idx(j);
            K(ii,jj) = K(ii,jj) + K_e(i,j);
            M(ii,jj) = M(ii,jj) + M_e(i,j);
        end
    end
end

%% Boundary Conditions

constNod = [];
for i=1:numel(P)
    if P{i}(1) == 0 || P{i}(1) == L
        constNod = [constNod i];
    end
end
bc = reshape(ID(1,constNod),numel(ID(1,constNod)),1);
bc = sort(bc,'descend');
for i=1:numel(bc)
    if bc(i) < 3
        K(:,bc(i):bc(i)+2) = [];
        K(bc(i):bc(i)+2,:) = [];
        M(:,bc(i):bc(i)+2) = [];
        M(bc(i):bc(i)+2,:) = [];
    else
        K(:,bc(i)-2:bc(i)) = [];
        K(bc(i)-2:bc(i),:) = [];
        M(:,bc(i)-2:bc(i)) = [];
        M(bc(i)-2:bc(i),:) = [];
    end
end

K = sparse(K);
M = sparse(M);

%% Results

[autovector,omega] = eigs(K,M,length(K),'sm');
bc = sort(bc,'ascend');
[sz1,sz2] = size(autovector);
boundaries = zeros(1,sz2);
bc_autos = zeros(sz1+numel(bc),sz2);
for i=1:numel(bc)
    autovector = [autovector(1:bc,:); boundaries; autovector(bc+1:end,:)];
end
omega = sqrt(diag(omega));
[omega, idx] = sort(omega,'ascend');
autovector = autovector(:,idx);
