clear all
close all
clc
%% Pre-processing
% Construction of the Model
U = [0 0 1 1];
V = [0 0 1 1];
P{1,1} = [-5e-3 0 0 1]; P{1,2} = [-5e-3 0.1e-3 0 1];
P{2,1} = [0 0 0 1]; P{2,2} = [0 0.1e-3 0 1];
Patch1 = Geometry('surf',1,U,1,V,P);
P{1,1} = [0 0 0 1]; P{1,2} = [0 0.1e-3 0 1];
P{2,1} = [5e-3 0 0 1]; P{2,2} = [5e-3 0.1e-3 0 1];
Patch2 = Geometry('surf',1,U,1,V,P);

% k-refinement of the Models
Patch1.DegreeElevate(1,2);
Patch1.KnotRefine(1/10:1/10:1-1/10,1);
Patch2.DegreeElevate(1,2);
Patch2.KnotRefine(1/10:1/10:1-1/10,1);

% Material Properties
E = 4.5*(10^9); % Young Modulus
vu = 0.4; % Poisson Ratio
C = (E/(1-vu^2))*[1 vu 0; vu 1 0; 0 0 0.5*(1-vu)]; % Stiffness Tensor

%Body Force in y direction
f = @(x) 10000*[0; x.^2 - 25];
P = Patch1.get_point_cell;
B = Patch2.get_point_cell;
[INN1, IEN1, nel(1), nen(1)] = Patch1.get_connectivity;
[INN2, IEN2, nel(2), nen(2)] = Patch2.get_connectivity;
ID1 = reshape(1:max(max(IEN1))*2,2,max(max(IEN1)));
ID2 = reshape(1:max(max(IEN2))*2,2,max(max(IEN2)));

%% Boundary Conditions and Multipatch MasterSlave indexing
constNod = [];
const2Nod = [];
overlap1Nod = [];
overlap2Nod = [];

for i=1:numel(P)
    if (P{i}(1)==-5e-3)
        constNod = [constNod i];
    end
    if (P{i}(1) ==0)
        overlap1Nod = [overlap1Nod i];
    end
end
bc1 = reshape(ID1(:,constNod),numel(ID1(:,constNod)),1);
for i=1:numel(B)
    if (B{i}(1) ==5e-3)
        const2Nod = [const2Nod i];
    end
    if (B{i}(1) == 0)
        overlap2Nod = [overlap2Nod i];
    end
end
overlap(:,1) = reshape(ID1(:,overlap1Nod),numel(ID1(:,overlap1Nod)),1);
overlap(:,2) = reshape(ID2(:,overlap2Nod),numel(ID2(:,overlap2Nod)),1);
% Update Slave's ID vector with correct indexes
ID2(:,overlap2Nod) = ID1(:,overlap1Nod);
maxadd = 1;
for j=1:length(ID2)
    % Check if element is a slave element
    if ismember(j,overlap2Nod)
        continue
    elseif maxadd == 1
        % If it is not a slave, start the indexes from the last ID's
        % highest DOF
        ID2(:,j) = ID1(:,end) +2;
        maxadd = 0;
        continue
    end
    % Check if previous element is a slave element
    if ismember(j-1,overlap2Nod)
        ID2(:,j) = ID2(:,j-2)+2;
    else
        ID2(:,j) = ID2(:,j-1)+2;
    end
end
bc2 = reshape(ID2(:,const2Nod),numel(ID2(:,const2Nod)),1);
%% Assembly
[sz1, ~] = size(overlap);
K = zeros(numel(INN1)+numel(INN2)-sz1);
F = zeros(length(K),1);
Model = {Patch1, Patch2};
ID = {ID1, ID2};
IEN = {IEN1, IEN2};
LM = cell(size(Model));
    for i=1:numel(Model)
        LM{i} = zeros(2*nen(i),nel(i));
        for j=1:nel(i)
            LM{i}(:,j) = reshape(ID{i}(:,IEN{i}(:,j)),2*nen(i),1);
        end
    end

for z=1:length(Model)
    [INN, ~, nel, nen] = Model{z}.get_connectivity;
    % Model parameters
    pu = Model{z}.pu;
    pv = Model{z}.pv;
    U = Model{z}.U;
    pu = Model{z}.pu;
    nu = length(U)-pu-2;
    V = Model{z}.V;
    pv = Model{z}.pv;
    nv = length(V)-pv-2;
    P = Model{z}.get_point_cell;

    % Gauss-Legendre Quadrature Points
    [u, wu] = getGP(pu);
    [v, wv] = getGP(pv);

    N_QUAD_U = length(u);
    N_QUAD_V = length(v);

    % Alocate Memory for Stiffness and Force Arrays
    N_DOF = numel(INN);
    N_ELE_DOF = 2*nen;
    for e=1:nel
        ni = INN(IEN{z}(1,e),1);
        nj = INN(IEN{z}(1,e),2);
        % Check if element has measure zero
        if (U(ni+1) == U(ni) || V(nj+1) == V(nj))
            continue
        end
        K_e = zeros(N_ELE_DOF);
        F_e = zeros(2,nen);
        for i=1:N_QUAD_U
            for j=1:N_QUAD_V
                [R, dR, J] = Shape2D(Model{z},u(i),v(j),e,P,IEN{z},INN);
                Jmod = abs(J*wu(i)*wv(j));
                K_e = K_e + BuildKLocal2D(dR,Jmod,C);
                F_e = F_e+BuildFLocal2D(R,Jmod,Model{z},u(i),ni,v(j),nj,f);
            end
        end
        F_e = reshape(F_e,[nen*2,1]);
        % Assembly
        idx = LM{z}(:,e);
        for i=1:N_ELE_DOF
            ii = idx(i);
            F(ii) = F(ii) + F_e(i);
            for j=1:N_ELE_DOF
                jj = idx(j);
                K(ii,jj) = K(ii,jj) +K_e(i,j);
            end
        end
    end
end




%% Applying Constraints
for i=1:numel(bc1)
    K(bc1(i),bc1(i)) = 1e30;
end
for i=1:numel(bc2)
    K(bc2(i),bc2(i)) = 1e30;
end

%% Solving
d = K\F;
figure(1)
subplot(1,2,1)
for z=1:2
    B = Model{z}.get_point_cell;
    u = cell(size(B));
    comb = u;
    scaling_factor = 10000;

    for i=1:size(ID{z},2)
        u{i} = scaling_factor*[full(d(ID{z}(:,i)))' 0 0];
        comb{i} = B{i} +u{i};
    end

    DeformedModel{z} = Geometry('surf',Model{z}.pu,Model{z}.U,Model{z}.pv,Model{z}.V,comb);
    DeformedModel{z}.plot_geo;
end
xlim([-5e-3 5e-3])
ylim([-5e-3 5e-3])
