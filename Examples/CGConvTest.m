clear all
close all
clc

% Membrane Geometry
a = 1;
b = 1.5;
P1 = [0 0 0 1];
P2 = [0 a 0 1];
P3 = [b 0 0 1];
P4 = [b a 0 1];
Refinements = [1,1; 2,2; 3,3; 4,4; 5,5; 6,6];
strings = {'GL11', 'GL22', 'GL33','GL44','GL55', 'GL66'};
for pp=1:6
m1 = Geometry('curve',1,[0 0 1 1],{P1,P2});
m2 = Geometry('curve',1,[0 0 1 1],{P3, P4});
Model = geo_ruled(m1,m2);
% Refinement
Model.DegreeElevate(Refinements(pp,1)+1,1);
Model.DegreeElevate(Refinements(pp,2),2);
interval = linspace(0,1,10);
interval = interval(2:end-1);
Model.KnotRefine(interval,1);
Model.KnotRefine(interval,2);
tic;

%% Assembly
  [INN, IEN, nel, nen] = Model.get_connectivity;
    ID = reshape(1:max(max(IEN)),1,max(max(IEN)));
    LM = zeros(nen,nel);
    for i=1:nel
        LM(:,i) = reshape(ID(:,IEN(:,i)),nen,1);
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
    N_DOF = numel(INN)/2;
    K = zeros(N_DOF);
    M = zeros(N_DOF);
    N_ELE_DOF = nen;
    check1 = min(INN);
    check2 = max(INN);    
    for e=1:nel
        ni = INN(IEN(1,e),1);
        nj = INN(IEN(1,e),2);
        % Check if element has measure zero
        if (U(ni+1) == U(ni) || V(nj+1) == V(nj))
            continue
        end
        K_e = zeros(N_ELE_DOF);
        M_e = zeros(N_ELE_DOF);
%         cond1 = (ni == check1(1) || nj == check1(2));
%         cond2 = (ni == check2(1) || nj == check2(2));
%         if cond1 || cond2
%             % If on the boundary, use Gauss-Legendre Quadrature
%             [u, wu] = getGP(pu);
%             [v, wv] = getGP(pv);
%         else
% % %             % If on the interior, use Cauchy-Galerkin Points for Quadrature
%             [u, wu] = getCG(pu);
%             [v, wv] = getCG(pv);
% 
%         end
% N_QUAD_U = length(u);
% N_QUAD_V = length(v);
        for i=1:N_QUAD_U
            for j=1:N_QUAD_V
                [R, dR, J] = Shape2D(Model,u(i),v(j),e,P,IEN,INN);
                Jmod = J*wu(i)*wv(j);
                K_e = K_e + Jmod*(dR(:,1)*dR(:,1)' +dR(:,2)*dR(:,2)');
                M_e = M_e+ R*R'*Jmod;
            end
        end

        % Assembly
        idx = LM(:,e);
        for i=1:N_ELE_DOF
            ii = idx(i);
            for j=1:N_ELE_DOF
                jj = idx(j);
                K(ii,jj) = K(ii,jj) +K_e(i,j);
                M(ii,jj) = M(ii,jj) +M_e(i,j);
            end
        end
    end

%% Boundary Conditions
constNod1 = [];
constNod2 = [];

P = Model.get_point_cell;
for i=1:numel(P)
    if P{i}(1) == 0 || abs(P{i}(1) - b) <= sqrt(eps)
        constNod1 = [constNod1 i];
    end
    if P{i}(2) == 0 || abs(P{i}(2) - a) <= sqrt(eps)
        constNod2 = [constNod2 i];
    end
end

bc = reshape(ID(1,[constNod1 constNod2]),numel(ID(1,[constNod1 constNod2])),1);
bc = sort(unique(bc),'descend');
for i=1:numel(bc)
    K(bc(i),:) = [];
    K(:,bc(i)) = [];
    M(bc(i),:) = [];
    M(:,bc(i)) = [];
end

% Results
K = sparse(K);
M = sparse(M);
[V,O] = eigs(K,M,length(K),'sm');
omega = sqrt(diag(O));
t = toc
save(strings{pp},'omega','t');
end

% Put back removed DOFs
% bc = sort(

% u = cell(size(P));
% comb = u;
% d = V(:,1);
% 
% for i=1:size(ID,2)
%     u{i} = [0 0 d(ID(:,i)) 0];
%     comb{i} = P{i} +u{i};
% end
ww = @(n,m) pi*sqrt((n/a)^2 +(m/b)^2);