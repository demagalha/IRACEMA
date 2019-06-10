clear all
close all
clc

%% ISOGEOMETRIC ANALYSIS VIBRATION OF A ROD

%% Load Data / Create Model
P = {[0 0 0 1],[1 0 0 1]}; % Control Points for 1D Rod
Model = Geometry('curve',1,[0 0 1 1],P); % 1D Rod Constructor
% Refinement Phase: First Degree Elevate, then h refine
Model.DegreeElevate(1,1); % Quadratic Elements
Model.KnotRefine([0.1:0.1:0.9],1); % Add 9 Knot Spans
P = Model.get_point_cell;
U = Model.U;
% Connectivity Arrays
[INN, IEN, nel, nen] = Model.get_connectivity;
ID = reshape(1:max(max(IEN)),1,max(max(IEN)));
LM = zeros(nen,nel);
for i = 1:nel
    LM(:,i)=reshape(ID(:,IEN(:,i)),nen,1);
end
%% Assembly
[x,wx] = getGP(Model.pu); % Quadrature Rule (Gauss Legendre)
N_QUAD_X = length(x);
pu = Model.pu;
rho = 1; % Density
E = 1; % Elastic 
K = sparse(numel(INN),numel(INN));
M = K;
F = sparse(numel(INN),1);
N_DOF = numel(INN);
N_ELE_DOF = nen;

for e=1:nel % Loop Through Elements
    ni = INN(IEN(1,e),1);
    % Check if element has zero measure
    if U(ni+1) == U(ni)
         continue
    end
    K_e = zeros(nen,nen);
    M_e = K_e;
    F_e = zeros(nen,1);
    for i=1:N_QUAD_X % Loop through quadrature points
        [R, dR, J] = Shape1D(x(i),e,pu,P,U,INN,IEN);
        Jmod = J*wx(i);
        K_e = K_e +Jmod*dR'*E*dR;
        M_e = Jmod*rho*R'*R;
    end
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

%% Post-Processing
[autovector,omega] = eig(full(K),full(M));
omega = diag(omega);
[omega, idx] = sort(omega,'ascend');
autovector = autovector(:,idx);
omega_t = zeros(size(omega)); % Theoretical value
for i=1:length(omega)
    omega_t(i) = i*pi;
end
plot(omega./omega_t)
