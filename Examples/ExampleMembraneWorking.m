close all
clearvars
clc

%% Geometry Input
a = 1;
b = 1.5;
% Membrane
    P1 = [0 0 0 1];
    P2 = [0 a 0 1];
    P3 = [b 0 0 1];
    P4 = [b a 0 1];
    c1 = Geometry('curve',1,[0 0 1 1],{P1 P2});
    c2 = Geometry('curve',1,[0 0 1 1],{P3 P4});
    Model = geo_ruled(c1,c2);
    
clearvars c1 c2 P1 P2 P3 P4

%% Refinement of the Models
% p refinement
    p = 1; % Number of p-refinements
    Model.DegreeElevate(p,1);
    Model.DegreeElevate(p,2);
    
% h refinement
    r = 13;% Number of h-refinements
    interval = linspace(0,1,r+2);
    interval = interval(2:end-1); % This is setup so it has exactly r new entries
    Model.KnotRefine(interval,1);
    Model.KnotRefine(interval,2);
    
    
%% Assembly
[K, M, ID] = MembraneAssemble(Model);
%% Constraints
% Domain definition
%          
%u=1_________________________ 
% |                         |
% |                         |
% |                         |
% | v = 0               v=1 |Gamma
% |                         |
% |                         |
% |_________________________|
% u=0       Gamma

[INN, IEN, ~, ~] = Model.get_connectivity;
ut_zero = find(INN(:,1) == 1);
ut_one = find(INN(:,1) == max(INN(:,1)));
vt_zero = find(INN(:,2) == 1);
vt_one = find(INN(:,2) == max(INN(:,2)));
Gamma = union(ut_zero,ut_one);
Gamma = union(Gamma,vt_zero);
Gamma = union(Gamma,vt_one);
Gamma = reshape(ID(:,Gamma),numel(ID(:,Gamma)),1);

%% Apply Constraints
% First, we have to take out the constrained DOFs from K and M matrices
K(Gamma,:) = [];
K(:,Gamma) = [];
M(Gamma,:) = [];
M(:,Gamma) = [];
% Now we solve the Eigenvalue Problem
K = sparse(K);
M = sparse(M);

graphs = [3,2]; % How you want to subplot your modes 
n = prod(graphs); % Take the first eigenvalues
[autovectors, eigenvalues] = eigs(K,M,n,'sm');
omega = sqrt(diag(eigenvalues));
% And we add back the constrained dofs to the autovectors in ascending
% order
Gamma = sort(Gamma,'ascend');
autovectors = BoundariesPostProcess(autovectors,Gamma);

% Visualization of the Modes / Construction of a deformed model
Modes = VisualizeModes(Model,autovectors,ID);

figure(1)
for i=1:n
    subplot(graphs(1),graphs(2),i)
    Modes{i}.plot_geo
    title(strcat('Modo #',num2str(i),' , \omega_n = ', num2str(omega(i))));
end