close all
clearvars
clc

%% Geometry Input
a = 1;
b = 1.5;
overlap = 0.75;
% Omega 1:
    P1 = [0 0 0 1];
    P2 = [0 a 0 1];
    P3 = [b*(1+overlap)/2 0 0 1];
    P4 = [b*(1+overlap)/2 a 0 1];
    c1 = Geometry('curve',1,[0 0 1 1],{P1 P2});
    c2 = Geometry('curve',1,[0 0 1 1],{P3 P4});
    Omega1 = geo_ruled(c1,c2);

% Omega 2:
    P1 = [b*(1-overlap)/2 0 0 1];
    P2 = [b*(1-overlap)/2 a 0 1];
    P3 = [b 0 0 1];
    P4 = [b a 0 1];
    c1 = Geometry('curve',1,[0 0 1 1],{P1 P2});
    c2 = Geometry('curve',1,[0 0 1 1],{P3 P4});
    Omega2 = geo_ruled(c1,c2);

% Whole Domain
    P1 = [0 0 0 1];
    P2 = [0 a 0 1];
    P3 = [b 0 0 1];
    P4 = [b a 0 1];
    c1 = Geometry('curve',1,[0 0 1 1],{P1 P2});
    c2 = Geometry('curve',1,[0 0 1 1],{P3 P4});
    TheOmega = geo_ruled(c1,c2);
    
clearvars c1 c2 P1 P2 P3 P4

%% Refinement of the Models
% p refinement
    p = 1; % Number of p-refinements
    Omega1.DegreeElevate(p,1);
    Omega1.DegreeElevate(p,2);
    Omega2.DegreeElevate(p,1);
    Omega2.DegreeElevate(p,2);
    TheOmega.DegreeElevate(p,1);
    TheOmega.DegreeElevate(p,2);
    
% h refinement
    r = 13;% Number of h-refinements
    interval = linspace(0,1,r+2);
    interval = interval(2:end-1);
    Omega1.KnotRefine(interval,1);
    Omega1.KnotRefine(interval,2);
    Omega2.KnotRefine(interval,1);
    Omega2.KnotRefine(interval,2);
    TheOmega.KnotRefine(interval,1);
    TheOmega.KnotRefine(interval,2);
    
    
%% Inverse NURBS functions
    invNURBS1x = @(x) x/(b*(1+overlap)/2);
    conf1 = invNURBS1x(b*(1-overlap)/2);
    invNURBS1y = @(y) y;
    InverseFunOne = {invNURBS1x,invNURBS1y};
    
    invNURBS2x = @(x) invNURBS1x(x- b*(1-overlap)/2);
    conf2 = invNURBS2x(b*(1+overlap)/2);
    invNURBS2y = @(y) y;
    InverseFunTwo = {invNURBS2x, invNURBS2y};
    
%% Assembly
[K1,M1,ID1] = MembraneAssemble(Omega1);
K_one = K1; M_one = M1;  % Need these for later.

[K2,M2,ID2] = MembraneAssemble(Omega2);
K_two = K2; M_two = M2;

[KK, MM, IDD] = MembraneAssemble(TheOmega);
%% Constraints
% Domain definition
%% Omega 1

[INN1, IEN1, ~, ~] = Omega1.get_connectivity;
u1_zero = find(INN1(:,1) == 1);
u1_one = find(INN1(:,1) == max(INN1(:,1)));
v1_zero = find(INN1(:,2) == 1);
v1_one = find(INN1(:,2) == max(INN1(:,2)));

% Increase in v -> increase in x
% Increase in u -> increase in y 
%           Gamma1
%u=1_________________________ 
% G|                         |
% a|                         |
% m|                         |
% m| v = 0               v=1 |Gamma12
% a|                         |
% 1|                         |
%  |_________________________|
% u=0       Gamma1
Gamma1 = union(u1_one,u1_zero);
Gamma1 = union(Gamma1,v1_zero);
% Transform Gamma1 from basis function to dof ID
Gamma1 = reshape(ID1(:,Gamma1),numel(ID1(:,Gamma1)),1);

Gamma12 = v1_one;
% Gamma11 = union(Gamma1,Gamma12);
Gamma11 = Gamma1;
% Transform Gamma12 from basis function to dof ID
% Gamma12 = reshape(ID1(:,Gamma12),numel(ID1(:,Gamma12)),1);
%% Omega 2
[INN2, IEN2, ~, ~] = Omega2.get_connectivity;
u2_zero = find(INN2(:,1) == 1);
u2_one = find(INN2(:,1) == max(INN2(:,1)));
v2_zero = find(INN2(:,2) == 1);
v2_one = find(INN2(:,2) == max(INN2(:,2)));
% Increase in v -> increase in x
% Increase in u -> increase in y 
%           Gamma2
%u=1_________________________ 
% G|                         |
% a|                         |
% m|                         |
% m| v = 0               v=1 |Gamma2
% a|                         |
% 2|                         |
% 1|_________________________|
% u=0       Gamma2
Gamma2 = union(u2_one,u2_zero);
Gamma2 = union(Gamma2,v2_one);
% Transform Gamma2 from basis function to dof ID
Gamma2 = reshape(ID2(:,Gamma2),numel(ID2(:,Gamma2)),1);

Gamma21 = v2_zero;
% Gamma22 = union(Gamma2,Gamma21);
Gamma22 = Gamma2;
% Gamma22 = reshape(ID2(:,Gamma22),numel(ID2(:,Gamma22)),1);
% Transform Gamma21 from basis function to dof ID
% Gamma21 = reshape(ID2(:,Gamma21),numel(ID2(:,Gamma21)),1);

%% The Omega
[INNN, IENN, ~, ~] = TheOmega.get_connectivity;
ut_zero = find(INNN(:,1) == 1);
ut_one = find(INNN(:,1) == max(INNN(:,1)));
vt_zero = find(INNN(:,2) == 1);
vt_one = find(INNN(:,2) == max(INNN(:,2)));
TheGamma = union(ut_zero,ut_one);
TheGamma = union(TheGamma,vt_zero);
TheGamma = union(TheGamma,vt_one);
TheGamma = reshape(IDD(:,TheGamma),numel(IDD(:,TheGamma)),1);

%% Apply Constraints

% First, we have to take out the constrained DOFs from K and M matrices
K1(Gamma11,:) = [];
K1(:,Gamma11) = [];
M1(Gamma11,:) = [];
M1(:,Gamma11) = [];
K2(Gamma22,:) = [];
K2(:,Gamma22) = [];
M2(Gamma22,:) = [];
M2(:,Gamma22) = [];
KK(TheGamma,:) = [];
KK(:,TheGamma) = [];
MM(TheGamma,:) = [];
MM(:,TheGamma) = [];
% Now we solve the Eigenvalue Problem
K1 = sparse(K1);
M1 = sparse(M1);
K2 = sparse(K2);
M2 = sparse(M2);
KK = sparse(KK);
MM = sparse(MM);
[av1,O1] = eigs(K1,M1,1,'sm');
[av2,O2] = eigs(K2,M2,1,'sm');
[aa, OO] = eigs(KK,MM,1,'sm');

if abs(min(av1)) > abs(max(av1))
    av1 = -av1;
end
if abs(min(av2)) > abs(max(av2))
    av2 = -av2;
end
if abs(min(aa)) > abs(max(aa))
    aa = -aa;
end

% av1 = av1/max(av1);
% av2 = av2/max(av2);

% And we add back the constrained dofs to the autovectors
Gamma11 = sort(Gamma11,'ascend');
Gamma22 = sort(Gamma22,'ascend');
TheGamma = sort(TheGamma,'ascend');
av1 = BoundariesPostProcess(av1,Gamma11);
av2 = BoundariesPostProcess(av2,Gamma22);
aa = BoundariesPostProcess(aa,TheGamma);

% Visualization of the Modes / Construction of a deformed model
Omega1_Modes = VisualizeModes(Omega1,av1,ID1);
Omega2_Modes = VisualizeModes(Omega2,av2,ID2);
TheOmega_Modes = VisualizeModes(TheOmega,aa,IDD);

%% First Graphs
str1 = '\omega_1 = ';
str2 = '\omega_2 = ';
figure(1)
subplot(2,2,1)
Omega1_Modes{1}.plot_geo;
title(strcat('\Omega_1 Inicial, \omega_1 = ',num2str(sqrt(O1))))
xlim([0 1.5])
subplot(2,2,2)
Omega2_Modes{1}.plot_geo;
xlim([0 1.5])
title(strcat('\Omega_2 Inicial, \omega_2 = ',num2str(sqrt(O1))))
subplot(2,2,3)
Omega1_Modes{1}.plot_geo;
title('Overlap Inicial')
hold on
Omega2_Modes{1}.plot_geo;
xlim([0 1.5])
alpha(0.9)
subplot(2,2,4)
TheOmega_Modes{1}.plot_geo;
title(strcat('\Omega Objetivo, \omega = ',num2str(sqrt(OO))))
xlim([0 1.5])
%% Schwarz Alternating Algorithm
iter = 0;
graphs = [3,2];
max_iter = 6;
str3 = ', Iteration #';
while iter < max_iter
%% First Step
[K1,M1] = LuiBoundaryCondition(K_one,M_one,Gamma12,'u',1,Omega1_Modes{1},Omega2_Modes{1},InverseFunTwo);
% K_one = K1; M_one = M1;
K1(Gamma1,:) = [];
K1(:,Gamma1) = [];
M1(Gamma1,:) = [];
M1(:,Gamma1) = [];
K1 = sparse(K1);
M1 = sparse(M1);
[av1,O1] = eigs(K1,M1,1,'sm');
if abs(min(av1)) > abs(max(av1))
    av1 = -av1;
end
if abs(min(av2)) > abs(max(av2))
    av2 = -av2;
end

av2 = av2/max(av2);
av1 = av1/max(av1);
av1 = BoundariesPostProcess(av1,Gamma1);
Omega1_Modes = VisualizeModes(Omega1,av1,ID1);

% Plot Iteration
iter = iter+1;
if iter > max_iter - prod(graphs)
    tmp1 = strcat(str1,num2str(sqrt(O1)));
    tmp2 = strcat(str3,num2str(iter));
    figure(2)
    subplot(graphs(1),graphs(2),iter+prod(graphs)-max_iter)
    Omega1_Modes{1}.plot_geo('coarse',0,0, [0 1], [0 1]);
    hold on
    Omega2_Modes{1}.plot_geo('coarse',0,0, [0 1], [conf2 1]);
    title(strcat(tmp1,tmp2));
    alpha(0.9)
    xlabel('x [m]','FontWeight','bold')
    ylabel('y [m]','FontWeight','bold')
    zlabel('z [m]','FontWeight','bold')
end
%% Second Step
[K2, M2] = LuiBoundaryCondition(K_two,M_two,Gamma21,'u',0,Omega2_Modes{1},Omega1_Modes{1},InverseFunOne);
% K_two = K2; M_two = M2;
K2(Gamma2,:) = [];
K2(:,Gamma2) = [];
M2(Gamma2,:) = [];
M2(:,Gamma2) = [];
K2 = sparse(K2);
M2 = sparse(M2);
[av2,O2] = eigs(K2,M2,1,'sm');
if abs(min(av1)) > abs(max(av1))
    av1 = -av1;
end
if abs(min(av2)) > abs(max(av2))
    av2 = -av2;
end
av1 = av1/max(av1);
av2 = av2/max(av2);
av2 = BoundariesPostProcess(av2,Gamma2);
Omega2_Modes = VisualizeModes(Omega2,av2,ID2);

% Plot Iteration
iter = iter+1;
if iter > max_iter - prod(graphs)
    tmp1 = strcat(str2,num2str(sqrt(O2)));
    tmp2 = strcat(str3,num2str(iter));

    figure(2)
    subplot(graphs(1),graphs(2),iter+prod(graphs)-max_iter)
    Omega1_Modes{1}.plot_geo('coarse',0,0,[0 1],[0 conf1]);
    hold on
    Omega2_Modes{1}.plot_geo('coarse',0,0,[0 1], [0 1])
    title(strcat(tmp1,tmp2));
    xlabel('x [m]','FontWeight','bold')
    ylabel('y [m]','FontWeight','bold')
    zlabel('z [m]','FontWeight','bold')
    alpha(0.9)
end
end