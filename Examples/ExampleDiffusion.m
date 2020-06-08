close all
clearvars
clc

%% Model parameters

P1 = [0 0 0 1]; %Control Point 1, x y z weight
P2 = [0.6 0 0 1]; %Control Point 2, x y z weight
pu = 1; % Polynomial Degree
U = [0 0 1 1]; % Knot Vector

c1 = Geometry('curve',pu,U,{P1,P2}); % A straight line from y=0 to y=1

P3 = [0 1 0 1];
P4 = [0.6 1 0 1];

c2 = Geometry('curve',pu,U,{P3,P4}); % A straight line from (x,y) = (0,1) to (1,1)

Model = geo_ruled(c1,c2); % A ruled surface from curve 1 to curve 2

clearvars -except Model pu U

%% Refinement
p = 1; % Number of p refinements
    Model.DegreeElevate(p,1); % Elevate the degree of direction 1 p times
    Model.DegreeElevate(p,2); % Elevate the degree of direction 2 p times

h = 20; %Number of h refinements in each parametric direction
    interval = linspace(0,1,h+2);
    interval = interval(2:end-1);
    
    Model.KnotRefine(interval,1);% Add h equally spaced knots in direction 1
    Model.KnotRefine(interval,2);% Add h equally spaced knots in direction 2

%% Assembly
% Membrane Eigenvalue problem Stiffness = Heat equation stiffness
[K, ~, ID] = MembraneAssemble(Model);
alpha = 52; % Heat Diffusivity coefficient
K = alpha*K;
F = zeros(length(K),1);
%% BoundaryConditions
% Domain definition
%        Gamma2  
%u=1_________________________ 
%G |                         |
%a |                         |
%m |                         |
%m | v = 0               v=1 |Gamma4
%a |                         |
%3 |                         |
%  |_________________________|
% u=0       Gamma1

% Data for the problem:
% Boundary Conditions are 
% 373.15K at y = 0m -> Dirichlet Imposing
Gamma1 = GetBoundaryConditionArray(Model,2,0,373.15);
DirichletStrong = true;
DirichletBoundaryArray = Gamma1;

% Insulated at x = 0m -> Neumann Condition
Gamma2 = GetBoundaryElements(Model,1,0,[0, 0]);
NeumannEnforcement = true;
NeumannBoundaryElements = Gamma2;

% Convecting at 750W/m with Ta = 273.15 at x = 0.6m & at y = 1m -> Robin
% From Convection BC: -k dTdt = h(T_inf - T)
% So robin = -hT_inf/k
% BETA = -h/k
h = -750;
T_inf = 273.15; % Kelvin
robin = -h*T_inf/alpha;
BETA = -h/alpha;
Gamma3 = GetBoundaryElements(Model,1,1,[robin,BETA]);
Gamma4 = GetBoundaryElements(Model,2,1,[robin,BETA]);

RobinEnforcement = true;
RobinBoundaryElements = [Gamma3; Gamma4];

% OBS: We use
% GetBoundaryConditionArray(Model,direction,boundary_value,lift) for strong
% imposition of BCs and
% GetBoundaryElements(Model,direction,boundary_value,lift) for weak
% imposition of BCs.
if RobinEnforcement
    [K, F] = RobinBC(Model,K,F,RobinBoundaryElements);
end

% Neumann Enforcement of Boundary Conditions
if NeumannEnforcement
   [~, F] = RobinBC(Model,K,F,NeumannBoundaryElements);
end
% Enforcement of Robin Boundary Condition


%% Solution of the Model
% Strong Enforcement of Dirichlet Boundary Conditions
if DirichletStrong
    d = zeros(size(F));
    BoundaryDOFS = DirichletBoundaryArray(:,1);
    FreeDOFS = setdiff(1:length(d),BoundaryDOFS);
        
    d(BoundaryDOFS) = DirichletBoundaryArray(:,2);
    F(FreeDOFS) = F(FreeDOFS) - K(FreeDOFS,BoundaryDOFS)*DirichletBoundaryArray(:,2);
    d(FreeDOFS) = K(FreeDOFS,FreeDOFS)\F(FreeDOFS);
else
    d = K\F;
end

%% Plot
PlotDisplacement(d,ID,Model);
MM = VisualizeModes(Model,d,ID);
Convergence = MM{1}.eval_point(1,0.2)
