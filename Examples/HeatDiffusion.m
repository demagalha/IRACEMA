clearvars
run startup.m

%% Geometry Definition
point1 = [0 0 0 1]; % x y z and weight
point2 = [0.6 0 0 1];

Line1 = geo_line(point1, point2);

point3 = [0 1 0 1];
point4 = [0.6 1 0 1];
Line2 = geo_line(point3,point4); %Line2 is the same as Line1, but 1m up in y
Omega = geo_ruled(Line1,Line2);


% k-refinement
Omega.DegreeElevate(1,1); % Elevate u-direction by 2 degrees
Omega.DegreeElevate(1,2); % Elevate v-direction by 2 degrees

uKnots = linspace(0,1,15);
uKnots = setdiff(uKnots, [0 1]); % Making it an open interval

vKnots = linspace(0,1,11);
vKnots = setdiff(vKnots,[0 1]);

Omega.KnotRefine(uKnots,1);
Omega.KnotRefine(vKnots,2);

%% Processing
% Assembly of the Model
SOLUTION_DIMENSIONS = 1; % Scalar problem
[K, ~] = PoissonProblemAssembly(Omega,SOLUTION_DIMENSIONS);
alpha = 52; % Heat Diffusivity Coefficient
K = alpha*K;
F = zeros(length(K),1);
% Boundary Conditions
Boundaries = GetBoundaries(Omega);

% By definition, GetBoundaries returns the boundaries in the order:
% Gamma 1 -> u = 0
% Gamma 2 -> u = 1
% Gamma 3 -> v = 0
% Gamma 4 -> v = 1
%  In this plate model, 
% 
%     _______Gamma4__________
% G  |                       | G
% a  |                       | a
% m  |                       | m
% m  |                       | m
% a  |_______________________| a
% 1          Gamma 3           2
% 
% y,v
% Î -> x,u

Gamma1Basis = Boundaries{1,1};
Gamma1Elements = Boundaries{1,2};

Gamma2Basis = Boundaries{2,1};
Gamma2Elements = Boundaries{2,2};

Gamma3Basis = Boundaries{3,1};
Gamma3Elements = Boundaries{3,2};

Gamma4Basis = Boundaries{4,1};
Gamma4Elements = Boundaries{4,2};

% Neumann Boundary Conditions
% We have Homogeneous Neumann Boundary conditions at Gamma1
% Since this is the natural kind of BC, they don't need to be enforced.
% But, for didatic purposes, we will run the script
h = 0;
gamma1_values = zeros(length(Gamma1Elements),1);
gamma1_values(:) = h;

% F = ApplyNeumannBCs(Omega,F,Gamma1Elements,1,gamma1_values);

% Robin Boundary Conditions
% We have Gamma2 and Gamma4 convecting to 0C with Heat Transfer Coefficient
% of 750W/(mK) 

h = 750;
T_inf = 273.15; % Zero Celsius, in Kelvin
ROBIN = h*(T_inf)/alpha;
BETA = h/alpha;

gamma2_values = zeros(length(Gamma2Elements),2);
gamma4_values = zeros(length(Gamma4Elements),2);

gamma2_values(:,1) = ROBIN;
gamma2_values(:,2) = BETA;
gamma4_values(:,1) = ROBIN;
gamma4_values(:,2) = BETA;

[K,F] = ApplyRobinBCs(Omega,K,F,Gamma2Basis,2,gamma2_values);
[K,F] = ApplyRobinBCs(Omega,K,F,Gamma4Basis,4,gamma4_values);



% Dirichlet Boundary Conditions
% Our method to apply Dirichlet BCs also gives the solution of the problem,
% so we will have to bundle up every BC together in a single array.
% Fortunately for us, only Gamma3 has DirichletBCs and everything is neat
% and tidy already.
% Also, since we will apply the BCs directly to the solution space, we only
% need the Basis address and the values
lift = 373.15; % 100 Celsius, in Kelvin
gamma3_values = zeros(length(Gamma3Basis),1);
gamma3_values(:) = lift;

[d, F] = ApplyDirichletBCs(K,F,Gamma3Basis,gamma3_values);

%% Post-Processing
[~, element_local_mapping, ~] = GetConnectivityArrays(Omega);
[ID, ~] = BuildGlobalLocalMatrices(element_local_mapping,SOLUTION_DIMENSIONS);
PlotDisplacement(d,ID,Omega);
caxis([273.15 373.15])
view([-0.96 90.00])
MM = VisualizeModes(Omega,d,ID);
Convergence = MM{1}.eval_point(1,0.2);
str = 'Temperature at x = 0.6m, y = 0.2m:';
str2 = num2str(Convergence.z - 273.15, '%.2f');
a = strcat({str},{' '},{str2},{'C'});
display(a{1});