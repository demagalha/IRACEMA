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

% Boundary Conditions
Gamma3 = GetBoundaryConditionArray(Omega,2,0,373.15);
Gamma1Elements = GetBoundaryElements(Omega,1,0,0);
h = -750;
T_inf = 273.15; % Zero Celsius, in Kelvin
ROBIN = -h*T_inf/alpha;
BETA = -h/alpha;
Gamma2Elements = GetBoundaryElements(Omega,1,1,[ROBIN,BETA]);
Gamma4Elements = GetBoundaryElements(Omega,2,1,[ROBIN,BETA]);

F = zeros(length(K),1);
RobinElements = [Gamma2Elements; Gamma4Elements];
NeumannElements = Gamma1Elements;
[K,F] = RobinBC(Omega,K,F,RobinElements);
[~,F] = RobinBC(Omega,K,F,NeumannElements);

[d,~] = DirichletBC(K,F,Gamma3);
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