clearvars
run startup.m
clc

p = 1; % Number of p refinements
h = 1; % Number of h-refinements
Results = cell(8,1);
while p < 8
    current_results = zeros(50,2);
    while h < 50
%% Geometry Input
    a = 1;
    b = 1.5;
    % Ruled Surface Definition
        P1 = [0 0 0 1];
        P2 = [0 a 0 1];
        line1 = geo_line(P1,P2); % First line
        P3 = [b 0 0 1];
        P4 = [b a 0 1];
        line2 = geo_line(P3,P4); % Second line
        Omega = geo_ruled(line1,line2); % The Ruled Surface defined by lines 1 and 2
    % Subdomains Definition
        overlap = 0.75; % Percentage of intersection area
        P5 = [b*(1+overlap)/2 0 0 1];
        P6 = [b*(1+overlap)/2 a 0 1];
        line3 = geo_line(P5,P6);
        P7 = [b*(1-overlap)/2 0 0 1];
        P8 = [b*(1-overlap)/2 a 0 1];
        line4 = geo_line(P7,P8);
        Omega1 = geo_ruled(line1,line3); % Subdomain 1
        Omega2 = geo_ruled(line4,line2); % Subdomain 2
%% Refinement
    Omega.DegreeElevate(p,1);
    Omega.DegreeElevate(p,2);
    Omega1.DegreeElevate(p,1);
    Omega1.DegreeElevate(p,2);
    Omega2.DegreeElevate(p,1);
    Omega2.DegreeElevate(p,2);

    interval = linspace(0,1,h+2);
    interval = setdiff(interval,[0, 1]);
    Omega.KnotRefine(interval,1);
    Omega.KnotRefine(.5,2);
    Omega1.KnotRefine(interval,1);
    Omega1.KnotRefine(.5,2);
    Omega2.KnotRefine(interval,1);
    Omega2.KnotRefine(.5,2);
%% Assembly
% Stiffness and Mass matrices
    NUMBER_OF_SOLUTION_DIMENTIONS = 1;
    [K, M] = PoissonProblemAssembly(Omega,NUMBER_OF_SOLUTION_DIMENTIONS);
    [K1,M1] = PoissonProblemAssembly(Omega1,NUMBER_OF_SOLUTION_DIMENTIONS);
    [K2,M2] = PoissonProblemAssembly(Omega2,NUMBER_OF_SOLUTION_DIMENTIONS);
% Boundary Definitions
    DomainBoundaries = GetBoundaries(Omega);
    SubdomainOneBoundaries = GetBoundaries(Omega1);
    SubdomainTwoBoundaries = GetBoundaries(Omega2);
     % Gamma 1 and 2 -> Boundaries defined by u = 0 and u = 1
        Gamma1 = DomainBoundaries{1};
        Gamma11 = SubdomainOneBoundaries{1};
        Gamma12 = SubdomainTwoBoundaries{1};
      
        Gamma2 = DomainBoundaries{2};
        Gamma21 = SubdomainOneBoundaries{2};
        Gamma22 = SubdomainTwoBoundaries{2};
    % Gammas 3 and 4 -> Boundaries defined by v = 0 and v = 1
        Gamma3 = DomainBoundaries{3};
        Gamma31 = SubdomainOneBoundaries{3};
        Gamma32 = SubdomainTwoBoundaries{3};
      
        Gamma4 = DomainBoundaries{4};
        Gamma41 = SubdomainOneBoundaries{4};
        Gamma42 = SubdomainTwoBoundaries{4};
    
%OBS: The overlapping boundaries in this case are:
%Omega1 -> Gamma41
%Omega2 -> Gamma32

% Boundary Restraints
    % Constrained Boundaries in this case are each boundary that is not
    % overlapping another domain
    % We have to run a unique() command so we don't constrain the same boundary
    % twice
    ConstrainedBoundaries = unique([Gamma1;Gamma2;Gamma3;Gamma4]);
    ConstrainedOneBoundaries = unique([Gamma11; Gamma21; Gamma31]);
    ConstrainedTwoBoundaries = unique([Gamma12; Gamma22; Gamma42]);

% Constraining the DOFs
    % First, let's "save" the original matrices
    KK = K; MM = M; 
    % Now we delete the constrained rows and columns from the matrices
    K(ConstrainedBoundaries,:) = [];
    K(:,ConstrainedBoundaries) = [];
    M(ConstrainedBoundaries,:) = [];
    M(:,ConstrainedBoundaries) = [];
    
    % Repeating for the SubDomains
    KK1 = K1; MM1 = M1;
    K1(ConstrainedOneBoundaries,:) = [];
    K1(:,ConstrainedOneBoundaries) = [];
    M1(ConstrainedOneBoundaries,:) = [];
    M1(:,ConstrainedOneBoundaries) = [];

    KK2 = K2; MM2 = M2;
    K2(ConstrainedTwoBoundaries,:) = [];
    K2(:,ConstrainedTwoBoundaries) = [];
    M2(ConstrainedTwoBoundaries,:) = [];
    M2(:,ConstrainedTwoBoundaries) = [];

% Solving the Eigenvalue Problems
    K = sparse(K);
    M = sparse(M);
    [EigenVector, EigenValues] = eigs(K,M,1,'sm');
    
    K1 = sparse(K1);
    M1 = sparse(M1);
    [EigenVector1, EigenValues1] = eigs(K1,M1,1,'sm');
    K1 = KK1;
    
    K2 = sparse(K2);
    M2 = sparse(M2);
    [EigenVector2, EigenValues2] = eigs(K2,M2,1,'sm');
    K2 = KK2;
%% Post-Processing
% We have to re-add the constrained DOFs in the right position as 0s in the
% EigenVector. Fortunately, we have made a function for that.
EigenVector = BoundariesPostProcess(EigenVector,ConstrainedBoundaries);
EigenVector1 = BoundariesPostProcess(EigenVector1,ConstrainedOneBoundaries);
EigenVector2 = BoundariesPostProcess(EigenVector2,ConstrainedTwoBoundaries);

% Visualization
    OmegaVisualizationCell = VisualizeModes(Omega,EigenVector);
    Omega1VisualizationCell = VisualizeModes(Omega1,EigenVector1);
    Omega2VisualizationCell = VisualizeModes(Omega2,EigenVector2);

%% Schwarz Algorithm Loops
% With the inital problem setted up, we now cycle through a series of Lui
% Boundary Conditions in Omega1 and Omega2 until we satisfy a convergence
% condition

% For that, we need to build inverse functions to communicate between the
% subdomains

OmegaOneInverse = @(x,y) [y, 2*x/(b*(1+overlap))];
tmp = b*(1-overlap)/2;
OmegaTwoInverse = @(x,y) [y, (x-tmp)/(b-tmp)];
i = 1;
nrows = 3;
ncols = 5;
% figure(1)
while i<nrows*ncols
    % Lui Boundary Condition for Omega 1
    K1 = LuiBC_2D(Omega1VisualizationCell{1},Omega2VisualizationCell{1}, ...
        KK1,Gamma41,4,OmegaTwoInverse);
        M1 = MM1;
        K1(ConstrainedOneBoundaries,:) = [];
        K1(:,ConstrainedOneBoundaries) = [];
        M1(ConstrainedOneBoundaries,:) = [];
        M1(:,ConstrainedOneBoundaries) = [];
        K1 = sparse(K1);
        M1 = sparse(M1);
        [EigenVector1, EigenValues1] = eigs(K1,M1,1,'sm');
        if abs(min(EigenVector1)) > abs(max(EigenVector1))
            EigenVector1 = -EigenVector1;
        end
        K1 = KK1;
    EigenVector1 = BoundariesPostProcess(EigenVector1,ConstrainedOneBoundaries);
    Omega1VisualizationCell = VisualizeModes(Omega1,EigenVector1);
    subplot(nrows,ncols,i)
%     Omega1VisualizationCell{1}.plot_geo
    i = i+1;
    % Lui Boundary Condition for Omega 2
    K2 = LuiBC_2D(Omega2VisualizationCell{1},Omega1VisualizationCell{1}, ...
        KK2,Gamma32,3,OmegaOneInverse);
        M2 = MM2;
        K2(ConstrainedTwoBoundaries,:) = [];
        K2(:,ConstrainedTwoBoundaries) = [];
        M2(ConstrainedTwoBoundaries,:) = [];
        M2(:,ConstrainedTwoBoundaries) = [];
        K2 = sparse(K2);
        M2 = sparse(M2);
        [EigenVector2, EigenValues2] = eigs(K2,M2,1,'sm');
        if abs(min(EigenVector2)) > abs(max(EigenVector2))
            EigenVector2 = -EigenVector2;
        end
        K2 = KK2;
    EigenVector2 = BoundariesPostProcess(EigenVector2,ConstrainedTwoBoundaries);
    Omega2VisualizationCell = VisualizeModes(Omega2,EigenVector2);
    subplot(nrows,ncols,i)
%     Omega2VisualizationCell{1}.plot_geo
end
    current_results(h,1) = EigenValues1;
    current_results(h,2) = EigenValues2;
    current_results(h,3) = EigenValues;
    h = h+1;
    end
    Results{p} = current_results;
p = p+1;
end
